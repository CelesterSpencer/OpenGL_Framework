#include "Texture.h"

Texture::Texture(std::string path) {

    textureHandle = load(path);
}

Texture::Texture() {
}

GLuint Texture::getInternalFormat() const
{
    return internalFormat;
}

GLuint Texture::getFormat() const
{
    return format;
}

GLuint Texture::getType() const
{
    return type;
}

GLuint Texture::getTarget() const
{
    return target;
}

bool Texture::getIsImageTex() const
{
    return isImageTex;
}

GLuint Texture::getIsLayered() const
{
    return isLayered;
}

int Texture::getX() const
{
    return x;
}

int Texture::getY() const
{
    return y;
}

int Texture::getZ() const
{
    return z;
}

void Texture::reset()
{
    glBindTexture(target, textureHandle);
    //glBindImageTexture(0, textureHandle, 0, GL_TRUE, 0, GL_READ_WRITE, internalFormat);
    glClearTexImage( textureHandle, 0, format, type, pixels);
}

Texture::~Texture() {

}

Texture::Texture(int w, int h) {
    this->w = w;
    this->h = h;
    textureHandle = this->genTexture(w, h);
}

Texture::Texture(GLuint internalFormat, GLuint format, GLuint type):
    type(type),
    internalFormat(internalFormat),
    format(format),
    isImageTex(false)
{

}

GLuint Texture::getHandle() {
    return textureHandle;
}

GLuint Texture::load(std::string path) {
    if (!devILInitialized) {
        ilInit();
        iluInit();
        ilutRenderer(ILUT_OPENGL);
        devILInitialized = true;
    }

    ILuint iid;
    ilGenImages(1, &iid);
    ilBindImage(iid);
    ilLoadImage(path.c_str());

    int w = ilGetInteger(IL_IMAGE_WIDTH);
    int h = ilGetInteger(IL_IMAGE_HEIGHT);
    const int byteCount = sizeof(float) * 3 * w * h;
    unsigned char* pixels = new unsigned char[byteCount];
    ilCopyPixels(0, 0, 0, w, h, 1, IL_RGB, IL_FLOAT, pixels);

    GLuint t;
    glGenTextures(1, &t);
    glActiveTexture(GL_TEXTURE0);
    glBindTexture(GL_TEXTURE_2D, t);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
    glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA16F, w, h, 0, GL_RGB, GL_FLOAT, pixels);
    glBindImageTexture(0, t, 0, GL_FALSE, 0, GL_READ_ONLY, GL_RGBA16F);

    return t;
}

void Texture::clear() {
    //glActiveTexture(GL_TEXTURE0);
    glBindTexture(GL_TEXTURE_2D, textureHandle);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
    glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA16F, w, h, 0, GL_RGBA, GL_FLOAT, pixels);
    glBindImageTexture(0, textureHandle, 0, GL_FALSE, 0, GL_READ_WRITE, GL_RGBA16F);
}

GLuint Texture::genTexture(int w, int h) {

    byteCount = sizeof(float) * 4 * w * h;
    pixels = new unsigned char[byteCount];

    float f = -1000.0f;

    unsigned char const * p = reinterpret_cast<unsigned char const *>(&f);

    for (int i = 0; i < byteCount; i+=4)
    {
        pixels[i] = p[0];
        pixels[i+1] = p[1];
        pixels[i+2] = p[2];
        pixels[i+3] = p[3];
    }

    GLuint t;
    glGenTextures(1, &t);
    glActiveTexture(GL_TEXTURE0);
    glBindTexture(GL_TEXTURE_2D, t);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
    glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA16F, w, h, 0, GL_RGBA, GL_FLOAT, pixels);
    glBindImageTexture(0, t, 0, GL_FALSE, 0, GL_READ_WRITE, GL_RGBA16F);

    return t;
}

GLuint Texture::gen2DTexture(int w, int h)
{
    target = GL_TEXTURE_2D;
    isImageTex = true;
    isLayered = GL_FALSE;

    byteCount = sizeof(unsigned int) * w*h;
    pixels = new unsigned char[byteCount];

    unsigned int f = 0;

    unsigned char const * p = reinterpret_cast<unsigned char const *>(&f);

    for (int i = 0; i < byteCount; i+=4)
    {
        pixels[i] = p[0];
        pixels[i+1] = p[1];
        pixels[i+2] = p[2];
        pixels[i+3] = p[3];
    }

    GLuint t;
    glGenTextures(1, &t);
    glActiveTexture(GL_TEXTURE0);
    glBindTexture(target, t);
    glTexParameteri(target, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
    glTexParameteri(target, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
    glTexImage2D(target, 0, internalFormat, w, h, 0, format, type, pixels);
    glBindImageTexture(0, t, 0, GL_FALSE, 0, GL_READ_WRITE, internalFormat);

    textureHandle = t;

    return t;
}

GLuint Texture::genUimageBuffer(int size)
{
    target = GL_TEXTURE_1D;
    isImageTex = true;
    isLayered = GL_FALSE;

    byteCount = sizeof(unsigned int) * size;
    pixels = new unsigned char[byteCount];

    unsigned int f = 0;

    unsigned char const * p = reinterpret_cast<unsigned char const *>(&f);

    for (int i = 0; i < byteCount; i+=4)
    {
        pixels[i] = p[0];
        pixels[i+1] = p[1];
        pixels[i+2] = p[2];
        pixels[i+3] = p[3];
    }

    GLuint t;
    glGenTextures(1, &t);
    glActiveTexture(GL_TEXTURE0);
    glBindTexture(target, t);
    glTexParameteri(target, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
    glTexParameteri(target, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
    glTexImage1D(target, 0, internalFormat, size, 0, format, type, pixels);
    glBindImageTexture(0, t, 0, GL_FALSE, 0, GL_READ_WRITE, internalFormat);

    textureHandle = t;

    return t;
}

GLuint Texture::gen3DTexture(int x, int y, int z)
{
    target = GL_TEXTURE_3D;
    isImageTex = true;
    isLayered = GL_TRUE;
    this->x =x;
    this->y =y;
    this->z =z;

    byteCount = sizeof(unsigned int) * x*y*z;
    pixels = new unsigned char[byteCount];

    float f = 0;

    unsigned char const * p = reinterpret_cast<unsigned char const *>(&f);

    for (int i = 0; i < byteCount; i+=4)
    {
        pixels[i] = p[0];
        pixels[i+1] = p[1];
        pixels[i+2] = p[2];
        pixels[i+3] = p[3];
    }

    GLuint t;
    glGenTextures(1, &t);
    glActiveTexture(GL_TEXTURE0);
    glBindTexture(target, t);
    glTexParameteri(target, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
    glTexParameteri(target, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
    glTexParameteri(target, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
    glTexParameteri(target, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
    glTexImage3D(target, 0, internalFormat, x, y, z, 0, format, type, 0);
    glBindImageTexture(0, t, 0, GL_TRUE, 0, GL_READ_WRITE, internalFormat);
    glClearTexImage( t, 0, format, type, pixels);

    textureHandle = t;

    return t;
}
