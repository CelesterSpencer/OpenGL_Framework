#include "PerfectSurfaceDetection.h"

#include "CPPImplementation.h"

#include "OrbitCamera.h"
#include "ShaderTools/Renderer.h"
#include "Molecule/MDtrajLoader/MdTraj/MdTrajWrapper.h"
#include "Molecule/MDtrajLoader/Data/Protein.h"

#include <glm/gtx/component_wise.hpp>

#include <iostream>

// ### Class implementation ###

PerfectSurfaceDetection::PerfectSurfaceDetection()
{
    // Setup members
    mRotateCamera = false;
    mSurfaceAtomCount = 0;

    // Create window
    mpWindow = generateWindow();

    // Register keyboard callback
    std::function<void(int, int, int, int)> kC = [&](int k, int s, int a, int m)
    {
        this->keyCallback(k, s, a, m);
    };
    setKeyCallback(mpWindow, kC);

    // Register mouse button callback
    std::function<void(int, int, int)> kB = [&](int b, int a, int m)
    {
        this->mouseButtonCallback(b, a, m);
    };
    setMouseButtonCallback(mpWindow, kB);

    // Register scroll callback
    std::function<void(double, double)> kS = [&](double x, double y)
    {
        this->scrollCallback(x,y);
    };
    setScrollCallback(mpWindow, kS);

    // Clear color
    glClearColor(0.f, 0.f, 0.f, 1.f);

    // Path to protein molecule
    std::vector<std::string> paths;
    paths.push_back(std::string(RESOURCES_PATH) + "/molecules/PDB/1crn.pdb");

    // Load protein
    MdTrajWrapper mdwrap;
    mupProtein = mdwrap.load(paths);

    // Atom count
    mAtomCount = (int) mupProtein->getAtoms()->size();

    // Output atom count
    std::cout << "Atom count: " << mAtomCount << std::endl;

    // Get min/max extent of protein
    mupProtein->minMax(); // first, one has to calculate min and max value of protein
    glm::vec3 proteinMinExtent = mupProtein->getMin();
    glm::vec3 proteinMaxExtent = mupProtein->getMax();

    // Test protein extent
    std::cout
        << "Min extent of protein: "
        << proteinMinExtent.x << ", "
        << proteinMinExtent.y << ", "
        << proteinMinExtent.z << std::endl;
    std::cout
        << "Max extent of protein: "
        << proteinMaxExtent.x << ", "
        << proteinMaxExtent.y << ", "
        << proteinMaxExtent.z << std::endl;

    /*
    // Test atom radii
    std::vector<Atom*>* pAtoms = mupProtein->getAtoms();
    for(int i = 0; i < (int)pAtoms->size(); i++)
    {
        std::string element = pAtoms->at(i)->getElement();
        std::cout << "Atom: " <<  element << " Radius: " << atomLUT.vdW_radii_picometer.at(element) << std::endl;
    }
    */

    // Create camera
    glm::vec3 cameraCenter = (proteinMinExtent + proteinMaxExtent) / 2.f;
    float cameraRadius = glm::compMax(proteinMaxExtent - cameraCenter);
    mupCamera = std::unique_ptr<OrbitCamera>(new OrbitCamera(cameraCenter, 90.f, 90.f, cameraRadius, cameraRadius / 2.f, 3.f * cameraRadius));

    // Create query to measure execution time
    glGenQueries(1, &mQuery);

     // Variable to measure elapsed time
    GLuint timeElapsed = 0;

    // ### Set up compute shader ###

    // Start query for time measurement
    glBeginQuery(GL_TIME_ELAPSED, mQuery);

    // # Prepare atoms input (position + radius)

    // Vector which is used as data for SSBO or CPP implementation
    for(Atom const * pAtom : *(mupProtein->getAtoms()))
    {
        // Push back all atoms
        mAtomStructs.push_back(
            AtomStruct(
                pAtom->getPosition(),
                mAtomLUT.vdW_radii_picometer.at(
                    pAtom->getElement())));
    }

    // Fill into ssbo
    glGenBuffers(1, &mAtomsSSBO);
    glBindBuffer(GL_SHADER_STORAGE_BUFFER, mAtomsSSBO);
    glBufferData(GL_SHADER_STORAGE_BUFFER, sizeof(AtomStruct) * mAtomStructs.size(), mAtomStructs.data(), GL_STATIC_DRAW);
    glBindBuffer(GL_SHADER_STORAGE_BUFFER, 0);

    // # Prepare uint image to write indices of surface atoms. Lets call it image list in shader for easier understanding

    // Buffer
    glGenBuffers(1, &mSurfaceAtomBuffer);
    glBindBuffer(GL_TEXTURE_BUFFER, mSurfaceAtomBuffer);
    glBufferData(GL_TEXTURE_BUFFER, sizeof(GLuint) * mAtomCount, 0, GL_STATIC_DRAW);
    glBindBuffer(GL_TEXTURE_BUFFER, 0);

    // Texture (which will be bound as image)
    glGenTextures(1, &mSurfaceAtomTexture);
    glBindTexture(GL_TEXTURE_BUFFER, mSurfaceAtomTexture);
    glTexBuffer(GL_TEXTURE_BUFFER, GL_R32UI, mSurfaceAtomBuffer);
    glBindTexture(GL_TEXTURE_BUFFER, 0);

    // Print time for data transfer
    glEndQuery(GL_TIME_ELAPSED);
    glGetQueryObjectuiv(mQuery, GL_QUERY_RESULT, &timeElapsed);
    std::cout << "Time for data transfer on GPU: " << std::to_string(timeElapsed) << "ns" << std::endl;

    // Run implementation
    #ifdef USE_GLSL_IMPLEMENTATION
        runGLSLImplementation();
    #else
        runCPPImplementation();
    #endif

    // Output count
    std::cout << "Surface atom count: " << mSurfaceAtomCount << std::endl;
}

PerfectSurfaceDetection::~PerfectSurfaceDetection()
{
    // TODO: Delete OpenGL stuff
    // - ssbo
    // -

    // Delete query object
    glDeleteQueries(1, &mQuery);
}

void PerfectSurfaceDetection::renderLoop()
{
    // Point size for rendering
    glPointSize(15.f);

    // Cursor
    float prevCursorX, prevCursorY = 0;

    // Prepare shader programs for rendering
    ShaderProgram proteinPointProgram = ShaderProgram("/PerfectSurfaceDetection/proteinPoint.vert", "/PerfectSurfaceDetection/point.frag");
    ShaderProgram surfacePointProgram = ShaderProgram("/PerfectSurfaceDetection/surfacePoint.vert", "/PerfectSurfaceDetection/point.frag");

    // Bind SSBO with atoms
    glBindBufferBase(GL_SHADER_STORAGE_BUFFER, 0, mAtomsSSBO);

    // Bind image where indices of surface atoms are read from
    glBindImageTexture(1,
       mSurfaceAtomTexture,
       0,
       GL_TRUE,
       0,
       GL_READ_ONLY,
       GL_R32UI);

    // Projection matrix (hardcoded viewport size)
    proteinPointProgram.use();
    proteinPointProgram.update("projection", glm::perspective(glm::radians(45.f), (GLfloat)1280 / (GLfloat)720, 0.1f, 1000.f));
    surfacePointProgram.use();
    surfacePointProgram.update("projection", glm::perspective(glm::radians(45.f), (GLfloat)1280 / (GLfloat)720, 0.1f, 1000.f));

    // Call render function of Rendering.h with lambda function
    render(mpWindow, [&] (float deltaTime)
    {
        // Clear buffers
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

        // Calculate cursor movement
        double cursorX, cursorY;
        glfwGetCursorPos(mpWindow, &cursorX, &cursorY);
        GLfloat cursorDeltaX = (float)cursorX - prevCursorX;
        GLfloat cursorDeltaY = (float)cursorY - prevCursorY;
        prevCursorX = cursorX;
        prevCursorY = cursorY;

        // Orbit camera
        if(mRotateCamera)
        {
            mupCamera->setAlpha(mupCamera->getAlpha() + 0.25f * cursorDeltaX);
            mupCamera->setBeta(mupCamera->getBeta() - 0.25f * cursorDeltaY);
        }
        mupCamera->update();

        // Draw complete protein
        proteinPointProgram.use();
        proteinPointProgram.update("view", mupCamera->getViewMatrix());
        glDrawArrays(GL_POINTS, 0, mAtomCount);

        // Draw surface atoms
        surfacePointProgram.use();
        surfacePointProgram.update("view", mupCamera->getViewMatrix());
        glDrawArrays(GL_POINTS, 0, mSurfaceAtomCount);
    });
}

void PerfectSurfaceDetection::keyCallback(int key, int scancode, int action, int mods)
{
    if (action == GLFW_PRESS)
    {
        switch (key)
        {
        case GLFW_KEY_ESCAPE: { glfwSetWindowShouldClose(mpWindow, GL_TRUE); break; } // Does not work, but method gets called?!
        }
    }
}

void PerfectSurfaceDetection::mouseButtonCallback(int button, int action, int mods)
{
    if (button == GLFW_MOUSE_BUTTON_LEFT && action == GLFW_PRESS)
    {
        mRotateCamera = true;
    }
    else if (button == GLFW_MOUSE_BUTTON_LEFT && action == GLFW_RELEASE)
    {
        mRotateCamera = false;
    }
}

void PerfectSurfaceDetection::scrollCallback(double xoffset, double yoffset)
{
    mupCamera->setRadius(mupCamera->getRadius() - 0.5f * (float)yoffset);
}

GLuint PerfectSurfaceDetection::readAtomicCounter(GLuint atomicCounter) const
{
    // Read atomic counter
    glBindBuffer(GL_ATOMIC_COUNTER_BUFFER, atomicCounter);

    GLuint *mapping = (GLuint*)glMapBufferRange(GL_ATOMIC_COUNTER_BUFFER,
                                                0,
                                                sizeof(GLuint),
                                                GL_MAP_READ_BIT);

    glUnmapBuffer(GL_ATOMIC_COUNTER_BUFFER);
    glBindBuffer(GL_ATOMIC_COUNTER_BUFFER, 0);

    return mapping[0];
}

void PerfectSurfaceDetection::resetAtomicCounter(GLuint atomicCounter) const
{
    glBindBuffer(GL_ATOMIC_COUNTER_BUFFER, atomicCounter);

    // Map the buffer
    GLuint* mapping = (GLuint*)glMapBufferRange(GL_ATOMIC_COUNTER_BUFFER,
                                                0 ,
                                                sizeof(GLuint),
                                                GL_MAP_WRITE_BIT | GL_MAP_INVALIDATE_BUFFER_BIT | GL_MAP_UNSYNCHRONIZED_BIT);
    // Set memory to new value
    memset(mapping, 0, sizeof(GLuint));

    // Unmap the buffer
    glUnmapBuffer(GL_ATOMIC_COUNTER_BUFFER);
    glBindBuffer(GL_ATOMIC_COUNTER_BUFFER, 0);
}

void PerfectSurfaceDetection::runCPPImplementation()
{
    // # Input vector with atoms
}



void PerfectSurfaceDetection::runGLSLImplementation()
{
     // # Compile shader
    ShaderProgram computeProgram(GL_COMPUTE_SHADER, "/PerfectSurfaceDetection/surface.comp");

    // # Prepare atomic counter for writing results to unique position in image
    GLuint atomicCounter;
    glGenBuffers(1, &atomicCounter);
    glBindBuffer(GL_ATOMIC_COUNTER_BUFFER, atomicCounter);
    glBufferData(GL_ATOMIC_COUNTER_BUFFER, sizeof(GLuint), NULL, GL_DYNAMIC_DRAW);
    resetAtomicCounter(atomicCounter);

    // # Execute compute shader to determine surface atoms

    // Variable to measure elapsed time
    GLuint timeElapsed = 0;

    // Start query for time measurement
    glBeginQuery(GL_TIME_ELAPSED, mQuery);

    // Use compute shader program
    computeProgram.use();

    // Bind SSBO with atoms
    glBindBufferBase(GL_SHADER_STORAGE_BUFFER, 0, mAtomsSSBO);

    // Tell shader about count of atoms
    computeProgram.update("atomCount", mAtomCount);

    // Probe radius
    computeProgram.update("probeRadius", 140.f);

    // Bind atomic counter
    glBindBufferBase(GL_ATOMIC_COUNTER_BUFFER, 1, atomicCounter);

    // Bind image where indices of surface atoms are written to
    glBindImageTexture(2,
                       mSurfaceAtomTexture,
                       0,
                       GL_TRUE,
                       0,
                       GL_WRITE_ONLY,
                       GL_R32UI);

    // Print time for setup
    glEndQuery(GL_TIME_ELAPSED);
    glGetQueryObjectuiv(mQuery, GL_QUERY_RESULT, &timeElapsed);
    std::cout << "Time for setup on GPU: " << std::to_string(timeElapsed) << "ns" << std::endl;

    // Start query for time measurement
    glBeginQuery(GL_TIME_ELAPSED, mQuery);

    // Dispatch
    glDispatchCompute((mAtomCount / 16) + 1, 1, 1);
    glMemoryBarrier(GL_ALL_BARRIER_BITS);

    // Print time for execution
    glEndQuery(GL_TIME_ELAPSED);
    glGetQueryObjectuiv(mQuery, GL_QUERY_RESULT, &timeElapsed);
    std::cout << "Time for execution on GPU: " << std::to_string(timeElapsed) << "ns (= " << std::to_string(timeElapsed / 1000000.f) << "ms)" << std::endl;

    // Fetch count
    mSurfaceAtomCount = readAtomicCounter(atomicCounter);

    // TODO: Delete atomic counter
}

// ### Main function ###

int main()
{
    PerfectSurfaceDetection detection;
    detection.renderLoop();
    return 0;
}

// ### Snippets ###

/*
    // Read back SSBO
    glBindBuffer(GL_SHADER_STORAGE_BUFFER, atomsSSBO);
    AtomStruct *ptr;
    ptr = (AtomStruct *) glMapBuffer(GL_SHADER_STORAGE_BUFFER, GL_READ_ONLY);
    for(int i = 0; i < atomCount; i++)
    {
        std::cout << i << ". " << ptr[i].radius << std::endl;
    }
    glUnmapBuffer(GL_SHADER_STORAGE_BUFFER);

*/
