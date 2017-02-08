//============================================================================
// Distributed under the MIT License. Author: Adrian Derstroff
//============================================================================

// external includes
#include <GL/glew.h>
#include <GLFW/glfw3.h>
#include <imgui/imgui.h>
#include <imgui/examples/opengl3_example/imgui_impl_glfw_gl3.h>
#include <sstream>
#include <iomanip>
#include <Utils/PickingTexture.h>

// framework includes
#include "ShaderTools/Renderer.h"
#include "ShaderTools/ShaderProgram.h"
#include "Utils/Logger.h"
#include "NeighborSearch/NeighborhoodSearch.h"
#include "NeighborSearch/GPUHandler.h"

// local includes
#include "ProteinLoader.h"
#include "Utils/OrbitCamera.h"
#include "AmberForceFieldParameterLoader.h"
#include "Utils/PickingTexture.h"






/*
 * DEFINES
 */
#define WIDTH 1280
#define HEIGHT 720






/*
 * VARIABLES
 */
// window
GLFWwindow* mp_Window;

// interaction
std::unique_ptr<OrbitCamera>    mp_camera;
glm::vec2                       m_CameraDeltaMovement;
float                           m_CameraSmoothTime;
bool                            m_rotateCamera = false;
PickingTexture                  m_pickingTexture;
ShaderProgram                   m_idPickingShader;
int                             m_interactionMode = 0;
glm::vec3                       m_centerOfGravityCameraSpace;
glm::vec3                       m_reprojectedMousePositionCameraSpace;
bool                            m_resetStartMousePosition = true;
glm::vec2                       m_startMousePosition;
glm::vec2                       m_currentMousePosition;
float                           m_translationDelta = 0.0;
float                           m_rotationAngle = 0.0;
ShaderProgram                   m_linesShader;
glm::mat4                       m_tempModelMatrix;

// protein
ProteinLoader   m_proteinLoader;
int             m_selectedAtom = 0;
int             m_selectedProtein = 0;
float           m_proteinMoveSpeed = 2.f;
GLuint m_atomsSSBO;
GLuint m_indexCubeSSBO;
GLuint m_valuesVectorSSBO;
GLuint m_bondNeighborsSSBO;
ShaderProgram m_calcChargeProgram;
uint m_numberOfAtomSymbols;
GLuint m_modelMatricesBuffer;
std::vector<glm::vec4> m_proteinColors;
GLuint m_proteinColorsBuffer;

// neighborhood search
ShaderProgram m_extractAtomPositionsShader;
ShaderProgram m_drawSearchRadiusShader;
NeighborhoodSearch m_search;
GLuint m_positionsSSBO;
GLuint m_searchResultsSSBO;
glm::vec3 m_gridResolution = glm::vec3(10.0, 10.0, 10.0);
float m_searchRadius = 10.f;
GLuint m_gridVBO;
GLuint m_gridVAO;
bool m_isGridVisible = false;

// amber force calculation
bool useSimplifiedCalculation = false;

// rendering
glm::vec3   m_lightDirection;
ShaderProgram m_impostorProgram;
GLuint m_pointsVBO;
GLuint m_pointsVAO;
int    m_numVBOEntries;
bool m_drawSelectedProtein = false;






/*
 * forward declarations
 */
void setup();
void compileShaderPrograms();

void keyCallback(int key, int scancode, int action, int mods);
void mouseButtonCallback(int button, int action, int mods);
void scrollCallback(double xoffset, double yoffset);
void mousePositionCallback(double xpos, double ypos);

void run();

void updateGUI();






/*
 * INITIALIZATION
 */
void setup()
{
    /*
     * init window
     */
    mp_Window = generateWindow("Test", WIDTH, HEIGHT);

    /*
     * init imgui
     */
    ImGui_ImplGlfwGL3_Init(mp_Window, true);
    ImGuiIO& io = ImGui::GetIO();

    /*
     * setup opengl
     */
    glEnable(GL_DEPTH_TEST);
    glDepthFunc(GL_LESS);
    glEnable(GL_CULL_FACE);

    /*
     * register callbacks
     */
    // key press callback
    std::function<void(int, int, int, int)> kC = [&](int k, int s, int a, int m)
    {
        // Check whether ImGui is handling this
        ImGuiIO& io = ImGui::GetIO();
        if(io.WantCaptureKeyboard)
        {
            return;
        }
        keyCallback(k, s, a, m);
    };
    setKeyCallback(mp_Window, kC);

    // mouse button callback
    std::function<void(int, int, int)> kB = [&](int b, int a, int m)
    {
        // Check whether ImGui is handling this
        ImGuiIO& io = ImGui::GetIO();
        if(io.WantCaptureMouse)
        {
            return;
        }
        mouseButtonCallback(b, a, m);
    };
    setMouseButtonCallback(mp_Window, kB);

    // mouse scroll callback
    std::function<void(double, double)> kS = [&](double x, double y)
    {
        // Check whether ImGui is handling this
        ImGuiIO& io = ImGui::GetIO();
        if(io.WantCaptureMouse)
        {
            return;
        }
        scrollCallback(x,y);
    };
    setScrollCallback(mp_Window, kS);

    // cursor position callback
    std::function<void(double, double)> kM = [&](double x, double y)
    {
        // Check whether ImGui is handling this
        ImGuiIO& io = ImGui::GetIO();
        if(io.WantCaptureMouse)
        {
            return;
        }
        mousePositionCallback(x,y);
    };
    setCursorPosCallback(mp_Window, kM);

    /*
     * init opengl
     */
    glClearColor(0.2f, 0.2f, 0.2f, 1.f);

    /*
     * init protein loader
     */
    m_proteinLoader = ProteinLoader();
}



/*
 * INTERACTION
 */
void keyCallback(int key, int scancode, int action, int mods)
{
    if (key == GLFW_KEY_ESCAPE && action == GLFW_PRESS)
    {
        glfwSetWindowShouldClose(mp_Window, GL_TRUE);
    }

    if (key == GLFW_KEY_T && action == GLFW_PRESS)
    {
        m_resetStartMousePosition = true;
        m_interactionMode = 1; // set interaction mode to TRANSLATE
    }

    if (key == GLFW_KEY_R && action == GLFW_PRESS)
    {
        m_resetStartMousePosition = true;
        m_interactionMode = 2; // set interaction mode to ROTATE
    }
}

void mouseButtonCallback(int button, int action, int mods)
{
    if (button == GLFW_MOUSE_BUTTON_LEFT && action == GLFW_PRESS)
    {
        m_interactionMode = 0; // set interaction mode to NORMAL
        m_rotateCamera = true;
    }
    else if (button == GLFW_MOUSE_BUTTON_LEFT && action == GLFW_RELEASE)
    {
        m_rotateCamera = false;
    }

    if (button == GLFW_MOUSE_BUTTON_RIGHT && action == GLFW_PRESS)
    {
        if (m_interactionMode == 1 || m_interactionMode == 2)
        {
            // selected protein must be valid
            if (m_selectedProtein >= 0 && m_selectedProtein < m_proteinLoader.getNumberOfProteins())
            {
                if (m_interactionMode == 1)
                {
                    m_proteinLoader.getProteinAt(m_selectedProtein)->applyTranslationMatrix(m_tempModelMatrix);
                }
                else
                {
                    m_proteinLoader.getProteinAt(m_selectedProtein)->applyRotationMatrix(m_tempModelMatrix);
                }

                m_tempModelMatrix = glm::mat4();
                m_interactionMode = 0;
            }
        }
        else
        {
            // get the selected atom
            int atomIdx = -1;
            double cursorX, cursorY;
            glfwGetCursorPos(mp_Window, &cursorX, &cursorY);
            PickingTexture::PixelInfo pixel = m_pickingTexture.ReadPixel((uint)cursorX, HEIGHT - (uint)(cursorY) - 1);
            if (pixel.ObjectID >= 1)
            {
                atomIdx = (int)pixel.ObjectID;
            }
            m_selectedAtom = atomIdx;
            std::cout << "Selected id " << m_selectedAtom << std::endl;

            // get the corresponding protein
            if (atomIdx >= 0)
            {
                int currentAtomsCount = 0;
                int newProteinIdx = -1;
                for (int i = 0; i < m_proteinLoader.getNumberOfProteins(); i++)
                {
                    currentAtomsCount += m_proteinLoader.getProteinAt(i)->getNumberOfAtoms();
                    if (atomIdx < currentAtomsCount)
                    {
                        newProteinIdx = i;
                        break;
                    }
                }

                if (newProteinIdx >= 0)
                {
                    m_selectedProtein = newProteinIdx;
                }
            }
            std::cout << "Selected protein " << m_selectedProtein << std::endl;
        }
    }
}

void scrollCallback(double xoffset, double yoffset)
{
    mp_camera->setRadius(mp_camera->getRadius() - 2.f * (float)yoffset);
}

void mousePositionCallback(double xpos, double ypos)
{
    m_currentMousePosition = glm::vec2(xpos, HEIGHT - ypos);
}


/*
 * UPDATE
 */
void resetNeighborhoodSearch()
{
    glm::vec3 min, max;
    m_proteinLoader.getCenteredBoundingBoxAroundProteins(min, max);
    m_search.update(m_proteinLoader.getNumberOfAllAtoms(), min, max, m_gridResolution, m_searchRadius);
}

/*
 * DRAWING
 */
void drawSelectionRadius()
{
    glEnable(GL_BLEND);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

    m_drawSearchRadiusShader.use();
    SimpleAtom selectedAtom = m_proteinLoader.getAllAtoms().at(m_selectedAtom);
    SimpleProtein* selectedProtein = m_proteinLoader.getProteinAt(m_selectedProtein);
    glm::mat4 tempModelMatrix = selectedProtein->rotationMatrix;
    if (m_interactionMode == 1)
    {
        tempModelMatrix = m_tempModelMatrix * selectedProtein->translationMatrix * tempModelMatrix;
    }
    else if (m_interactionMode == 2)
    {
        tempModelMatrix = selectedProtein->translationMatrix * m_tempModelMatrix * tempModelMatrix;
    }
    else
    {
        tempModelMatrix = selectedProtein->translationMatrix * tempModelMatrix;
    }
    m_drawSearchRadiusShader.update("selectedAtomPosition", glm::vec3(tempModelMatrix * glm::vec4(selectedAtom.pos, 1.0)));
    m_drawSearchRadiusShader.update("searchRadius", m_searchRadius);
    m_drawSearchRadiusShader.update("view", mp_camera->getViewMatrix());
    m_drawSearchRadiusShader.update("projection", mp_camera->getProjectionMatrix());
    glDrawArrays(GL_POINTS, 0, 1);
    glBindVertexArray(0);

    glDisable(GL_BLEND);
}
void drawGizmo()
{
    /*
  * get the center of gravity for the selected protein in camera space
  */
    SimpleProtein* selectedProtein = m_proteinLoader.getProteinAt(m_selectedProtein);
    glm::vec3 proteinCenterOfGravity = selectedProtein->getCenterOfGravity();
    glm::vec4 proteinCenterOfGravityCameraSpace = mp_camera->getViewMatrix() * glm::vec4(proteinCenterOfGravity, 1.0);

    /*
     * unproject the cursor position into camera space
     */
    glm::vec4 viewport = glm::vec4(0.0, 0.0, WIDTH, HEIGHT);

    // project the center of gravity
    glm::vec4 proteinCenterOfGravityProjected = mp_camera->getProjectionMatrix() * proteinCenterOfGravityCameraSpace;
    glm::vec2 proteinCenterOfGravity2D;
    proteinCenterOfGravity2D.x = round(((proteinCenterOfGravityProjected.x/proteinCenterOfGravityProjected.w+1) /2.0) * WIDTH);
    proteinCenterOfGravity2D.y = HEIGHT - round(((proteinCenterOfGravityProjected.y/proteinCenterOfGravityProjected.w+1) /2.0) * HEIGHT);
    float proteinCenterOfGravityDepth = proteinCenterOfGravityProjected.z / proteinCenterOfGravityProjected.w;

    // reset mouse position if necessary
    if (m_resetStartMousePosition)
    {
        m_startMousePosition = m_currentMousePosition;
        m_resetStartMousePosition = false;
    }

    // current cursor position
    glm::vec3 pos3DCurrentCursor = mp_camera->getPositionAtPixel((int)m_currentMousePosition.x, (int)m_currentMousePosition.y);
    glm::vec3 cursorPosition = glm::vec3(m_currentMousePosition, proteinCenterOfGravityDepth);
    glm::vec3 unprojectedCursorPosition = glm::unProject(cursorPosition, mp_camera->getViewMatrix(), mp_camera->getProjectionMatrix(), viewport);
    glm::vec4 unprojectedCursorPositionCameraSpace = mp_camera->getViewMatrix() * glm::vec4(unprojectedCursorPosition, 1.0);
    // start cursor position
    glm::vec3 pos3DStartCursor = mp_camera->getPositionAtPixel((int)m_startMousePosition.x, (int)m_startMousePosition.y);
    glm::vec3 cursorStartPosition = glm::vec3(m_startMousePosition, proteinCenterOfGravityDepth);
    glm::vec3 unprojectedCursorStartPosition = glm::unProject(cursorStartPosition, mp_camera->getViewMatrix(), mp_camera->getProjectionMatrix(), viewport);
    glm::vec4 unprojectedCursorStartPositionCameraSpace = mp_camera->getViewMatrix() * glm::vec4(unprojectedCursorStartPosition, 1.0);

    // for the translation
    glm::vec4 unprojectedRestrictedPositionCameraSpace;
    glm::vec2 deltaDir = m_currentMousePosition - m_startMousePosition;
    if (abs(deltaDir.x) > abs(deltaDir.y))
    {
        glm::vec3 restrictedPosition = glm::vec3(m_currentMousePosition.x, m_startMousePosition.y, 0.2);
        glm::vec3 unprojectedRestrictedPosition = glm::unProject(restrictedPosition, mp_camera->getViewMatrix(), mp_camera->getProjectionMatrix(), viewport);
        unprojectedRestrictedPositionCameraSpace = mp_camera->getViewMatrix() * glm::vec4(unprojectedRestrictedPosition, 1.0);
    }
    else
    {
        glm::vec3 restrictedPosition = glm::vec3(m_startMousePosition.x, m_currentMousePosition.y, 0.2);
        glm::vec3 unprojectedRestrictedPosition = glm::unProject(restrictedPosition, mp_camera->getViewMatrix(), mp_camera->getProjectionMatrix(), viewport);
        unprojectedRestrictedPositionCameraSpace = mp_camera->getViewMatrix() * glm::vec4(unprojectedRestrictedPosition, 1.0);
    }

    /*
     * upload points to gpu
     */
    std::vector<glm::vec4> points;
    if (m_interactionMode == 1)
    {
        points.push_back(unprojectedCursorStartPositionCameraSpace);
        points.push_back(unprojectedRestrictedPositionCameraSpace);
    }
    else if (m_interactionMode == 2)
    {
        points.push_back(proteinCenterOfGravityCameraSpace);
        points.push_back(unprojectedCursorStartPositionCameraSpace);
        points.push_back(proteinCenterOfGravityCameraSpace);
        points.push_back(unprojectedCursorPositionCameraSpace);
    }

    GLuint pointsVBO = 0;
    glGenBuffers(1, &pointsVBO);
    glBindBuffer(GL_ARRAY_BUFFER, pointsVBO);
    glBufferData(GL_ARRAY_BUFFER, points.size()*4*sizeof(float), points.data(), GL_STATIC_DRAW);
    GLuint pointsVAO = 0;
    glGenVertexArrays(1, &pointsVAO);
    glBindVertexArray(pointsVAO);
    glVertexAttribPointer(0, 4, GL_FLOAT, GL_FALSE, 0, 0);
    // disable vbo
    glEnableVertexAttribArray(0);

    /*
     * draw a line between the center of gravity and the unprojected mouse position
     */
    glEnable(GL_BLEND);
    glDisable(GL_DEPTH_TEST);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
    glLineWidth(3.0);
    m_linesShader.use();
    m_linesShader.update("projMat", mp_camera->getProjectionMatrix());
    glDrawArrays(GL_LINES, 0, points.size());
    glLineWidth(1.0);
    glEnable(GL_DEPTH_TEST);
    glDisable(GL_BLEND);

    // delete buffers
    glBindVertexArray(0);
    glDeleteBuffers(1, &pointsVBO);
    glDeleteVertexArrays(1, &pointsVAO);


    /*
     * handle the interaction
     */
    // apply the right interaction
    if (m_interactionMode == 1) // TRANSLATION
    {
        glm::vec3 deltaDir = unprojectedCursorPositionCameraSpace - unprojectedCursorStartPositionCameraSpace;
        if (abs(deltaDir.x) > abs(deltaDir.y))
        {
            glm::vec3 at = glm::normalize(mp_camera->getCenter() - mp_camera->getPosition());
            glm::vec3 right = glm::normalize(glm::cross(at, glm::vec3(0.0, 1.0, 0.0)));
            m_tempModelMatrix = glm::translate(right*deltaDir.x);
        }
        else
        {
            glm::vec3 at = glm::normalize(mp_camera->getCenter() - mp_camera->getPosition());
            glm::vec3 right = glm::normalize(glm::cross(at, glm::vec3(0.0, 1.0, 0.0)));
            glm::vec3 up = glm::normalize(glm::cross(right, at));
            m_tempModelMatrix = glm::translate(up*deltaDir.y);
        }
    }

    if (m_interactionMode == 2) // ROTATION
    {
        glm::vec2 v1 = m_startMousePosition - proteinCenterOfGravity2D;
        v1 = glm::normalize(v1);
        glm::vec2 v2 = m_currentMousePosition - proteinCenterOfGravity2D;
        v2 = glm::normalize(v2);

        float angle = glm::orientedAngle(v1, v2); // * 180 / M_PI;

        m_tempModelMatrix = glm::rotate(angle, mp_camera->getPosition() - mp_camera->getCenter());
    }
}
void updateGrid()
{
    /*
     * update the grid dimensions
     */
    glm::vec3 min, max;
    m_search.getGridMinMax(min, max);

    float cellSize = m_search.getCellSize();
    glm::ivec3 gridRes = m_search.getGridResolution();



    /*
     * update the visual grid representation
     */
    /*
     * grid hull
     */
    std::vector<glm::vec4> points;
    glm::vec4 p1 = glm::vec4(min.x, min.y, min.z, 1);
    glm::vec4 p2 = glm::vec4(max.x, min.y, min.z, 1);
    glm::vec4 p3 = glm::vec4(min.x, min.y, max.z, 1);
    glm::vec4 p4 = glm::vec4(max.x, min.y, max.z, 1);
    glm::vec4 p5 = glm::vec4(min.x, max.y, min.z, 1);
    glm::vec4 p6 = glm::vec4(max.x, max.y, min.z, 1);
    glm::vec4 p7 = glm::vec4(min.x, max.y, max.z, 1);
    glm::vec4 p8 = glm::vec4(max.x, max.y, max.z, 1);
    // bottom
    points.push_back(p1);
    points.push_back(p2);
    points.push_back(p1);
    points.push_back(p3);
    points.push_back(p2);
    points.push_back(p4);
    points.push_back(p3);
    points.push_back(p4);
    // side
    points.push_back(p1);
    points.push_back(p5);
    points.push_back(p2);
    points.push_back(p6);
    points.push_back(p3);
    points.push_back(p7);
    points.push_back(p4);
    points.push_back(p8);
    // top
    points.push_back(p5);
    points.push_back(p6);
    points.push_back(p5);
    points.push_back(p7);
    points.push_back(p6);
    points.push_back(p8);
    points.push_back(p7);
    points.push_back(p8);

    m_numVBOEntries = 24;

    /*
     * grid cell divisions
     */
    float xs[2] = {min.x, max.x};
    float ys[2] = {min.y, max.y};
    float zs[2] = {min.z, max.z};

    // x resolution
    for (int i = 0; i < gridRes.x; i++) {
        float x = min.x + (i * cellSize);
        for (int j = 0; j < 2; j++) {
            for (int k = 0; k < 2; k++) {
                float y = ys[j];
                float z = zs[k];
                points.push_back(glm::vec4(x, y, z, 1));
                m_numVBOEntries++;
            }
        }
        for (int j = 0; j < 2; j++) {
            for (int k = 0; k < 2; k++) {
                float y = ys[k];
                float z = zs[j];
                points.push_back(glm::vec4(x, y, z, 1));
                m_numVBOEntries++;
            }
        }
    }

    // y resolution
    for (int i = 0; i < gridRes.y; i++) {
        float y = min.y + (i * cellSize);
        for (int j = 0; j < 2; j++) {
            for (int k = 0; k < 2; k++) {
                float x = xs[j];
                float z = zs[k];
                points.push_back(glm::vec4(x, y, z, 1));
                m_numVBOEntries++;
            }
        }
        for (int j = 0; j < 2; j++) {
            for (int k = 0; k < 2; k++) {
                float x = xs[k];
                float z = zs[j];
                points.push_back(glm::vec4(x, y, z, 1));
                m_numVBOEntries++;
            }
        }
    }

    // z resolution
    for (int i = 0; i < gridRes.z; i++) {
        float z = min.z + (i * cellSize);
        for (int j = 0; j < 2; j++) {
            for (int k = 0; k < 2; k++) {
                float x = xs[j];
                float y = ys[k];
                points.push_back(glm::vec4(x, y, z, 1));
                m_numVBOEntries++;
            }
        }
        for (int j = 0; j < 2; j++) {
            for (int k = 0; k < 2; k++) {
                float x = xs[k];
                float y = ys[j];
                points.push_back(glm::vec4(x, y, z, 1));
                m_numVBOEntries++;
            }
        }
    }

    if (m_pointsVBO != 0) glDeleteBuffers(1, &m_pointsVBO);
    m_pointsVBO = 0;
    glGenBuffers(1, &m_gridVBO);
    glBindBuffer(GL_ARRAY_BUFFER, m_gridVBO);
    glBufferData(GL_ARRAY_BUFFER, points.size() * 4 * sizeof(float), points.data(), GL_STATIC_DRAW);

    if (m_pointsVAO != 0) glDeleteVertexArrays(1, &m_gridVAO);
    m_gridVAO = 0;
    glGenVertexArrays(1, &m_gridVAO);
    glBindVertexArray(m_gridVAO);
    glVertexAttribPointer(0, 4, GL_FLOAT, GL_FALSE, 0, 0);

    // disable the vao and vbo
    glEnableVertexAttribArray(0);
    glBindVertexArray(0);
}
void drawGrid()
{
    glEnable(GL_BLEND);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

    glBindVertexArray(m_gridVAO);

    m_linesShader.use();
    //m_linesShader.update("viewMat", mp_camera->getViewMatrix());
    m_linesShader.update("projMat", mp_camera->getProjectionMatrix() * mp_camera->getViewMatrix());
    glDrawArrays(GL_LINES, 0, m_numVBOEntries);

    glBindVertexArray(0);

    glDisable(GL_BLEND);
}


/*
 * SETUP AND MAINLOOP
 */
void run()
{
    /*
     * setup camera
     */
    glm::vec3 cameraCenter = glm::vec3(0.0, 0.0, 0.0);
    float cameraRadius = 0.0f;
    glm::vec3 globalMin;
    glm::vec3 globalMax;
    m_proteinLoader.getCenteredBoundingBoxAroundProteins(globalMin, globalMax);
    cameraCenter = (globalMax + globalMin) / 2;
    float radius = glm::length(globalMax - globalMin);
    cameraRadius = (radius > cameraRadius) ? radius : cameraRadius;

    mp_camera = std::unique_ptr<OrbitCamera>(
            new OrbitCamera(
                    cameraCenter,
                    90.0f,
                    90.0f,
                    cameraRadius,
                    cameraRadius / 2.0f,
                    5.0f * cameraRadius,
                    45.0f,
                    0.05f
            )
    );


    /*
     * setup cursor position
     */
    float prevCursorX, prevCursorY = 0;


    /*
     * setup neighborhood
     */
    Neighborhood neighborhood;



    /*
     * renderloop
     * call render function of Renderer.h and provide it with a function
     */
    render(mp_Window, [&] (float deltaTime)
    {
        /*
         * clear buffer
         */
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);


        /*
         * request new frame for ImGui
         */
        ImGui_ImplGlfwGL3_NewFrame();

        /*
         * update viewport
         */
        glm::vec2 resolution = getResolution(mp_Window);
        glViewport(0,0,resolution.x, resolution.y);

        /*
         * Calculate cursor movement
         */
        double cursorX, cursorY;
        glfwGetCursorPos(mp_Window, &cursorX, &cursorY);
        GLfloat cursorDeltaX = (float)cursorX - prevCursorX;
        GLfloat cursorDeltaY = (float)cursorY - prevCursorY;
        prevCursorX = cursorX;
        prevCursorY = cursorY;

        /*
         * update camera
         */
        if(m_rotateCamera) {
            m_CameraDeltaMovement = glm::vec2(cursorDeltaX, cursorDeltaY);
            m_CameraSmoothTime = 1.f;
        }
        glm::vec2 cameraMovement = glm::lerp(glm::vec2(0), m_CameraDeltaMovement, m_CameraSmoothTime);
        mp_camera->setAlpha(mp_camera->getAlpha() + 0.25f * cameraMovement.x);
        mp_camera->setBeta(mp_camera->getBeta() - 0.25f * cameraMovement.y);
        mp_camera->update(resolution.x, resolution.y, true);

        /*
         * set light direction
         */
        //m_lightDirection = glm::normalize(-mp_camera->getPosition() + mp_camera->getCenter()); // light direction equals view direction of the camera
        m_lightDirection = glm::normalize(glm::vec3(0, -1, 0));

        /*
         * update model matrices
         */
        std::vector<glm::mat4> modelMatrices;
        for (int i = 0; i < m_proteinLoader.getNumberOfProteins(); i++)
        {
            if (m_selectedProtein == i && m_interactionMode == 1)
            {
                modelMatrices.push_back(m_tempModelMatrix * m_proteinLoader.getProteinAt(i)->translationMatrix * m_proteinLoader.getProteinAt(i)->rotationMatrix);
            }
            else if (m_selectedProtein == i && m_interactionMode == 2)
            {
                modelMatrices.push_back(m_proteinLoader.getProteinAt(i)->translationMatrix * m_tempModelMatrix *  m_proteinLoader.getProteinAt(i)->rotationMatrix);
            }
            else
            {
                modelMatrices.push_back(m_proteinLoader.getProteinAt(i)->translationMatrix * m_proteinLoader.getProteinAt(i)->rotationMatrix);
            }
        }

        glBindBuffer(GL_SHADER_STORAGE_BUFFER, m_modelMatricesBuffer);
        GLvoid* p_modelMatrices = glMapBuffer(GL_SHADER_STORAGE_BUFFER, GL_WRITE_ONLY);
        memcpy(p_modelMatrices, modelMatrices.data(), sizeof(glm::mat4)*modelMatrices.size());
        glUnmapBuffer(GL_SHADER_STORAGE_BUFFER);

        /*
         * update atom positions
         */
        glBindBuffer(GL_SHADER_STORAGE_BUFFER, m_atomsSSBO);
        GLvoid* p = glMapBuffer(GL_SHADER_STORAGE_BUFFER, GL_WRITE_ONLY);
        memcpy(p, m_proteinLoader.getAllAtoms().data(), sizeof(SimpleAtom)*m_proteinLoader.getNumberOfAllAtoms());
        glUnmapBuffer(GL_SHADER_STORAGE_BUFFER);

        /*
         * extract positions from atoms
         */
        int numBlocks = ceil((float)m_proteinLoader.getNumberOfAllAtoms() / BLOCK_SIZE);
        glBindBufferBase(GL_SHADER_STORAGE_BUFFER, 0, m_atomsSSBO);
        glBindBufferBase(GL_SHADER_STORAGE_BUFFER, 1, m_positionsSSBO);
        glBindBufferBase(GL_SHADER_STORAGE_BUFFER, 2, m_modelMatricesBuffer);
        m_extractAtomPositionsShader.use();
        m_extractAtomPositionsShader.update("pnum",m_proteinLoader.getNumberOfAllAtoms());
        glDispatchCompute(numBlocks,1,1);
        glMemoryBarrier(GL_ALL_BARRIER_BITS);
        glBindBufferBase(GL_SHADER_STORAGE_BUFFER, 0, 0);
        glBindBufferBase(GL_SHADER_STORAGE_BUFFER, 1, 0);
        glBindBufferBase(GL_SHADER_STORAGE_BUFFER, 2, 0);

        /*
         * find neighborhood
         */
        resetNeighborhoodSearch();
        m_search.run(&m_positionsSSBO, neighborhood);



        /*
         * bind relevant buffers
         */
        glBindBufferBase(GL_SHADER_STORAGE_BUFFER, 0, m_atomsSSBO);
        glBindBufferBase(GL_SHADER_STORAGE_BUFFER, 1, m_indexCubeSSBO);
        glBindBufferBase(GL_SHADER_STORAGE_BUFFER, 2, m_valuesVectorSSBO);
        glBindBufferBase(GL_SHADER_STORAGE_BUFFER, 3, m_bondNeighborsSSBO);

        glBindBufferBase(GL_SHADER_STORAGE_BUFFER, 4, *neighborhood.dp_particleCell);
        glBindBufferBase(GL_SHADER_STORAGE_BUFFER, 5, *neighborhood.dp_particleCellIndex);
        glBindBufferBase(GL_SHADER_STORAGE_BUFFER, 6, *neighborhood.dp_gridCellCounts);
        glBindBufferBase(GL_SHADER_STORAGE_BUFFER, 7, *neighborhood.dp_gridCellOffsets);
        glBindBufferBase(GL_SHADER_STORAGE_BUFFER, 8, *neighborhood.dp_grid);
        glBindBufferBase(GL_SHADER_STORAGE_BUFFER, 9, *neighborhood.dp_particleOriginalIndex);

        glBindBufferBase(GL_SHADER_STORAGE_BUFFER, 10, m_modelMatricesBuffer);
        glBindBufferBase(GL_SHADER_STORAGE_BUFFER, 11, m_proteinColorsBuffer);


        /*
         * update charge
         */
        m_calcChargeProgram.use();
        m_calcChargeProgram.update("numAtoms", m_proteinLoader.getNumberOfAllAtoms());
        m_calcChargeProgram.update("numAtomSymbols", (int)m_numberOfAtomSymbols);
        m_calcChargeProgram.update("searchRadius2", m_searchRadius*m_searchRadius);
        m_calcChargeProgram.update("gridAdjCnt", neighborhood.numberOfSearchCells);
        m_calcChargeProgram.update("searchCellOff", neighborhood.startCellOffset);
        m_calcChargeProgram.update("useSimplifiedCalculation", useSimplifiedCalculation);
        m_calcChargeProgram.update("selectedAtomID", m_selectedAtom);
        glDispatchCompute(ceil(m_proteinLoader.getNumberOfAllAtoms()/256.0), 1, 1);
        glUniform1iv(glGetUniformLocation(m_calcChargeProgram.getProgramHandle(),"gridAdj"), 216, neighborhood.p_searchCellOffsets);
        glMemoryBarrier (GL_ALL_BARRIER_BITS);

        /*
         * draw proteins as impostor
         */
        glBindBufferBase(GL_SHADER_STORAGE_BUFFER, 0, m_atomsSSBO);
        m_impostorProgram.use();
        m_impostorProgram.update("view", mp_camera->getViewMatrix());
        m_impostorProgram.update("projection", mp_camera->getProjectionMatrix());
        m_impostorProgram.update("cameraWorldPos", mp_camera->getPosition());
        m_impostorProgram.update("probeRadius", 0.f);
        m_impostorProgram.update("lightDir", m_lightDirection);
        m_impostorProgram.update("selectedProtein", m_selectedProtein);
        m_impostorProgram.update("drawSelectedProteinOnly", m_drawSelectedProtein);
        glDrawArrays(GL_POINTS, 0, (GLsizei)m_proteinLoader.getNumberOfAllAtoms());


        /*
         * fill id picking texture with atom ids
         */
        m_pickingTexture.EnableWriting();
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
        m_idPickingShader.use();
        m_idPickingShader.update("view", mp_camera->getViewMatrix());
        m_idPickingShader.update("projection", mp_camera->getProjectionMatrix());
        m_idPickingShader.update("probeRadius", 0.f);
        glDrawArrays(GL_POINTS, 0, (GLsizei)m_proteinLoader.getNumberOfAllAtoms());
        m_pickingTexture.DisableWriting();

        /*
         * unbind buffers
         */
        glBindBufferBase(GL_SHADER_STORAGE_BUFFER, 0, 0);
        glBindBufferBase(GL_SHADER_STORAGE_BUFFER, 1, 0);
        glBindBufferBase(GL_SHADER_STORAGE_BUFFER, 2, 0);
        glBindBufferBase(GL_SHADER_STORAGE_BUFFER, 3, 0);
        glBindBufferBase(GL_SHADER_STORAGE_BUFFER, 4, 0);
        glBindBufferBase(GL_SHADER_STORAGE_BUFFER, 5, 0);
        glBindBufferBase(GL_SHADER_STORAGE_BUFFER, 6, 0);
        glBindBufferBase(GL_SHADER_STORAGE_BUFFER, 7, 0);
        glBindBufferBase(GL_SHADER_STORAGE_BUFFER, 8, 0);
        glBindBufferBase(GL_SHADER_STORAGE_BUFFER, 9, 0);
        glBindBufferBase(GL_SHADER_STORAGE_BUFFER, 10, 0);
        glBindBufferBase(GL_SHADER_STORAGE_BUFFER, 11, 0);



        /*
         * draw radius around selected atom
         */
        if (m_selectedAtom >= 0 && m_selectedAtom < m_proteinLoader.getNumberOfAllAtoms() && m_selectedProtein >= 0 && m_selectedProtein < m_proteinLoader.getNumberOfProteins())
        {
            drawSelectionRadius();
        }



        /*
         * drawing the gizmo
         */
        // only do something when a valid protein is selected
        if (m_selectedProtein >= 0 && m_selectedProtein < m_proteinLoader.getNumberOfProteins() && (m_interactionMode == 1 || m_interactionMode == 2))
        {
            drawGizmo();
        }

        /*
         * draw grid
         */
        if (m_isGridVisible)
        {
            updateGrid();
            drawGrid();
        }


        /*
         * update gui
         */
        updateGUI();
    });
}



/*
 * IMGUI
 */
void updateGUI()
{
    // Main menu bar
    if (ImGui::BeginMainMenuBar())
    {
        /*
         * General menu
         */
        if (ImGui::BeginMenu("Menu"))
        {
            if(ImGui::MenuItem("Quit", "Esc", false, true)) { glfwSetWindowShouldClose(mp_Window, GL_TRUE); }
            ImGui::EndMenu();
        }

        /*
         * Shortcut infos
         */
        if (ImGui::BeginMenu("Controls"))
        {
            ImGui::Text("Rotate the camera by holding left mouse button while moving the mouse");
            ImGui::Text("Select an atom with right click");
            ImGui::Text("T: Translate selected protein");
            ImGui::Text("R: Rotate selected protein");
            ImGui::EndMenu();
        }

        /*
         * Protein infos
         */
        if (ImGui::BeginMenu("Proteins"))
        {
            if (m_proteinLoader.getNumberOfProteins() > 0) {
                for (int i = 0; i < m_proteinLoader.getNumberOfProteins(); i++) {
                    SimpleProtein* protein = m_proteinLoader.getProteinAt(i);
                    std::string text = std::to_string(i) + ": " + protein->name + " - atom number: " + std::to_string(protein->atoms.size());
                    ImGui::Text(text.c_str());
                }
            } else {
                ImGui::Text("No proteins loaded!");
            }
            ImGui::Separator();
            std::string atomText = "Total number of atoms: " + std::to_string(m_proteinLoader.getNumberOfAllAtoms());
            ImGui::Text(atomText.c_str());
            std::string drawSelectedProteinText = "Draw only selected Protein";
            ImGui::Checkbox(drawSelectedProteinText.c_str(), &m_drawSelectedProtein);

            ImGui::EndMenu();
        }

        /*
         * Neighborhood search
         */
        if (ImGui::BeginMenu("NeighborhoodSearch"))
        {
            std::string isGridVisibleText = "Grid visible";
            ImGui::Checkbox(isGridVisibleText.c_str() ,&m_isGridVisible);
            std::string searchRadiusSize = "Suchradius in Angström";
            ImGui::SliderFloat(searchRadiusSize.c_str(), &m_searchRadius, 0.1, 10.0);
            ImGui::EndMenu();
        }

        /*
         * Amber force field
         */
        if (ImGui::BeginMenu("AmberForcefield"))
        {
            std::string useSimplifiedCalculationText = "Use simplified force calculation";
            ImGui::Checkbox(useSimplifiedCalculationText.c_str(), &useSimplifiedCalculation);

            std::string selectedAtomText = "Selected atom: " + std::to_string(m_selectedAtom);
            std::string selectedProteinText = "Selected protein: " + std::to_string(m_selectedProtein);
            ImGui::Text(selectedAtomText.c_str());
            ImGui::Text(selectedProteinText.c_str());

            ImGui::EndMenu();
        }



        /*
         * Frametime
         */
        float framerate = ImGui::GetIO().Framerate;
        std::stringstream stream;
        stream << std::fixed << std::setprecision(0) << framerate;
        std::string fps = "FPS: " + stream.str();
        ImGui::MenuItem(fps.c_str(), "", false, false);

        /*
         * End main menu bar
         */
        ImGui::EndMainMenuBar();
    }

    ImGui::Render();
}



int main()
{
    Logger::instance().changeTab("     ");
    Logger::instance().print("Start Amber force field visualization!"); Logger::instance().tabIn();

    /*
     * load general amber force field parameter
     */
    AmberForceFieldParameter amberForceFieldParameter;
    std::string basePath = RESOURCES_PATH;
    std::string filepath = basePath + "/forcefieldParameters/" + "gaff.dat";
    AmberForceFieldParameterLoader amberForceFieldParameterLoader;
    amberForceFieldParameterLoader.load(filepath, amberForceFieldParameter);
    m_numberOfAtomSymbols = amberForceFieldParameter.getNumberOfAtomSymbols();

    /*
     * setup OpenGL
     */
    setup();

    /*
     * compile shader programs
     */
    m_impostorProgram            = ShaderProgram("/AmberForceFieldVisualization/renderAtoms/impostor.vert", "/AmberForceFieldVisualization/renderAtoms/impostor.geom", "/AmberForceFieldVisualization/renderAtoms/impostor.frag");
    m_calcChargeProgram          = ShaderProgram("/AmberForceFieldVisualization/calculateForceField/amberff.comp");
    m_extractAtomPositionsShader = ShaderProgram("/AmberForceFieldVisualization/calculateForceField/extractAtomPositions.comp");
    m_drawSearchRadiusShader     = ShaderProgram("/NeighborSearch/renderSearchRadius/radius.vert", "/NeighborSearch/renderSearchRadius/radius.geom", "/NeighborSearch/renderSearchRadius/radius.frag");
    m_idPickingShader            = ShaderProgram("/AmberForceFieldVisualization/pickingTexture/impostor.vert", "/AmberForceFieldVisualization/pickingTexture/impostor.geom", "/AmberForceFieldVisualization/pickingTexture/impostor.frag");
    m_linesShader                = ShaderProgram("/AmberForceFieldVisualization/renderLinesCameraSpace/lines.vert", "/AmberForceFieldVisualization/renderLinesCameraSpace/lines.frag");
    GLenum err = glGetError();
    if (err != GL_NO_ERROR) {
        //Logger::instance().print("GLerror after init shader programs: " + std::to_string(err), Logger::Mode::ERROR);
    }

    /*
     * load proteins
     */
    SimpleProtein* proteinA = m_proteinLoader.loadProtein("PDB/3ah8 protein.pdb", amberForceFieldParameter.getAtomSymbolsMap());
    m_proteinColors.push_back(glm::vec4(0.8, 0.0, 0.4, 1.0));
    SimpleProtein* proteinB = m_proteinLoader.loadProtein("PDB/3ah8 ligand.pdb", amberForceFieldParameter.getAtomSymbolsMap());
    m_proteinColors.push_back(glm::vec4(0.0, 0.2, 8.0, 1.0));

    proteinA->center();
    proteinB->center();

    /*
     * setup neighborhood search
     */
    glm::fvec3 min, max;
    m_proteinLoader.getCenteredBoundingBoxAroundProteins(min, max);
    m_search.init(m_proteinLoader.getNumberOfAllAtoms(),min,max,m_gridResolution,m_searchRadius);

    /*
     * setup buffers
     */
    m_atomsSSBO;
    glGenBuffers(1, &m_atomsSSBO);
    glBindBuffer(GL_SHADER_STORAGE_BUFFER, m_atomsSSBO);
    glBufferData(GL_SHADER_STORAGE_BUFFER, sizeof(SimpleAtom) * m_proteinLoader.getNumberOfAllAtoms(), m_proteinLoader.getAllAtoms().data(), GL_STATIC_DRAW);
    glBindBufferBase(GL_SHADER_STORAGE_BUFFER, 0, m_atomsSSBO);
    glBindBuffer(GL_SHADER_STORAGE_BUFFER, 0);

    glGenBuffers(1, &m_indexCubeSSBO);
    glBindBuffer(GL_SHADER_STORAGE_BUFFER, m_indexCubeSSBO);
    glBufferData(GL_SHADER_STORAGE_BUFFER, sizeof(uint) * amberForceFieldParameter.getIndexCube().size(), amberForceFieldParameter.getIndexCube().data(), GL_STATIC_DRAW);
    glBindBufferBase(GL_SHADER_STORAGE_BUFFER, 1, m_indexCubeSSBO);
    glBindBuffer(GL_SHADER_STORAGE_BUFFER, 0);

    glGenBuffers(1, &m_valuesVectorSSBO);
    glBindBuffer(GL_SHADER_STORAGE_BUFFER, m_valuesVectorSSBO);
    glBufferData(GL_SHADER_STORAGE_BUFFER, sizeof(uint) * amberForceFieldParameter.getValuesList().size(), amberForceFieldParameter.getValuesList().data(), GL_STATIC_DRAW);
    glBindBufferBase(GL_SHADER_STORAGE_BUFFER, 2, m_valuesVectorSSBO);
    glBindBuffer(GL_SHADER_STORAGE_BUFFER, 0);

    glGenBuffers(1, &m_bondNeighborsSSBO);
    glBindBuffer(GL_SHADER_STORAGE_BUFFER, m_bondNeighborsSSBO);
    glBufferData(GL_SHADER_STORAGE_BUFFER, sizeof(uint) * m_proteinLoader.getAllNeighbors().size(), m_proteinLoader.getAllNeighbors().data(), GL_STATIC_DRAW);
    glBindBufferBase(GL_SHADER_STORAGE_BUFFER, 3, m_bondNeighborsSSBO);
    glBindBuffer(GL_SHADER_STORAGE_BUFFER, 0);

    GPUHandler::initSSBO<glm::vec4>(&m_positionsSSBO, m_proteinLoader.getNumberOfAllAtoms());
    GPUHandler::initSSBO<int>(&m_searchResultsSSBO, m_proteinLoader.getNumberOfAllAtoms());
    GPUHandler::initSSBO<glm::mat4>(&m_modelMatricesBuffer, m_proteinLoader.getNumberOfProteins());
    GPUHandler::initSSBO<glm::vec4>(&m_proteinColorsBuffer, m_proteinLoader.getNumberOfProteins());

    glBindBuffer(GL_SHADER_STORAGE_BUFFER, m_proteinColorsBuffer);
    GLvoid* p_proteinColors = glMapBuffer(GL_SHADER_STORAGE_BUFFER, GL_WRITE_ONLY);
    memcpy(p_proteinColors, m_proteinColors.data(), sizeof(glm::vec4)*m_proteinColors.size());
    glUnmapBuffer(GL_SHADER_STORAGE_BUFFER);

    /*
     * setup id picking texture
     */
    m_pickingTexture.Init(WIDTH, HEIGHT);



    /*
     * run application
     */
    run();

    Logger::instance().tabOut(); Logger::instance().print("Exit Amber force field visualization!");

    return 0;
}
