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

// local includes
#include "ProteinLoader.h"
#include "Utils/OrbitCamera.h"
#include "AmberForceFieldParameterLoader.h"






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

// protein
ProteinLoader   m_proteinLoader;
int             m_selectedAtom = -1;
int             m_selectedProtein = 0;
float           m_proteinMoveSpeed = 2.f;
GLuint m_atomsSSBO;
GLuint m_indexCubeSSBO;
GLuint m_valuesVectorSSBO;
GLuint m_bondNeighborsSSBO;
ShaderProgram m_calcChargeProgram;
uint m_numberOfAtomSymbols;

// rendering
glm::vec3   m_lightDirection;
ShaderProgram m_impostorProgram;
GLuint m_pointsVBO;
GLuint m_pointsVAO;
int    m_numVBOEntries;






/*
 * forward declarations
 */
void setup();
void compileShaderPrograms();

void keyCallback(int key, int scancode, int action, int mods);
void mouseButtonCallback(int button, int action, int mods);
void scrollCallback(double xoffset, double yoffset);

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

}

void mouseButtonCallback(int button, int action, int mods)
{
    if (button == GLFW_MOUSE_BUTTON_LEFT && action == GLFW_PRESS)
    {
        m_rotateCamera = true;
    }
    else if (button == GLFW_MOUSE_BUTTON_LEFT && action == GLFW_RELEASE)
    {
        m_rotateCamera = false;
    }
}

void scrollCallback(double xoffset, double yoffset)
{
    mp_camera->setRadius(mp_camera->getRadius() - 2.f * (float)yoffset);
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
         * update atom positions
         */
        glBindBuffer(GL_SHADER_STORAGE_BUFFER, m_atomsSSBO);
        GLvoid* p = glMapBuffer(GL_SHADER_STORAGE_BUFFER, GL_WRITE_ONLY);
        memcpy(p, m_proteinLoader.getAllAtoms().data(), sizeof(SimpleAtom)*m_proteinLoader.getNumberOfAllAtoms());
        glUnmapBuffer(GL_SHADER_STORAGE_BUFFER);

        /*
         * update charge
         */
        m_calcChargeProgram.use();
        m_calcChargeProgram.update("numAtoms", m_proteinLoader.getNumberOfAllAtoms());
        //m_calcChargeProgram.update("numAtomSymbols", (int)m_numberOfAtomSymbols);
        glDispatchCompute(ceil(m_proteinLoader.getNumberOfAllAtoms()/256.0), 1, 1);
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
        //m_impostorProgram.update("proteinNum", (int)m_proteinLoader.getNumberOfProteins());
        glDrawArrays(GL_POINTS, 0, (GLsizei)m_proteinLoader.getNumberOfAllAtoms());

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
            ImGui::Text("P: Switch between proteins");
            ImGui::Text("O: Switch between atoms");
            ImGui::Text("W: Move selected protein up");
            ImGui::Text("A: Move selected protein left");
            ImGui::Text("S: Move selected protein down");
            ImGui::Text("D: Move selected protein right");
            ImGui::Text("Q: Move selected protein back");
            ImGui::Text("E: Move selected protein forth");
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
     * general setup
     */
    setup();

    /*
     * compile shader programs
     */
    m_impostorProgram     = ShaderProgram("/AmberForceFieldVisualization/impostor.vert", "/AmberForceFieldVisualization/impostor.geom", "/AmberForceFieldVisualization/impostor.frag");
    m_calcChargeProgram   = ShaderProgram("/AmberForceFieldVisualization/amberff.comp");
    GLenum err = glGetError();
    if (err != GL_NO_ERROR) {
        //Logger::instance().print("GLerror after init shader programs: " + std::to_string(err), Logger::Mode::ERROR);
    }

    /*
     * load proteins
     */
    SimpleProtein* proteinA = m_proteinLoader.loadProtein("PDB/1a19.pdb", amberForceFieldParameter.getAtomSymbolsMap());
    SimpleProtein* proteinB = m_proteinLoader.loadProtein("PDB/5bs0.pdb", amberForceFieldParameter.getAtomSymbolsMap());
    proteinA->center();
    proteinB->center();
    proteinB->move(glm::vec3(proteinA->extent().x/2 + proteinB->extent().x/2, 0, 0));

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

    /*
     * run application
     */
    run();

    Logger::instance().tabOut(); Logger::instance().print("Exit Amber force field visualization!");

    return 0;
}
