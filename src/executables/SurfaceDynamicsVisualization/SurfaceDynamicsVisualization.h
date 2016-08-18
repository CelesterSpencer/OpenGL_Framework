//============================================================================
// Distributed under the MIT License. Author: Raphael Menges
//============================================================================

#ifndef SURFACE_DYNAMICS_VISUALIZATION_H
#define SURFACE_DYNAMICS_VISUALIZATION_H

#include <memory>
#include <GL/glew.h>
#include <GLFW/glfw3.h>
#include <glm/glm.hpp>
#include <set>

#include "ShaderTools/ShaderProgram.h"
#include "SurfaceExtraction/GPUProtein.h"
#include "SurfaceExtraction/GPUSurfaceExtraction.h"
#include "SurfaceExtraction/SurfaceValidation.h"
#include "Framebuffer.h"
#include "Path.h"
#include "SurfaceExtraction/GPUHullSamples.h"

// Notes:
// - Calculations done in angstrom

// Forward declaration
class Protein;
class GPUProtein;
class OrbitCamera;

// Class
class SurfaceDynamicsVisualization
{
public:

    // Constructor
    SurfaceDynamicsVisualization(std::string filepathPDB, std::string filepathXTC = "");

    // Destructor
    virtual ~SurfaceDynamicsVisualization();

    // Render
    void renderLoop();

    // Set window title displaying progress (hides processing part when progress is equal one)
    void setProgressDisplay(std::string task, float progress = 1.0f);

private:

    // Enumeration for surface rendering
    enum SurfaceRendering
    {
        HULL, ASCENSION, ELEMENTS, AMINOACIDS, ANALYSIS
    };

    // Enumeration for background cubemaps
    enum Background
    {
        SCIENTIFIC, COMPUTERVISUALISTIK, BEACH
    };

    // Keyboard callback for GLFW
    void keyCallback(int key, int scancode, int action, int mods);

    // Mouse button callback for GLFW
    void mouseButtonCallback(int button, int action, int mods);

    // Scroll callback for GLFW
    void scrollCallback(double xoffset, double yoffset);

    // Render GUI
    void renderGUI();

    // Update computation information
    void updateComputationInformation(std::string device, float computationTime);

    // Set frame. Returns whether frame has been changed
    bool setFrame(int frame);

    // Compute layers
    void computeLayers(bool useGPU);

    // Compute hull samples
    void computeHullSamples();

    // Compute ascension
    void computeAscension();

    // Get atom beneath cursor. Returns -1 when fails
    int getAtomBeneathCursor() const;

    // Calculate approximated surface of molecule
    float approximateSurfaceArea(std::vector<GLuint> indices, int frame) const;

    // Update global analysis
    void updateGlobalAnalysis();

    // Update group analysis
    void updateGroupAnalysis();

    // Get whether frame was computed (otherwise prohibit doing thing which would go wrong)
    bool frameComputed() const { return (mFrame >= mComputedStartFrame) && (mFrame <= mComputedEndFrame); }

    // Reset path
    void resetPath(std::string& rPath, std::string appendage = "") const;

    // Load cubemap texture. Returns texture handle
    GLuint createCubemapTexture(
        std::string filepathPosX,
        std::string filepathNegX,
        std::string filepathPosY,
        std::string filepathNegY,
        std::string filepathPosZ,
        std::string filepathNegZ) const;

    // Setup
    const float mCameraSmoothDuration = 1.5f;
    const float mSamplePointSize = 2.f;
    const float mPathPointSize = 3.f;
    const float mMinDrawingExtentOffset = -5.f;
    const float mMaxDrawingExtentOffset = 5.f;
    const float mCameraDefaultAlpha = 90.f;
    const float mCameraDefaultBeta = 45.f;
    const glm::vec3 mInternalAtomColor = glm::vec3(0.75f, 0.75f, 0.75f);
    const glm::vec3 mSurfaceAtomColor = glm::vec3(1.f, 0.25f, 0.f);
    const glm::vec3 mInternalValidationSampleColor = glm::vec3(1.f, 0.9f, 0.0f);
    const glm::vec3 mSurfaceValidationSampleColor = glm::vec3(0.f, 1.0f, 0.2f);
    const float mClippingPlaneMin = 0.f;
    const float mClippingPlaneMax = 200.f;
    const std::string mWindowTitle = "Surface Dynamics Visualization";
    const float mDepthDarkeningMaxEnd = 1000.f;
    const int mInitialWindowWidth = 1280;
    const int mInitialWindowHeight = 720;
    const float mOutlineWidth = 0.15f;
    const glm::vec4 mOutlineColor = glm::vec4(1.f, 1.f, 0.f, 0.9f);
    const glm::vec3 mPastPathColor = glm::vec3(1.f, 0.f, 0.f);
    const glm::vec3 mFuturePathColor = glm::vec3(0.f, 1.f, 0.f);
    const int mCameraAutoCenterSmoothFrameRadius = 10;
    const glm::vec3 mInternalHullSampleColor = glm::vec3(0.0f, 0.0f, 0.0f);
    const glm::vec3 mSurfaceHullSampleColor = glm::vec3(1.0f, 1.0f, 1.0f);
    const glm::vec3 mFallbackAtomColor = glm::vec3(0.7f, 0.7f, 0.7f);
    const float mAscensionColorOffsetAngle = 1.25f * glm::pi<float>();
    const float mSurfaceMarkPointSize = 5.f;
    const glm::vec3 mSelectionColor = glm::vec3(0.2f, 1.0f, 0.0f);
    const bool mFrameLogging = false;
    const std::string mNoComputedFrameMessage = "Frame was not computed.";

    // Controllable parameters
    bool mShowCameraWindow = true;
    bool mShowInformationWindow = true;
    bool mShowValidationWindow = false;
    bool mShowRenderingWindow = false;
    bool mShowVisualizationWindow = true;
    bool mShowComputationWindow = true;

    bool mRotateCamera = false;
    bool mMoveCamera = false;
    int mSelectedAtom = 0;
    bool mRenderWithProbeRadius = false;
    bool mUsePerspectiveCamera = false; // removed from GUI since both spheres
                                        // and cut of spheres is only correct
                                        // for orthographic projection
    bool mShowInternal = true;
    bool mShowSurface = true;
    float mComputationProbeRadius = 1.4f;
    int mCPUThreads = 8;
    int mSurfaceValidationAtomSampleCount = 20;
    bool mShowValidationSamples = true;
    float mClippingPlane = 0.f;
    int mSurfaceValidationSeed = 0;
    bool mShowAxesGizmo = false;
    bool mPlayAnimation = false;
    int mPlayAnimationRate = 15;
    bool mShowInternalSamples = true;
    bool mShowSurfaceSamples = true;
    int mStartFrame = 0;
    int mEndFrame = 0;
    int mComputationStartFrame = 0;
    int mComputationEndFrame = 0;
    bool mExtractLayers = true;
    bool mRepeatAnimation = false;
    int mSmoothAnimationRadius = 0;
    float mSmoothAnimationMaxDeviation = 5;
    float mDepthDarkeningStart = 100.f;
    float mDepthDarkeningEnd = 500.f;
    bool mShowAnalysisWindow = true;
    bool mShowPath = true;
    int mPathFrameRadius = 5; // radius of frames which are visualized
    int mPathSmoothRadius = 0; // radius of frames which are used for smoothing the path
    SurfaceRendering mSurfaceRendering = HULL;
    Background mBackground = SCIENTIFIC;
    int mHullSampleCount = 250; // sample count per atom
    bool mRenderHullSamples = false;
    bool mRenderOutline = true;
    bool mShowTooltips = true;
    float mAscensionUpToHotFrameCount = 10.f;
    float mAscensionBackToHotFrameCount = 10.f;
    float mAscensionUpToColdFrameCount = 10.f;
    float mAscensionBackToColdFrameCount = 10.f;
    bool mMarkSurfaceAtoms = false;
    bool mRepeatOnlyComputed = false;
    bool mSuperSampling = true;
    int mPathLengthStartFrame = 0;
    int mPathLengthEndFrame = 0;

    // Report output
    std::string mComputeInformation = "No computation info available.";
    std::string mValidationInformation = "No validation info available.";

    // Members
    GLFWwindow* mpWindow;
    std::unique_ptr<OrbitCamera> mupCamera; // camera for visualization
    glm::vec2 mCameraDeltaRotation;
    float mCameraRotationSmoothTime;
    glm::vec3 mLightDirection;
    std::unique_ptr<GPUProtein> mupGPUProtein; // protein on GPU
    std::unique_ptr<GPUSurfaceExtraction> mupGPUSurfaceExtraction;  // factory for GPUSurfaces
                                                                    // (unique pointer because has to be constructed after OpenGL initialization)
    std::vector<std::unique_ptr<GPUSurface> > mGPUSurfaces; // vector with surfaces
    int mFrame = 0; // do not set it directly, let it be done by setFrame() method!
    int mLayer = 0;
    float mFramePlayTime = 0; // time of displaying a molecule state at playing the animation
    int mComputedStartFrame = -1;
    int mComputedEndFrame = -1;
    float mComputedProbeRadius = 0.f;
    bool mFramebuffersExist = false;
    int mWindowWidth;
    int mWindowHeight;
    GLuint mScientificCubemapTexture;
    GLuint mCVCubemapTexture;
    GLuint mBeachCubemapTexture;
    std::unique_ptr<Framebuffer> mupMoleculeFramebuffer;
    std::unique_ptr<Framebuffer> mupSelectedAtomFramebuffer;
    std::unique_ptr<Framebuffer> mupOverlayFramebuffer;
    std::unique_ptr<GPUTextureBuffer> mupOutlineAtomIndices;
    std::unique_ptr<GPUBuffer<GLfloat> > mupAscension; // values have range [0..2Pi]
    std::unique_ptr<GPUHullSamples> mupHullSamples;

    // Analysis
    std::set<GLuint> mAnalyseGroup;
    int mNextAnalyseAtomIndex = 0;
    std::vector<float> mAnalysisSurfaceAmount;
    std::vector<float> mAnalysisSurfaceArea;
    std::vector<float> mAnalysisGroupMinLayers;
    std::vector<float> mAnalysisGroupAvgLayers;
    std::vector<float> mAnalysisGroupSurfaceAmount;
    std::vector<float> mAnalysisGroupSurfaceArea;
    std::unique_ptr<Path> mupPath;
    std::string mSurfaceIndicesFilePath = "";
    std::string mGlobalAnalysisFilePath = "";
    std::string mGroupAnalysisFilePath = "";

    // Surface validation
    std::unique_ptr<SurfaceValidation> mupSurfaceValidation;
    int mSurfaceValidationSampleCount = 0;
};

#endif // SURFACE_DYNAMICS_VISUALIZATION_H
