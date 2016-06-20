#ifndef CPP_IMPLEMENTATION_H
#define CPP_IMPLEMENTATION_H

// Notes
// - Face is defined by vec4(Normal, Distance from origin)

#include "SurfaceExtraction/GPUProtein.h"
#include <glm/glm.hpp>
#include <vector>

class CPPImplementation
{
public:

    void execute(
        int executionIndex,
        int atomCount,
        float probeRadius,
        const std::vector<GPUAtom>& atoms,
        std::vector<unsigned int>& internalIndices,
        std::vector<unsigned int>& surfaceIndices);

private:

    void setup();

    bool checkParallelism(
        glm::vec4 plane,
        glm::vec4 otherPlane) const;

    bool pointInHalfspaceOfPlane(
        glm::vec4 plane,
        glm::vec3 point) const;

    void intersectPlanes(
        glm::vec4 plane,
        glm::vec4 otherPlane,
        glm::vec3 &linePoint,
        glm::vec3 &lineDir) const;

    float underSQRT(
        glm::vec3 linePoint,
        glm::vec3 lineDir,
        glm::vec3 sphereCenter,
        float sphereRadius) const;

    bool testEndpoint(
        glm::vec3 endpoint,
        int excludeA,
        int excludeB) const;

    // Members
    static const int neighborsMaxCount = 200;
    static const bool logging = false; // one has to remove /* */ before activating logging

    // All cutting faces, also those who gets cut away by others
    int cuttingFaceCount = 0;
    glm::vec3 cuttingFaceCenters[neighborsMaxCount];
    glm::vec4 cuttingFaces[neighborsMaxCount]; // Normal + Distance

    // Selection of cutting faces which get intersected pairwaise and produce endpoints
    int cuttingFaceIndicators[neighborsMaxCount]; // Indicator whether cutting face was cut away by other (1 == not cut away)
    int cuttingFaceIndicesCount = 0; // Count of not cut away cutting faces
    int cuttingFaceIndices[neighborsMaxCount]; // Indices of cutting faces which are not cut away by other
};

#endif // CPP_IMPLEMENTATION_H
