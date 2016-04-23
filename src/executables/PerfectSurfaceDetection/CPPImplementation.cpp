#include "CPPImplementation.h"

#include <iostream>

void CPPImplementation::setup()
{
    cuttingFaceCount = 0;
    cuttingFaceIndicesCount = 0;
}

// ## Intersection line of two planes
// http://stackoverflow.com/questions/6408670/line-of-intersection-between-two-planes
void CPPImplementation::intersectPlanes(
    glm::vec3 faceCenter,
    glm::vec3 faceNormal,
    glm::vec3 otherFaceCenter,
    glm::vec3 otherFaceNormal,
    glm::vec3 &linePoint,
    glm::vec3 &lineDir) const
{
    // Direction of line (length already should be one)
    lineDir = cross(faceNormal, otherFaceNormal);

    // No determinant necessary, because check for parallelity already done and normals are unit vectors

    // Point on line
    linePoint =
        (cross(lineDir, otherFaceNormal) * length(faceCenter)
        + (cross(faceNormal, lineDir) * length(otherFaceCenter)));
}

// ## Part under square root of intersection line and sphere
// https://en.wikipedia.org/wiki/Line%E2%80%93sphere_intersection
float CPPImplementation::underSQRT(
    glm::vec3 linePoint,
    glm::vec3 lineDir,
    glm::vec3 atomCenter,
    float atomRadius) const
{
    float underSQRT1 = glm::dot(lineDir, (linePoint - atomCenter));
    underSQRT1 = underSQRT1 * underSQRT1;
    float underSQRT2 = glm::length(linePoint - atomCenter);
    underSQRT2 = underSQRT2 * underSQRT2;
    return (underSQRT1 - underSQRT2 + (atomRadius * atomRadius));
}

// ## Function to test whether endpoint is NOT cut away. Called after cutting face list is minimized
bool CPPImplementation::testEndpoint(glm::vec3 endpoint) const
{
    // Just some reservation
    glm::vec3 vector;

    // Iterate over cuttingFaceIndices entries
    for(int i = 0; i < cuttingFaceIndicesCount; i++)
    {
        // Index of cutting face
        int index = cuttingFaceIndices[i];

        // Vector from cutting face center to endpoint
        vector = endpoint - cuttingFaceCenters[index];

        // Dot product should have positive sign if in same half space
        if(glm::dot(vector, cuttingFaceNormals[index]) > 0)
        {
            return false;
        }
    }
    return true;
}

// ## Execution function
void CPPImplementation::execute(
    int executionIndex,
    int atomCount,
    float probeRadius,
    const std::vector<AtomStruct>& atoms,
    std::vector<unsigned int>& surfaceAtomsIndices)
{
    std::cout << "Execution for atom: " << executionIndex << std::endl;

    // Reset members for new execution
    setup();

    // Index
    int atomIndex = executionIndex;

    // Check whether in range
    if(atomIndex >= atomCount) { return; }

    // When no endpoint was generated at all, atom is surface (value is false then)
    bool endpointGenerated = false;

    // When one endpoint survives cutting, atom is surface (value is true then)
    bool endpointSurvivesCut = false;

    // ### OWN VALUES ###

    // Own extended radius
    float atomExtRadius = atoms[atomIndex].radius + probeRadius;

    // Own center
    glm::vec3 atomCenter = atoms[atomIndex].center;

    // ### BUILD UP OF CUTTING FACE LIST ###

    // Go over other atoms and build cutting face list
    for(int i = 0; i < atomCount; i++)
    {
        // Do not cut with itself
        if(i == atomIndex) { continue; }

        // ### OTHER'S VALUES ###

        // Get values from other atom
        glm::vec3 otherAtomCenter = atoms[i].center;
        float otherAtomExtRadius = atoms[i].radius + probeRadius;

        // ### INTERSECTION TEST ###

        // Vector from center to other's
        glm::vec3 connection = otherAtomCenter - atomCenter;

        // Distance between atoms
        float atomDistance = glm::length(connection);

        // Do they intersect with extended radii?
        if(atomDistance >= (atomExtRadius + otherAtomExtRadius)) { continue; }

        // NODE: Following cases are NOT considered
        // - both spheres touch each other in one point
        // - one sphere lies completely inside other
        // - one sphere lies completely inside other and touches other's surface in one point

        // ### INTERSECTION ###

        // Squared atom distance
        float atomDistanceSquared = atomDistance * atomDistance;

        // Calculate center of intersection
        // http://gamedev.stackexchange.com/questions/75756/sphere-sphere-intersection-and-circle-sphere-intersection
        float h =
            0.5
            + ((atomExtRadius * atomExtRadius)
            - (otherAtomExtRadius * otherAtomExtRadius))
            / atomDistanceSquared;

        // ### CUTTING FACE LIST ###

        // Use connection between centers as line
        cuttingFaceCenters[cuttingFaceCount] = atomCenter + h * connection;

        // Calculate radius of intersection
        cuttingFaceRadii[cuttingFaceCount] =
            sqrt((atomExtRadius * atomExtRadius)
            - (h * h * atomDistanceSquared));

        // Calculate normal of intersection
        cuttingFaceNormals[cuttingFaceCount] = normalize(connection);

        // Initialize cutting face indicator with: 1 == was not cut away (yet)
        cuttingFaceIndicators[cuttingFaceCount] = 1;

        // Increment cutting face list index and break if max count of neighbors reached
        if((++cuttingFaceCount) == neighborsMaxCount) { break; }
    }

    // ### GET RID OF CUTTING FACES WHICH ARE COMPLETELY CUT AWAY BY OTHER ###

    for(int i = 0; i < cuttingFaceCount - 1; i++)
    {
        // Already cut away (think about where this is save to test)
        // if(cuttingFaceIndicators[i] == 0) { continue; }

        // Values of cutting face
        glm::vec3 faceCenter = cuttingFaceCenters[i];
        float faceRadius = cuttingFaceRadii[i];
        glm::vec3 faceNormal = cuttingFaceNormals[i];

        // Test every cutting face for intersection line with other
        for(int j = i+1; j < cuttingFaceCount; j++)
        {
            // Already cut away (think about where this is save to test)
            // if(cuttingFaceIndicators[j] == 0) { continue; }

            // Values of other cutting face
            glm::vec3 otherFaceCenter = cuttingFaceCenters[j];
            float otherFaceRadius = cuttingFaceRadii[j];
            glm::vec3 otherFaceNormal = cuttingFaceNormals[j];

            // ### CHECK THAT FACES TO NOT CUT EACH OTHER WITHIN ATOM ###

            // Check for parallelity, first
            bool notCutEachOther = (1.0 == glm::abs(dot(faceNormal, otherFaceNormal))); // If already parallel, they do not cut

            // Do further checking if not already parallel
            if(!notCutEachOther)
            {
                // Intersection of planes, resulting in line
                glm::vec3 lineDir; glm::vec3 linePoint;
                intersectPlanes(
                    faceCenter,
                    faceNormal,
                    otherFaceCenter,
                    otherFaceNormal,
                    linePoint,
                    lineDir);

                // Intersection of line with sphere, resulting in two, one or no endpoints
                // https://en.wikipedia.org/wiki/Line%E2%80%93sphere_intersection
                float valueUnderSQRT = underSQRT(linePoint, lineDir, atomCenter, atomExtRadius);

                // Only interesting case is for zero endpoints, because then there is no cut on atom surface
                notCutEachOther = valueUnderSQRT < 0;
            }

            // ### CHECK WHETHER CUTTING FACE CAN BE FORGOT ###

            // Faces do not cut each other, so they produce not later endpoints. Check them now
            if(notCutEachOther)
            {
                // Connection between faces' center
                glm::vec3 connection = otherFaceCenter - faceCenter;

                // Booleans for check
                bool faceConnectionPositive = dot(faceNormal, connection) > 0;
                bool otherFaceConnectionPositive = dot(otherFaceNormal, connection) > 0;

                // Atom is completely cut away by faces
                if(faceConnectionPositive && !otherFaceConnectionPositive)
                {
                    // No surface atom since cut away by the two cutting faces
                    return;
                }

                // Check for direction of normals (bigger zero is enough, since there is not intersection of both)
                if(dot(faceNormal, otherFaceNormal) > 0)
                {
                    if(!faceConnectionPositive && !otherFaceConnectionPositive)
                    {
                        // Face is cut away by other face
                        cuttingFaceIndicators[i] = 0;

                    }
                    else if(faceConnectionPositive && otherFaceConnectionPositive)
                    {
                        // Other face is cut away by face
                        cuttingFaceIndicators[j] = 0;
                    }
                }
            }
        }
    }

    // ### GO OVER CUTTING FACES AND COLLECT NOT CUT AWAY ONES ###

    for(int i = 0; i < cuttingFaceCount; i++)
    {
        // Check whether cutting face is still there after preprocessing
        if(cuttingFaceIndicators[i] == 1)
        {
            // Save index of that cutting face
            cuttingFaceIndices[cuttingFaceIndicesCount] = i;

            // Increase count of those cutting faces
            cuttingFaceIndicesCount++;
        }
    }

    // ### GO OVER OK CUTTING FACES AND TEST ENDPOINTS ###

    for(int i = 0; i < cuttingFaceIndicesCount; i++)
    {
        // Values of cutting face
        int index = cuttingFaceIndices[i];
        glm::vec3 faceCenter = cuttingFaceCenters[index];
        float faceRadius = cuttingFaceRadii[index];
        glm::vec3 faceNormal = cuttingFaceNormals[index];

        // Test every cutting face for intersection line with other
        for(int j = i+1; j < cuttingFaceIndicesCount; j++)
        {
            // Values of other cutting face
            int otherIndex = cuttingFaceIndices[j];
            glm::vec3 otherFaceCenter = cuttingFaceCenters[otherIndex];
            float otherFaceRadius = cuttingFaceRadii[otherIndex];
            glm::vec3 otherFaceNormal = cuttingFaceNormals[otherIndex];

            // Intersection of planes, resulting in line
            glm::vec3 lineDir; glm::vec3 linePoint;
            intersectPlanes(
                faceCenter,
                faceNormal,
                otherFaceCenter,
                otherFaceNormal,
                linePoint,
                lineDir);

            // Intersection of line with sphere, resulting in two, one or no endpoints
            // https://en.wikipedia.org/wiki/Line%E2%80%93sphere_intersection
            float valueUnderSQRT = underSQRT(linePoint, lineDir, atomCenter, atomExtRadius);

            // Left part of equation
            float left = -(glm::dot(lineDir, (linePoint - atomCenter)));

            // Check value under square root
            if(valueUnderSQRT > 0)
            {
                // Some endpoint was at least generated
                endpointGenerated = true;

                // Right part of equation
                float right = sqrt(valueUnderSQRT);

                // First endpoint
                float d = left + right;
                if(testEndpoint(linePoint + d * lineDir))
                {
                    // Break out of for loop (and outer)
                    endpointSurvivesCut = true;
                    break;
                }

                // Second endpoint
                d = left - right;
                if(testEndpoint(linePoint + d * lineDir))
                {
                    // Break out of for loop (and outer)
                    endpointSurvivesCut = true;
                    break;
                }
            }
            else if(valueUnderSQRT == 0)
            {
                // Some endpoint was at least generated
                endpointGenerated = true;

                // Just test the one endpoint
                float d = left;
                if(testEndpoint(linePoint + d * lineDir))
                {
                    // Break out of for loop (and outer)
                    endpointSurvivesCut = true;
                    break;
                }
            }
            // else, no endpoint is generated since cutting faces do not intersect

            if(endpointSurvivesCut) { break; }
        }
    }

    // ### ATOM IS SURFACE ATOM ###

    // If no endpoint was generated at all or one or more survived cutting, add this atom to surface
    if(!endpointGenerated || endpointSurvivesCut)
    {
        surfaceAtomsIndices.push_back((unsigned int) atomIndex);
    }
}