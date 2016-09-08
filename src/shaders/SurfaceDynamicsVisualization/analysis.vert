//============================================================================
// Distributed under the MIT License. Author: Raphael Menges
//============================================================================

#version 430

// Color and radius of impostor
out vec3 vertColor;
out float vertRadius;

// Index of atom
out int vertIndex;

// Radii
layout(std430, binding = 0) restrict readonly buffer RadiiBuffer
{
   float radii[];
};

// Trajectory
struct Position
{
    float x,y,z;
};

layout(std430, binding = 1) restrict readonly buffer TrajectoryBuffer
{
   Position trajectory[];
};

// Indices of group atoms
layout(binding = 2, r32ui) readonly restrict uniform uimageBuffer GroupAtomsIndices;

// Uniforms
uniform float probeRadius;
uniform int selectedIndex = 0;
uniform vec3 hotColor;
uniform vec3 coldColor;
uniform vec3 internalColor;
uniform int frame;
uniform int atomCount;
uniform int smoothAnimationRadius;
uniform float smoothAnimationMaxDeviation;
uniform int frameCount;
uniform int groupAtomCount;
uniform vec3 selectionColor;

// Global
int atomIndex;
vec3 centerAtFrame;

// Accumulate center at frame
void accumulateCenter(
    int accFrame,
    inout vec3 accCenters,
    inout int accCount)
{
    // Extract center at that frame
    Position position = trajectory[(accFrame*atomCount) + atomIndex];
    vec3 center = vec3(position.x, position.y, position.z);

    // Check whether center is not too far away
    float distanceToFramesCenter = distance(centerAtFrame, center);
    if(distanceToFramesCenter < smoothAnimationMaxDeviation)
    {
        accCenters += center;
        accCount++;
    }
}

// Main function
void main()
{
    // Extract center at frame which is given. Unlike hull shader, here are atom indices directly given by vertex id
    atomIndex = int(gl_VertexID); // write it to global variable
    Position position = trajectory[(frame*atomCount) + atomIndex];
    centerAtFrame = vec3(position.x, position.y, position.z); // write it to global variable

    // Calculate loop bounds for smoothing
    int lowerBound = max(0, frame - smoothAnimationRadius);
    int upperBound = min(frameCount - 1, frame + smoothAnimationRadius);
    int accCount = 0;
    vec3 accCenters = vec3(0,0,0);

    // Accumulate centers in frames below current one
    for(int i = lowerBound; i < frame; i++)
    {
        accumulateCenter(i, accCenters, accCount);
    }

    // Accumulate centers in frames above current one
    for(int i = frame + 1; i <= upperBound; i++)
    {
        accumulateCenter(i, accCenters, accCount);
    }

    // Extract center from accumulated ones
    vec3 center = (accCenters + centerAtFrame) / (accCount + 1);
    gl_Position = vec4(center, 1);

    // Extract radius
    vertRadius = radii[atomIndex] + probeRadius;

    // Determine whether this atom is in group
    bool inGroup = false;
    for(int i = 0; i < groupAtomCount; i++)
    {
        if(int(imageLoad(GroupAtomsIndices, i).x) == atomIndex)
        {
            inGroup = true;
            break;
        }
    }

    // Set color
    if(atomIndex == selectedIndex)
    {
        vertColor = selectionColor;
    }
    else if(inGroup)
    {
        vertColor = vec3(1, 1, 0);
    }
    else
    {
        vertColor = vec3(0.2, 0.2, 0.2);
    }

    // Set index
    vertIndex = atomIndex + 1; // plus one to distinguish from nothing!
}