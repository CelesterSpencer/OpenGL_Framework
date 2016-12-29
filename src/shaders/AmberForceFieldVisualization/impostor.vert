//============================================================================
// Distributed under the MIT License. Author: Adrian Derstroff
//============================================================================

#version 430

// Color of impostor
out vec3 vertColor;
out float vertRadius;

// Struct for atom
struct AtomStruct
{
    vec3 center;
    float radius;
    float atomSymbolIndex;
    float bondNeighborsStart;
    float bondNeighborsSize;
    float proteinID;
    vec4 charge;
};

// SSBOs
layout(std430, binding = 0) restrict readonly buffer AtomBuffer { AtomStruct atoms[]; };

// Uniforms
uniform vec3 cameraWorldPos;
uniform float probeRadius;
uniform int proteinNum;

// Main function
void main()
{
    /*
     * get atom center for every vertex id
     */
    int index = int(gl_VertexID);
    gl_Position = vec4(atoms[index].center, 1);

    /*
     * get the radius for the corresponding atom
     */
    vertRadius = atoms[index].radius + probeRadius;

    /*
     * color selected atom different to all atoms

    float proteinIdx = atoms[index].proteinID;
    float PI = 3.1415926;
    float colorF = PI*(float(proteinIdx)/proteinNum);
    float sinC = sin(colorF);
    float cosC = cos(colorF);
    vec3 proteinColor = vec3(sinC, cosC, min(1,sinC+cosC));
    vertColor = proteinColor;
    */

    //vertColor = normalize(atoms[index].charge.xyz) * 255;
    vertColor = normalize(atoms[index].charge.xyz);
}
