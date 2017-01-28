//============================================================================
// Distributed under the MIT License. Author: Adrian Derstroff
//============================================================================

#version 430

// Color of impostor
out int   atomID;
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
uniform float probeRadius;

// Main function
void main()
{
    /*
     * get atom center for every vertex id
     */
    int index = int(gl_VertexID);
    atomID = index;
    gl_Position = vec4(atoms[index].center, 1);

    /*
     * get the radius for the corresponding atom
     */
    vertRadius = atoms[index].radius + probeRadius;
}
