//============================================================================
// Distributed under the MIT License. Author: Adrian Derstroff
//============================================================================

#version 430

// Color of impostor
out vec4 vertColor;
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
layout(std430, binding = 10) restrict readonly buffer ModelMatrixBuffer { mat4 modelMatrix[]; };

// Uniforms
uniform float probeRadius;


// Main function
void main()
{
    /*
     * get atom center for every vertex id
     */
    int index = int(gl_VertexID);
    AtomStruct atom = atoms[index];
    gl_Position = modelMatrix[int(atom.proteinID)] * vec4(atom.center, 1);

    /*
     * get the radius and charge for the corresponding atom
     */
    vertRadius = atom.radius + probeRadius;
    vertColor = vec4(normalize(atoms[index].charge.xyz), atoms[index].charge.w);
}
