#version 430

layout(local_size_x = 256, local_size_y = 1, local_size_z = 1) in;

// Struct for atom
struct AtomStruct
{
    vec3  center;
    float radius;
    float atomSymbolIndex;
    float bondNeighborsStart;
    float bondNeighborsSize;
    float charge;
    vec4  proteinID;
};

layout(std430, binding = 0) buffer AtomsBuffer      { AtomStruct atoms[];  };
layout(std430, binding = 1) buffer IndexCubeBuffer  { uint indexCube[];    };
layout(std430, binding = 2) buffer ValuesBuffer     { vec4 valuesVector[]; };

void main() {
    // get element index
    uint i = gl_GlobalInvocationID.x;
    AtomStruct currentAtom = atoms[i];


}
