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

uniform int numAtoms;
uniform int numAtomSymbols;

uint getIndex(int x, int y, int z, int w) {
    return x + y*numAtomSymbols + z*numAtomSymbols*numAtomSymbols + w*numAtomSymbols*numAtomSymbols*numAtomSymbols;
}

void calcBond(
    in AtomStruct a1,
    in AtomStruct a2,
    in int isI,
    in int isJ,
    inout resDir)
{
    vec4 entry = valuesVector[indexCube[getIndex(a1.atomSymbolIndex, a2.atomSymbolIndex, 0, 0)]];
    float k = entry.x;
    float r0 = entry.y;

    vec3 r = a2.center - a1.center;

    vec3 f1 = 2*k*(length(r) - r0) * r / length(r);
    vec3 f2 = -f1;

    resDir += (isI*f1) + (isJ*f2);
}

void calcAngle(
    in AtomStruct a1,
    in AtomStruct a2,
    in AtomStruct a3,
    in int isI,
    in int isJ,
    in int isK,
    inout resDir)
{
    vec4 entry = valuesVector[indexCube[getIndex(a1.atomSymbolIndex, a2.atomSymbolIndex, a3.atomSymbolIndex, 0)]];
    float k = entry.x;
    float a0 = entry.y;

    vec3 v1 = a2.center - a1.center;
    vec3 v2 = a3.center - a2.center;
    float alpha = acos(dot(v1,v2) / (length(v1) * length(v2)));

    vec3 f1 = -2*k*(a-a0) / (length(v1) * sin(alpha)) * (v2 - (cos(alpha)*v1));
    vec3 f3 = -2*k*(a-a0) / (length(v2) * sin(alpha)) * (v1 - (cos(alpha)*v2));
    vec3 f2 = -f1 -f3;

    resDir += (isI*f1) + (isJ*f2) + (isK*f3);
}

void calcDihedral(
    in AtomStruct a1,
    in AtomStruct a2,
    in AtomStruct a3,
    in AtomStruct a4,
    in int isI,
    in int isJ,
    in int isK,
    in int isL,
    inout vec3 resDir)
{
    vec4 entry = valuesVector[indexCube[getIndex(a1.atomSymbolIndex, a2.atomSymbolIndex, a3.atomSymbolIndex, a4.atomSymbolIndex)]];
    float k = entry.x;
    float n = entry.w;
    float delta = entry.z;

    vec3 r1 = a2.center - a1.center;
    vec3 r2 = a3.center - a2.center;
    vec3 r3 = a4.center - a3.center;

    vec3 a = cross(r1,r2);
    vec3 b = cross(r2,r3);
    vec3 c = cross(r2,a);

    float cosa = dot(a,b) / (length(a) * length(b));
    float sina = dot(c,b) / (length(c) * length(b));
    float alpha = tan(sina, cosa);

    vec3 ds  = (b - (cosa * a)) / length(a);
    vec3 db  = (a - (cosa * b)) / length(b);
    float mu = -k * n * sin(n*alpha - delta) / sina;

    vec3 f1 = -mu * cross(r2,ds);
    vec3 f2 =  mu * (cross(r2,ds) - cross(db,r2));
    vec3 f3 =  mu * (cross(db,r2) - cross(ds,r1) - cross(r3,db));
    vec3 f4 =  mu * (cross(ds,r1) + cross(r3,db));

    resDir += (isI*f1) + (isJ*f2) + (isK*f3) + (isL*f4);
}

void calcImpropers(
    in AtomStruct a1,
    in AtomStruct a2,
    in AtomStruct a3,
    in AtomStruct a4,
    in int isI,
    in int isJ,
    in int isK,
    in int isL,
    inout vec3 resDir)
{
    vec4 entry = valuesVector[indexCube[getIndex(a1.atomSymbolIndex, a2.atomSymbolIndex, a3.atomSymbolIndex, a4.atomSymbolIndex)]];
    float s0 = entry.z;
    float k  = entry.y;
    float n  = entry.w;

    vec3 r1 = a2.center - a1.center;
    vec3 r2 = a3.center - a2.center;
    vec3 r3 = a4.center - a3.center;
    vec3 r4 = a4.center - a1.center;

    vec3 a  = cross(r1,r2);
    vec3 b  = cross(r2,r3);
    vec3 c  = cross(r2,a);
    vec3 d  = cross(r4,r1);

    float cosa = dot(a,b) / (length(a) * length(b));
    float sina = dot(c,b) / (length(c) * length(b));
    float alpha = atan(sina,cosa);

    vec3 ds  = (b - (cosa * a)) / length(a);
    vec3 db  = (a - (cosa * b)) / length(b);
    float s   = dot(a,d) / (length(a) * length(d));
    float mu = 2 * k * (s - s0);

    vec3 f1 = -mu * cross(r2,ds);
    vec3 f2 =  mu * (cross(r2,ds) - cross(db,r2));
    vec3 f3 =  mu * (cross(db,r2) - cross(ds,r1) - cross(r3,db));
    vec3 f4 =  mu * (cross(ds,r1) + cross(r3,db));

    resDir += (isI*f1) + (isJ*f2) + (isK*f3) + (isL*f4);
}

void calcVdW(
    in AtomStruct a1,
    in AtomStruct a2,
    in int isI,
    in int isJ,
    inout vec3 resDir)
{
    vec4 entry1 = valuesVector[indexCube[getIndex(0, a1.atomSymbolIndex, 0, 0)]];
    vec4 entry2 = valuesVector[indexCube[getIndex(0, a2.atomSymbolIndex, 0, 0)]];
    float sigma = entry1.x + entry2.x;
    float epsilon = sqrt(entry1.y * entry2.y);
    vec3 r = a2.center - a1.center;

    float scale = 48*epsilon*(pow(sigma/length(r)),12) - 0.5*pow(sigma/length(r),6) / (length(r)*length(r));
    vec3 f1 = scale*r;
    vec3 f2 = -f1;

    resDir += (isI*f1) + (isJ*f2);
}

void main() {
    // get element index
    uint i = gl_GlobalInvocationID.x;
    if (i > numAtoms) return;

    AtomStruct currentAtom = atoms[i];


}
