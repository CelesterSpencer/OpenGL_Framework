//============================================================================
// Distributed under the MIT License. Author: Adrian Derstroff
//============================================================================

#ifndef OPENGL_FRAMEWORK_SIMPLEATOM_H
#define OPENGL_FRAMEWORK_SIMPLEATOM_H

#include <glm/glm.hpp>

class SimpleAtom {
public:
    //_____________________________________//
    //            VARIABLES                //
    //_____________________________________//
    glm::vec3 pos;
    float radius;
    float atomSymbolIndex;
    float bondNeighborsStart;
    float bondNeighborsSize;
    float proteinID;
    glm::vec4 charge;
};

#endif //OPENGL_FRAMEWORK_SIMPLEATOM_H
