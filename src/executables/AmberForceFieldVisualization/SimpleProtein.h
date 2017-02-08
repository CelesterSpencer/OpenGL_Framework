//============================================================================
// Distributed under the MIT License. Author: Adrian Derstroff
//============================================================================

#ifndef OPENGL_FRAMEWORK_SIMPLEPROTEIN_H
#define OPENGL_FRAMEWORK_SIMPLEPROTEIN_H

#define FLOAT_MIN std::numeric_limits<float>::min()
#define FLOAT_MAX std::numeric_limits<float>::max()

// external includes
#include <GL/glew.h>
#include <glm/glm.hpp>

// project specific includes
#include "SimpleAtom.h"

class SimpleProtein {

public:
    //_____________________________________//
    //            VARIABLES                //
    //_____________________________________//
    std::string name;
    std::vector<SimpleAtom> atoms;
    glm::vec3 bbMin;
    glm::vec3 bbMax;
    glm::mat4 rotationMatrix;
    glm::mat4 translationMatrix;

    //_____________________________________//
    //           CONSTRUCTOR               //
    //_____________________________________//
    SimpleProtein(std::string name = "Protein") : name(name)
    {
        bbMin = glm::vec3(0, 0, 0);
        bbMax = glm::vec3(0, 0, 0);
        rotationMatrix = glm::mat4();
        translationMatrix = glm::mat4();
    }

    void setName(std::string name)
    {
        this->name = name;
    }

    int getNumberOfAtoms()
    {
        return atoms.size();
    }

    void move(glm::vec3 offset)
    {
        glm::mat4 translationMat = glm::translate(offset);
        translationMatrix = translationMat * translationMatrix;
        recalculateBB();
    }

    void rotate(glm::vec3 axis, float angle)
    {
        glm::mat4 rotationMat = glm::rotate(angle, axis);
        rotationMatrix = rotationMat * rotationMatrix;
        recalculateBB();
    }

    void applyRotationMatrix(glm::mat4 mat)
    {
        rotationMatrix = mat * rotationMatrix;
        recalculateBB();
    }

    void applyTranslationMatrix(glm::mat4 mat)
    {
        translationMatrix = mat * translationMatrix;
        recalculateBB();
    }

    glm::vec3 getCenterOfGravity()
    {
        glm::vec3 cog = glm::vec3(0,0,0);
        for (int i = 0; i < atoms.size(); i++) {
            cog += glm::vec3(translationMatrix * rotationMatrix * glm::vec4(atoms.at(i).pos, 1.0));
        }
        if (atoms.size() > 0) cog /= atoms.size();
        return cog;
    }

    void center()
    {
        glm::vec3 cog = getCenterOfGravity();
        move(cog);
        for (int i = 0; i < atoms.size(); i++) {
            atoms.at(i).pos -= cog;
        }
        recalculateBB();
    }

    void recalculateBB()
    {
        if (atoms.size() > 0) {
            bbMin = glm::vec3(FLOAT_MAX, FLOAT_MAX, FLOAT_MAX);
            bbMax = glm::vec3(FLOAT_MIN, FLOAT_MIN, FLOAT_MIN);
            for (int i = 0; i < atoms.size(); i++) {
                SimpleAtom atom = atoms.at(i);
                float radius = atom.radius;
                glm::vec4 transformedPos =  translationMatrix * rotationMatrix * glm::vec4(atom.pos, 1.0);
                bbMin.x = glm::min(bbMin.x, transformedPos.x-atom.radius);
                bbMin.y = glm::min(bbMin.y, transformedPos.y-atom.radius);
                bbMin.z = glm::min(bbMin.z, transformedPos.z-atom.radius);
                bbMax.x = glm::max(bbMax.x, transformedPos.x+atom.radius);
                bbMax.y = glm::max(bbMax.y, transformedPos.y+atom.radius);
                bbMax.z = glm::max(bbMax.z, transformedPos.z+atom.radius);
            }
        } else {
            Logger::instance().print("No atoms to recalculate BB!", Logger::Mode::WARNING);
        }
    }

    glm::vec3 extent()
    {
        //recalculateBB();
        return (bbMax - bbMin);
    }

};

#endif //OPENGL_FRAMEWORK_SIMPLEPROTEIN_H
