//============================================================================
// Distributed under the MIT License. Author: Adrian Derstroff
//============================================================================

#ifndef OPENGL_FRAMEWORK_PROTEINLOADER_H
#define OPENGL_FRAMEWORK_PROTEINLOADER_H

#define FLOAT_MIN std::numeric_limits<float>::min()
#define FLOAT_MAX std::numeric_limits<float>::max()

// standard includes
#include <limits>
#include <map>
#include <algorithm>
#include <regex>

// framework includes
#include "Molecule/MDtrajLoader/MdTraj/MdTrajWrapper.h"
#include "Utils/Logger.h"

// project specific includes
#include "SimpleProtein.h"

class ProteinLoader {

public:
    //__________________PUBLIC__________________//
    //_____________________________________//
    //           CONSTRUCTOR               //
    //_____________________________________//
    ProteinLoader();
    ~ProteinLoader();

    int getNumberOfProteins();
    std::vector<SimpleProtein*> getProteins();
    SimpleProtein* getProteinAt(int i);
    std::vector<SimpleAtom> &getAllAtoms();
    std::vector<uint> &getAllNeighbors();
    void updateAtoms();
    int  getNumberOfAllAtoms();
    void getBoundingBoxAroundProteins(glm::vec3& min, glm::vec3& max);
    void getCenteredBoundingBoxAroundProteins(glm::vec3& min, glm::vec3& max);


    //_____________________________________//
    //             METHODS                 //
    //_____________________________________//
    SimpleProtein* loadProtein(std::string fileName, std::map<std::string, uint> atomSymbolsMap);
    void loadPDB(std::string filePath, SimpleProtein &protein, std::map<std::string, uint> atomSymbolsMap, glm::vec3 &minPosition, glm::vec3 &maxPosition);
private:
    std::vector<SimpleProtein*> m_proteins;
    std::vector<SimpleAtom> m_allAtoms;
    std::vector<uint> m_allNeighbors;
    float m_currentProteinIdx;

    // tokenize input based on the given delimiter (regex)
    void tokenize(std::vector<std::string> & tokens, const std::string & in, const std::string & delimiter)
    {
        std::regex regex(delimiter);
        std::sregex_token_iterator iter(in.begin(), in.end(), regex, -1);
        std::sregex_token_iterator end;
        tokens = std::vector<std::string>(iter, end);
    }
};


#endif //OPENGL_FRAMEWORK_PROTEINLOADER_H
