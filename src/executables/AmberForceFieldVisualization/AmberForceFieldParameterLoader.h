//
// Created by ubundrian on 29.11.16.
//

#ifndef MOLECULARDYNAMICSVISUALIZATION_AMBERFORCEFIELDPARAMETERLOADER_H
#define MOLECULARDYNAMICSVISUALIZATION_AMBERFORCEFIELDPARAMETERLOADER_H

#include <stdlib.h>
#include <string>
#include <iostream>
#include <fstream>
#include <regex>
#include <limits>
#include <glm/glm.hpp>

#include "AmberForceFieldParameter.h"

typedef unsigned int uint;



class AmberForceFieldParameterLoader {
private:
    // trim from both ends
    void trim(std::string &s)
    {
        s.erase(s.begin(), std::find_if(s.begin(), s.end(),std::not1(std::ptr_fun<int, int>(std::isspace))));
        s.erase(std::find_if(s.rbegin(), s.rend(),std::not1(std::ptr_fun<int, int>(std::isspace))).base(), s.end());
    }

    // tokenize input based on the given delimiter (regex)
    void tokenize(std::vector<std::string> & tokens, const std::string & in, const std::string & delimiter)
    {
        std::regex regex(delimiter);
        std::sregex_token_iterator iter(in.begin(), in.end(), regex, -1);
        std::sregex_token_iterator end;
        tokens = std::vector<std::string>(iter, end);
    }

public:
    void load(std::string & filename, AmberForceFieldParameter & amberParameters);
};


#endif //MOLECULARDYNAMICSVISUALIZATION_AMBERFORCEFIELDPARAMETERLOADER_H
