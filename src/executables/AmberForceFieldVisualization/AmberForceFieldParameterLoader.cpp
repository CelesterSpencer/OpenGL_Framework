//
// Created by ubundrian on 29.11.16.
//

#include <Utils/Logger.h>
#include "AmberForceFieldParameterLoader.h"



void AmberForceFieldParameterLoader::load(std::string & filename, AmberForceFieldParameter & amberParameters)
{
    std::ifstream file;
    file.open(filename);

    if (file.is_open())
    {
        std::string line;
        std::vector<std::string> tokens;

        // atom symbols map
        std::map<std::string, uint> atomSymbolsMap;
        uint atomSymbolIndex = 0;

        // skip the first line
        getline(file, line);

        /*
         * collect atom masses and polarities
         */
        std::vector<std::string> atomSymbols;
        std::vector<glm::vec4> atomValues;
        while (getline(file, line) && line.length() > 0)
        {
            // tokenize the input
            tokenize(tokens, line, "\\s+");

            // parse properties
            glm::vec4 atomProperties;
            atomSymbols.push_back(tokens[0]);           // atom symbol
            atomProperties.x = std::stod(tokens[1]);    // atom mass
            atomProperties.y = std::stod(tokens[2]);    // atom polarizability
            atomProperties.z = 0;
            atomProperties.w = 0;
            atomValues.push_back(atomProperties);

            // fill atom symbol map
            atomSymbolsMap[tokens[0]] = atomSymbolIndex++;
            Logger::instance().print(tokens[0] + ": " + std::to_string(atomSymbolsMap[tokens[0]]));
        }

        // skip list of atom symbols that are hydrophilic in solution
        getline(file, line);

        // initialize amber Parameters with number of atom symbols as side length of the 4D cube
        int cubeSide = atomSymbolIndex; // cubeSize is number of atom symbols+1 because cube is 1-based
        Logger::instance().print("Cube size: " + std::to_string(cubeSide));
        amberParameters = AmberForceFieldParameter(atomSymbolsMap, cubeSide);

        // set all atom properties
        for (uint i = 0; i < atomValues.size(); i++)
        {
            // insert properties into 4D cube
            uint x = atomSymbolsMap[atomSymbols.at(i)]+1; // x+1 because cube is 1-based
            uint status = amberParameters.set(atomValues.at(i), x , 0, 0, 0);
            if (status == 1)
                std::cerr << "Wrong insertion, values array index is non null! For x=" << x << std::endl;
        }

        /*
         * collect atom bonds
         */
        while (getline(file, line) && line.length() > 0)
        {
            // tokenize the input
            tokenize(tokens, line, "(\\s?-|\\s+)");

            // parse properties
            glm::vec4 bondProperties;
            std::string atomA = tokens[0];
            std::string atomB = tokens[1];

            // parse properties
            bondProperties.x = std::stod(tokens[2]); // bond force constant
            bondProperties.y = std::stod(tokens[3]);  // bond equilibrium
            bondProperties.z = 0;
            bondProperties.w = 0;

            // insert properties into 4D cube
            uint x = atomSymbolsMap[atomA]+1;
            uint y = atomSymbolsMap[atomB]+1;
            uint status = amberParameters.set(bondProperties, x, y, 0, 0);
            if (status == 1)
                std::cerr << "Wrong insertion, values array index is non null! For x=" << x << ", y=" << y << std::endl;
        }

        /*
         * collect atom angles
         */
        while (getline(file, line) && line.length() > 0)
        {
            // tokenize the input
            tokenize(tokens, line, "(\\s?-|\\s+|\\t)");

            // parse properties
            glm::vec4 angleProperties;
            std::string atomA = tokens[0];
            std::string atomB = tokens[1];
            std::string atomC = tokens[2];

            // parse properties
            angleProperties.x = std::stod(tokens[3]);   // angle force constant
            angleProperties.y = std::stod(tokens[4]);   // angle equilibrium
            angleProperties.z = 0;
            angleProperties.w = 0;

            // insert properties into 4D cube
            uint x = atomSymbolsMap[atomA]+1;
            uint y = atomSymbolsMap[atomB]+1;
            uint z = atomSymbolsMap[atomC]+1;
            uint status = amberParameters.set(angleProperties, x, y, z, 0);
            if (status == 1)
                std::cerr << "Wrong insertion, values array index is non null! For x=" << x << ", y=" << y << ", z=" << z << std::endl;
        }

        /*
         * collect dihedral angles
         */
        while (getline(file, line) && line.length() > 0)
        {
            // tokenize the input
            tokenize(tokens, line, "(\\s\\s+)");

            // parse properties
            glm::vec4 dihedralProperties;
            dihedralProperties.x = std::stod(tokens[1]);    // torsional barrier factor
            dihedralProperties.y = std::stod(tokens[2]);    // barrier height divided by 2
            dihedralProperties.z = std::stod(tokens[3]);    // dihedral angle
            dihedralProperties.w = std::stod(tokens[4]);    // torsional barrier periodicity

            // tokenize atom symbols
            tokenize(tokens, tokens[0], "(\\s?-)");
            std::string atomA = tokens[0];
            std::string atomB = tokens[1];
            std::string atomC = tokens[2];
            std::string atomD = tokens[3];

            // insert properties into 4D cube
            uint x = atomSymbolsMap[atomA]+1;
            uint y = atomSymbolsMap[atomB]+1;
            uint z = atomSymbolsMap[atomC]+1;
            uint w = atomSymbolsMap[atomD]+1;
            uint status = amberParameters.set(dihedralProperties, x, y, z, w);
            if (status == 1)
                std::cerr << "Wrong insertion, values array index is non null! For x=" << x << ", y=" << y << ", z=" << z << ", w=" << w << std::endl;
        }

        /*
         * collect impropers dihedral angles
         */
        while (getline(file, line) && line.length() > 0)
        {
            // tokenize the input
            tokenize(tokens, line, "(\\s\\s+)");

            // parse properties
            glm::vec4 improperDihedralProperties;
            improperDihedralProperties.x = 0;
            improperDihedralProperties.y = std::stod(tokens[1]);    // barrier height divided by 2
            improperDihedralProperties.z = std::stod(tokens[2]);    // dihedral angle
            improperDihedralProperties.w = std::stod(tokens[3]);    // torsional barrier periodicity

            // parse properties
            tokenize(tokens, tokens[0], "(\\s?-)");
            std::string atomA = tokens[0];
            std::string atomB = tokens[1];
            std::string atomC = tokens[2];
            std::string atomD = tokens[3];

            // insert properties into 4D cube
            uint x = atomSymbolsMap[atomA]+1;
            uint y = atomSymbolsMap[atomB]+1;
            uint z = atomSymbolsMap[atomC]+1;
            uint w = atomSymbolsMap[atomD]+1;
            uint status = amberParameters.set(improperDihedralProperties, x, y, z, w);
            if (status == 1)
                std::cerr << "Wrong insertion, values array index is non null! For x=" << x << ", y=" << y << ", z=" << z << ", w=" << w << std::endl;
        }

        /*
         * collect 10-12 hydrogen-bonds
         */
        while (getline(file, line) && line.length() > 0)
        {
            // trim leading whitespaces
            trim(line);

            // tokenize the input
            tokenize(tokens, line, "\\s\\s+");

            // get properties
            std::string atomA = tokens[0];
            std::string atomB = tokens[1];

            glm::vec4 hBondPotential;
            hBondPotential.x = std::stod(tokens[2]);  // 12th coefficient
            hBondPotential.y = std::stod(tokens[3]);  // 10th coefficient

            // not sure if it is necessary
            // TODO: CHECK and maybe remove the lines above
        }

        /*
         * collect equivalent atoms for vdW potential
         */
        while (getline(file, line) && line.length() > 0)
        {
            // not used
        }

        // skip one line (reason?)
        getline(file, line);

        /*
         * collect van der Waals parameters
         */
        while (getline(file, line) && line.length() > 0)
        {
            // trim leading whitespace
            trim(line);

            // tokenize the input
            tokenize(tokens, line, "\\s\\s+");

            std::string atomSymbol = tokens[0];

            // parse properties
            glm::vec4 vdWPotential;
            vdWPotential.x = std::stod(tokens[1]);  // van der Waals radius
            vdWPotential.y = std::stod(tokens[2]); // Lennard Jones potential well depth (epsilon)
            vdWPotential.z = 0;
            vdWPotential.w = 0;

            uint y = atomSymbolsMap[atomSymbol]+1;
            uint status = amberParameters.set(vdWPotential, 0, y, 0, 0);
            if (status == 1)
                std::cerr << "Wrong insertion, values array index is non null! For y=" << y << std::endl;
        }

        file.close();
    }
    else {
        std::cerr << "Error opening " << filename << std::endl;
    }
}