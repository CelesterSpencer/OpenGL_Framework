//
// Created by ubundrian on 29.11.16.
//

#ifndef MOLECULARDYNAMICSVISUALIZATION_AMBERFORCEFIELDPARAMETER_H
#define MOLECULARDYNAMICSVISUALIZATION_AMBERFORCEFIELDPARAMETER_H

#include <vector>
#include <glm/glm.hpp>

typedef unsigned int uint;



class AmberForceFieldParameter
{
private:
    std::vector<uint> m_indexCube;
    std::vector<glm::vec4> m_valuesList;
    std::map<std::string, uint> m_atomSymbolsMap;
    uint m_size;

    uint getCubeIndex (uint x, uint y, uint z, uint w) { return x + y*m_size + z*m_size*m_size + w*m_size*m_size*m_size; }

    void updateCubeIndexY(uint index, uint y)
    {
        m_indexCube.at(getCubeIndex(0,y,0,0)) = index;
    }
    void updateCubeIndex(uint index, uint x)
    {
        m_indexCube.at(getCubeIndex(x,0,0,0)) = index;
    }
    void updateCubeIndex(uint index, uint x, uint y)
    {
        m_indexCube.at(getCubeIndex(x,y,0,0)) = index;
        m_indexCube.at(getCubeIndex(y,x,0,0)) = index;
    }
    void updateCubeIndex(uint index, uint x, uint y, uint z)
    {
        m_indexCube.at(getCubeIndex(x,y,z,0)) = index;
        m_indexCube.at(getCubeIndex(x,z,y,0)) = index;
        m_indexCube.at(getCubeIndex(y,x,z,0)) = index;
        m_indexCube.at(getCubeIndex(y,z,x,0)) = index;
        m_indexCube.at(getCubeIndex(z,x,y,0)) = index;
        m_indexCube.at(getCubeIndex(z,y,x,0)) = index;
    }
    void updateCubeIndex(uint index, uint x, uint y, uint z, uint w)
    {
        m_indexCube.at(getCubeIndex(x,y,z,w)) = index;
        m_indexCube.at(getCubeIndex(x,y,w,z)) = index;
        m_indexCube.at(getCubeIndex(x,z,y,w)) = index;
        m_indexCube.at(getCubeIndex(x,z,w,y)) = index;
        m_indexCube.at(getCubeIndex(x,w,y,z)) = index;
        m_indexCube.at(getCubeIndex(x,w,z,y)) = index;

        m_indexCube.at(getCubeIndex(y,x,z,w)) = index;
        m_indexCube.at(getCubeIndex(y,x,w,z)) = index;
        m_indexCube.at(getCubeIndex(y,z,x,w)) = index;
        m_indexCube.at(getCubeIndex(y,z,w,x)) = index;
        m_indexCube.at(getCubeIndex(y,w,x,z)) = index;
        m_indexCube.at(getCubeIndex(y,w,z,x)) = index;

        m_indexCube.at(getCubeIndex(z,x,y,w)) = index;
        m_indexCube.at(getCubeIndex(z,x,w,y)) = index;
        m_indexCube.at(getCubeIndex(z,y,x,w)) = index;
        m_indexCube.at(getCubeIndex(z,y,w,x)) = index;
        m_indexCube.at(getCubeIndex(z,w,x,y)) = index;
        m_indexCube.at(getCubeIndex(z,w,y,x)) = index;

        m_indexCube.at(getCubeIndex(w,x,y,z)) = index;
        m_indexCube.at(getCubeIndex(w,x,z,y)) = index;
        m_indexCube.at(getCubeIndex(w,y,x,z)) = index;
        m_indexCube.at(getCubeIndex(w,y,z,x)) = index;
        m_indexCube.at(getCubeIndex(w,z,x,y)) = index;
        m_indexCube.at(getCubeIndex(w,z,y,x)) = index;
    }
public:
    AmberForceFieldParameter()
    {
        m_size = 0;
    }
    AmberForceFieldParameter(std::map<std::string, uint>& atomSymbolsMap, uint size)
    {
        uint cubeSize = (size+1)*(size+1)*(size+1)*(size+1); // (#atoms+1)^4
        m_indexCube = std::vector<uint>(cubeSize,0);
        m_size = size;
        m_atomSymbolsMap = atomSymbolsMap;
    }

    /**
     * inserts entry into vector while keeping track of all
     * combinations for all indices xyzw
     * if there was no previous entry return -1 else return 1
     */
    int set(glm::vec4 entry, uint x, uint y, uint z, uint w)
    {
        // determine number of 0 indexes
        int numberOfZeros = 0;
        if (x == 0) numberOfZeros++;
        if (y == 0) numberOfZeros++;
        if (z == 0) numberOfZeros++;
        if (w == 0) numberOfZeros++;

        // get the index inside the values array for the given 4D index (xyzw)
        uint cubeIndex = getCubeIndex(x,y,z,w);
        uint valuesListIndex = m_indexCube.at(cubeIndex);

        // if array index is zero it means that entry hasn't been set yet
        if (valuesListIndex == 0)
        {
            m_valuesList.push_back(entry);
            uint insertedIndex = m_valuesList.size(); // caution actual index is (inserted index - 1) since its 1 based

            // assign inserted index to all possible cells within the cube
            switch(numberOfZeros)
            {
                case 1:
                    if (x != 0)
                        updateCubeIndex(insertedIndex,x);
                    else
                        updateCubeIndexY(insertedIndex,y);
                    break;
                case 2:
                    updateCubeIndex(insertedIndex,x,y);
                    break;
                case 3:
                    updateCubeIndex(insertedIndex,x,y,z);
                    break;
                case 4:
                    updateCubeIndex(insertedIndex,x,y,z,w);
                    break;
            }

            return -1;
        }
            // else get the index and change the value in the array
        else
        {
            m_valuesList.at(valuesListIndex-1) = entry;
        }
    }
    std::vector<uint> &getIndexCube() { return m_indexCube; }
    std::vector<glm::vec4> &getValuesList() { return m_valuesList; }
    std::map<std::string, uint> &getAtomSymbolsMap() { return m_atomSymbolsMap; };
    uint getNumberOfAtomSymbols() { return m_size; }
};

#endif //MOLECULARDYNAMICSVISUALIZATION_AMBERFORCEFIELDPARAMETER_H
