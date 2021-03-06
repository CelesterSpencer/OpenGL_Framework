//============================================================================
// Distributed under the MIT License. Author: Adrian Derstroff
//============================================================================

#ifndef OPENGL_FRAMEWORK_NEIGHBORHOODSEARCH_H
#define OPENGL_FRAMEWORK_NEIGHBORHOODSEARCH_H

#include <GL/glew.h>
#include <math.h>
#include <malloc.h>
#include <string.h>
#include <ShaderTools/ShaderProgram.h>
#include <Utils/Timer.h>

#include "GPUHandler.h"
#include "AssertData.h"
#include "Utils/Logger.h"
#include "NeighborhoodSearchDefines.h"
//#include "../../executables/NeighborSearchTest/SimpleProtein.h"


class NeighborhoodSearch {
public:
    ~NeighborhoodSearch();

    /*
     * Getter and setter
     */
    int getNumberOfGridCells();
    glm::vec3 getGridSize();
    glm::ivec3 getGridResolution();
    float getCellSize();
    void getGridMinMax(glm::vec3& min, glm::vec3& max);
    int getGridSearch();
    int getNumberOfBlocksForElementsComputation();
    int getNumberOfThreadsPerBlockForElementsComputation();
    int getNumberOfBlocksForGridComputation();
    int getNumberOfThreadsPerBlockForGridComputation();
    float getMaxSearchRadius();
    int getTotalGridNum();

    /*
     * neighbor search
     */
    void init(uint numElements, glm::fvec3 min, glm::fvec3 max, glm::ivec3 resolution, float searchRadius);
    void update(uint numElements, glm::fvec3 min, glm::fvec3 max, glm::ivec3 resolution, float searchRadius);
    void run(GLuint* positionsSSBO, Neighborhood& neighborhood);


private:
    // grid parameters
    uint*       m_grid;
    uint*       m_gridCnt;
    int         m_gridSearch;
    int         m_gridAdj[216];     // maximal size of the adjacency mask is 6x6x6
    int         m_gridAdjCnt;       // 3D search count =n^3 e.g. 2x2x2=8
    int         m_gridAdjOff;       // adjacency mask offset of the upper left cell
    glm::fvec3  m_gridMin;
    glm::fvec3  m_gridMax;
    glm::ivec3  m_gridRes;          // 3D grid resolution
    glm::fvec3  m_gridSize;         // 3D grid sizes
    glm::fvec3  m_gridDelta;        // delta translate from world space to cell space
    int         m_gridTotal;        // total number of cells in the grid
    Grid        m_gridDataGPU;
    float       m_cellSize;
    float       m_searchRadius;
    float       m_maxSearchRadius;
    int         m_numElements;      // number of particles

    // blocksums parameters
    uint        m_numLevelsAllocated;
    GLuint**    m_scanBlockSumsInt;

    // gpu
    GPUBuffers    m_gpuBuffers;
    uint          m_numBlocks;
    uint          m_numThreads;
    uint          m_gridBlocks;
    uint          m_gridThreads;

    // compute shader
    ShaderProgram m_insertElementsShader;
    ShaderProgram m_prescanIntShader;
    ShaderProgram m_uniformAddIntShader;
    ShaderProgram m_fillTempDataShader;
    ShaderProgram m_countingSortShader;



    /*
     * init helper functions
     */
    void setupComputeShaders();
    void allocateBuffers(uint numElements);
    void deallocateBuffers();
    void preallocBlockSumsInt(uint maxNumElements);
    void deallocBlockSumsInt();
    void setupGrid(glm::fvec3 min, glm::fvec3 max, glm::ivec3 resolution, float searchRadius);
    void freeGrid();
    void calculateNumberOfBlocksAndThreads(uint numElements);
    void computeNumBlocks(int numElements, int maxThreads, uint& numBlocks, uint &numThreads);

    /*
     * run helper functions
     */
    void insertElementsInGridGPU();

    void prefixSumCellsGPU();
    void prescanArrayRecursiveInt(GLuint* outArray, GLuint* inArray, int numElements, int level);
    void prescanInt(int numThreads, int numBlocks, int sharedMemSize, bool storeSum, bool isNP2, GLuint* outArray, GLuint* inArray, int level, int n, int blockIndex, int baseIndex);
    void uniformAddInt(int numThreads, int numBlocks, GLuint* outArray, int level, int n, int blockOffset, int baseIndex);
    bool isPowerOfTwo(int n);
    int floorPow2(int n);

    /*
     * sorting
     */
    void countingSort();
};


#endif //OPENGL_FRAMEWORK_NEIGHBORHOODSEARCH_H
