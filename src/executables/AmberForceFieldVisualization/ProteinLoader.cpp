//============================================================================
// Distributed under the MIT License. Author: Adrian Derstroff
//============================================================================

#include "ProteinLoader.h"

ProteinLoader::ProteinLoader()
{
    m_currentProteinIdx = 0;
}

ProteinLoader::~ProteinLoader()
{
    for (int i = 0; i < m_proteins.size(); i++) {
        delete m_proteins.at(i);
    }
    m_proteins.clear();
}

int ProteinLoader::getNumberOfProteins()
{
    return m_proteins.size();
}

std::vector<SimpleProtein*> ProteinLoader::getProteins()
{
    return m_proteins;
}

SimpleProtein* ProteinLoader::getProteinAt(int i)
{
    return m_proteins.at(i);
}

void ProteinLoader::updateAtoms()
{
    m_allAtoms.clear();
    for (int i = 0; i < getNumberOfProteins(); i++) {
        SimpleProtein* protein = getProteinAt(i);
        m_allAtoms.insert(std::end(m_allAtoms), std::begin(protein->atoms), std::end(protein->atoms));
    }
}

std::vector<SimpleAtom> &ProteinLoader::getAllAtoms()
{
    updateAtoms();
    return m_allAtoms;
}

int ProteinLoader::getNumberOfAllAtoms()
{
    return m_allAtoms.size();
}

std::vector<uint> &ProteinLoader::getAllNeighbors()
{
    return m_allNeighbors;
}

/*
 * LOADING PROTEIN
 */
SimpleProtein* ProteinLoader::loadProtein(std::string fileName, std::map<std::string, uint> atomSymbolsMap)
{
    /*
     * extracting protein name from file name
     */
    std::string proteinName = fileName;
    int lastSlash = proteinName.find_last_of("/");
    if (lastSlash >= 0) {
        proteinName = proteinName.substr(lastSlash+1, proteinName.size());
    }
    int lastDot = proteinName.find_last_of(".");
    if (lastDot >= 0) {
        proteinName = proteinName.substr(0, lastDot);
    }



    /*
     * concatenate full path
     */
    std::string subfolder = "/molecules/";
    std::string filePath = RESOURCES_PATH + subfolder + fileName;



    /*
     * load protein from pdb file
     */
    SimpleProtein* protein = new SimpleProtein;
    protein->name = proteinName;
    loadPDB(filePath, *protein, atomSymbolsMap, protein->bbMin, protein->bbMax);
    m_proteins.push_back(protein);

    return m_proteins.at(m_proteins.size()-1);
}

void ProteinLoader::loadPDB(std::string filePath, SimpleProtein &protein, std::map<std::string, uint> atomSymbolsMap, glm::vec3 &minPosition, glm::vec3 &maxPosition)
{

    /*
     * set start values of min and max position
     */
    double min = std::numeric_limits<double>::min();
    double max = std::numeric_limits<double>::max();
    minPosition = glm::vec3(max, max, max);
    maxPosition = glm::vec3(min, min, min);
    float maxRadius = 0;

    /*
     * the informations we need to create our atoms
     */
    std::vector<std::string> residueNames;
    std::vector<std::string> elementNames;
    std::vector<std::string> elementIds;
    std::vector<glm::vec3> positions;
    std::vector<float> radii;
    std::vector<std::string> bonds;
    std::map<std::string, uint> elementOrderMap;
    std::map<std::string, std::vector<std::string> > elementNeighborsMap;


    /*
     * init mdtraj wrapper and load pdb file as pyobject
     */
    MdTrajWrapper mdTrajWrapper;
    PyObject* pdbFile = mdTrajWrapper.loadFilePDB(filePath);


    /*
     * Atom positions
     */
    PyObject* xyz_py = mdTrajWrapper.getXYZ(pdbFile);
    PyObject* topo = mdTrajWrapper.getTopology(pdbFile);
    Py_DECREF(pdbFile);
    PyArrayObject* xyz_pyarray = reinterpret_cast<PyArrayObject*>(xyz_py);

    //get number of frames, atoms and components
    long long numAtom{ PyArray_SHAPE(xyz_pyarray)[1] };
    long long numComponents{ PyArray_SHAPE(xyz_pyarray)[2] };
    int numAtoms = (int)numAtom;

    // cast the 3D numpy array to a 1D c array
    float* xyz_carray;
    xyz_carray = reinterpret_cast<float*>(PyArray_DATA(xyz_pyarray));

    float positionX;
    float positionY;
    float positionZ;

    std::vector<glm::vec3> frameHolder;
    for (int a = 0; a < numAtoms; a++)
    {
        int id = a * numComponents;
        positionX = xyz_carray[id] * 10;
        positionY = xyz_carray[id + 1] * 10;
        positionZ = xyz_carray[id + 2] * 10;
        glm::vec3 position(positionX, positionY, positionZ);
        positions.push_back(position);

        /*
         * get min and max
         */
        minPosition.x = (positionX < minPosition.x) ? positionX : minPosition.x;
        minPosition.y = (positionY < minPosition.y) ? positionY : minPosition.y;
        minPosition.z = (positionZ < minPosition.z) ? positionZ : minPosition.z;

        maxPosition.x = (positionX > maxPosition.x) ? positionX : maxPosition.x;
        maxPosition.y = (positionY > maxPosition.y) ? positionY : maxPosition.y;
        maxPosition.z = (positionZ > maxPosition.z) ? positionZ : maxPosition.z;
    }
    Py_DECREF(xyz_py);


    /*
     * Atom properties
     */
    PyObject* atom_py = PyObject_GetAttrString(topo, "atoms");
    PyObject* atom_iterator_py = PyObject_GetIter(atom_py);
    Py_DECREF(atom_py);

    PyObject* atom;
    PyObject* residue_py;
    PyObject* residue_name_py;
    PyObject* element_py;
    PyObject* element_name_py;
    PyObject* element_specific_name_py;
    PyObject* atom_radius_py;

    atom = PyIter_Next(atom_iterator_py);
    while ((atom != NULL))
    {
        residue_py = PyObject_GetAttrString(atom, "residue");
        residue_name_py = PyObject_Str(residue_py);
        element_py = PyObject_GetAttrString(atom, "element");
        element_name_py = PyObject_GetAttrString(element_py, "symbol");
        element_specific_name_py = PyObject_GetAttrString(atom, "name");
        atom_radius_py = PyObject_GetAttrString(element_py, "radius");

        elementNames.push_back(PyUnicode_AsUTF8(element_name_py));
        std::string residueName = PyUnicode_AsUTF8(residue_name_py);
        std::string elementSpecificName = PyUnicode_AsUTF8(element_specific_name_py);
        std::string elementId = residueName + "-" + elementSpecificName;
        elementIds.push_back(elementId);
        residueNames.push_back(residueName);
        double radius = PyFloat_AsDouble(atom_radius_py) * 10;
        radii.push_back(radius);

        /*
         * find the atom with the biggest radius
         */
        maxRadius = std::max((float)radius, maxRadius);

        /*
         * create map between element id and element position
         */
        elementOrderMap[elementId] = elementIds.size()-1;

        Py_DECREF(residue_name_py);
        Py_DECREF(element_py);
        Py_DECREF(element_name_py);
        Py_DECREF(element_specific_name_py);
        Py_DECREF(atom_radius_py);

        atom = PyIter_Next(atom_iterator_py);
    }
    Py_DECREF(atom_iterator_py);


    /*
     * atom bonds
     */
    PyObject* bonds_py = PyObject_GetAttrString(topo, "bonds");
    Py_DECREF(topo);
    PyObject* bonds_iterator_py = PyObject_GetIter(bonds_py);

    PyObject* bond;

    bond = PyIter_Next(bonds_iterator_py);
    while ((bond != NULL))
    {
        PyObject* bondfi = PyObject_Str(bond);
        bonds.push_back(PyUnicode_AsUTF8(bondfi));
        Py_DECREF(bondfi);
        bond = PyIter_Next(bonds_iterator_py);

        /*
         * insert Bond into map
         * both combinations need to be added
         * since bonds have no direction
         */
        std::vector<std::string> tokens;
        tokenize(tokens, bonds.at(bonds.size()-1), "\\(|\\)|,\\s");
        if (elementNeighborsMap.find(tokens.at(1)) == elementNeighborsMap.end()) {
            elementNeighborsMap[tokens.at(1)] = std::vector<std::string>();
        }
        elementNeighborsMap.at(tokens.at(1)).push_back(tokens.at(2));
        if (elementNeighborsMap.find(tokens.at(2)) == elementNeighborsMap.end()) {
            elementNeighborsMap[tokens.at(2)] = std::vector<std::string>();
        }
        elementNeighborsMap.at(tokens.at(2)).push_back(tokens.at(1));
    }
    Py_DECREF(bonds_iterator_py);


    /*
     * extent the bounding box by the radius of the biggest atom
     */
    minPosition -= maxRadius;
    maxPosition += maxRadius;


    /*
     * Create atoms and add them to the protein
     * and to the allAtoms vector
     */
    if ((positions.size() == residueNames.size()) && (positions.size() == elementNames.size()) && (positions.size() == radii.size()))
    {
        uint atomsVectorStart = m_allAtoms.size();
        for (int i = 0; i < residueNames.size(); i++)
        {
            std::string elementName = elementNames.at(i);
            glm::vec3 position = positions.at(i);
            float radius = radii.at(i);

            // check if residueName is in map
            std::transform(elementName.begin(), elementName.end(), elementName.begin(), ::tolower); // turn element residueName to lowercase
            if (atomSymbolsMap.find(elementName) == atomSymbolsMap.end())
            {
                Logger::instance().print("Could not found element " + elementName + " in map!", Logger::Mode::ERROR);
            }

            /*
             * create atom
             */
            SimpleAtom atom;
            atom.pos = position;
            atom.radius = radius;
            atom.charge = glm::vec4(0.0,0.0,0.0,0.0);
            atom.atomSymbolIndex = atomSymbolsMap[elementName]+1;                   // +1 because the first atom should start with 1
            atom.proteinID = m_currentProteinIdx;
            //std::cout << elementName << ": " << atom.atomSymbolIndex << std::endl;

            /*
             * get all neighbors
             */
            std::string elementId = elementIds.at(i);
            atom.bondNeighborsStart = m_allNeighbors.size();
            if (elementNeighborsMap.find(elementId) != elementNeighborsMap.end())
            {
                std::vector<std::string> neighborIds = elementNeighborsMap.at(elementId);
                atom.bondNeighborsSize = neighborIds.size();
                for (int i = 0; i < neighborIds.size(); i++)
                {
                    /*
                     * atomsVectorStart is necessary if this is not the first protein that is
                     * going to be importet, so that the right atom positions within the
                     * allAtoms vector are referenced.
                     */
                    m_allNeighbors.push_back(atomsVectorStart + elementOrderMap.at(neighborIds.at(i)));
                }
            }
            else
            {
                atom.bondNeighborsSize = 0;
                Logger::instance().print(elementId + " has no bonds!");
            }

            /*
             * add atom to both protein and all atoms
             */
            protein.atoms.push_back(atom);
            m_allAtoms.push_back(atom);
        }
    } else {
        std::cerr << "Size of atom positions " << positions.size() << " and properties " << residueNames.size() << ", " << elementNames.size() << ", " << radii.size() << " dont match" << std::endl;
    }


    /*
     * increment protein idx
     */
    m_currentProteinIdx++;
}

void ProteinLoader::getBoundingBoxAroundProteins(glm::vec3& min, glm::vec3& max)
{
    if (m_proteins.size() > 0)
    {
        min = glm::vec3(FLOAT_MAX, FLOAT_MAX, FLOAT_MAX);
        max = glm::vec3(FLOAT_MIN, FLOAT_MIN, FLOAT_MIN);
        for (int i = 0; i < m_proteins.size(); i++)
        {
            SimpleProtein* protein = m_proteins.at(i);
            min.x = glm::min(min.x, protein->bbMin.x);
            min.y = glm::min(min.y, protein->bbMin.y);
            min.z = glm::min(min.z, protein->bbMin.z);
            max.x = glm::max(max.x, protein->bbMax.x);
            max.y = glm::max(max.y, protein->bbMax.y);
            max.z = glm::max(max.z, protein->bbMax.z);
        }
    }
    else
    {
        Logger::instance().print("No proteins there to calculate bounding box!", Logger::Mode::WARNING);
        min = glm::vec3(0,0,0);
        max = glm::vec3(0,0,0);
    }
}

void ProteinLoader::getCenteredBoundingBoxAroundProteins(glm::vec3& min, glm::vec3& max)
{
    getBoundingBoxAroundProteins(min, max);
    glm::vec3 extent = max - min;
    glm::vec3 center = (max + min) / 2;
    float longestSideHalf = std::max(std::max(extent.x, extent.y), extent.z) / 2;
    min = center - longestSideHalf;
    max = center + longestSideHalf;
}
