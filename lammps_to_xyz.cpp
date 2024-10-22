#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <iomanip> // For formatted output
#include <map>

// Structure to hold atom data
struct Atom {
    int id;
    int type;
    double x, y, z;
};

// Structure to hold box dimensions
struct Box {
    double xlo, xhi;
    double ylo, yhi;
    double zlo, zhi;
    double sidex, sidey, sidez;
};

// Function to process a single timestep and write it to an XYZ file
void writeXYZFrame(std::ofstream &xyzFile, const std::vector<Atom> &atoms, const Box &box, const std::map<int, std::string> &elementMap) {
    // Number of atoms
    xyzFile << atoms.size() << std::endl;

    // Comment line with box dimensions
    xyzFile << box.sidex << " " << box.sidey << " " << box.sidez << std::endl;

    // Write atom data in XYZ format
    for (const auto &atom : atoms) {
        // If elements are provided, use the element name; otherwise, use the type
        std::string element = elementMap.count(atom.type) ? elementMap.at(atom.type) : "AtomType" + std::to_string(atom.type);
        
        xyzFile << element << " " << std::fixed << std::setprecision(5) 
                << atom.x*(box.sidex) << " " << atom.y*(box.sidey) << " " << atom.z*(box.sidez) << std::endl;
    }
}

// Main function to convert LAMMPS dump to XYZ
int main(int argc, char *argv[]) {
    if (argc < 3) {
        std::cerr << "Usage: " << argv[0] << " <lammps_dump_file> <output_xyz_file> [element1 element2 ...]" << std::endl;
        return 1;
    }

    std::ifstream dumpFile(argv[1]);
    std::ofstream xyzFile(argv[2]);

    if (!dumpFile.is_open() || !xyzFile.is_open()) {
        std::cerr << "Error opening file!" << std::endl;
        return 1;
    }

    // Map atom types to elements (optional)
    std::map<int, std::string> elementMap;
    if (argc > 3) {
        for (int i = 3; i < argc; ++i) {
            elementMap[i - 2] = argv[i]; // map type 1 to element 1, type 2 to element 2, etc.
        }
    }

    std::string line;
    std::vector<Atom> atoms;
    Box box;
    bool readingAtoms = false;
    int numAtoms = 0;

    // Read the dump file
    while (std::getline(dumpFile, line)) {
        std::istringstream iss(line);

        if (line.find("ITEM: TIMESTEP") != std::string::npos) {
            // Ignore timestep line
            readingAtoms = false;
            std::getline(dumpFile, line);  // Skip the actual timestep number

        } else if (line.find("ITEM: NUMBER OF ATOMS") != std::string::npos) {
            std::getline(dumpFile, line);  // Get the number of atoms
            numAtoms = std::stoi(line);
            atoms.clear();
            atoms.reserve(numAtoms);

        } else if (line.find("ITEM: BOX BOUNDS") != std::string::npos) {
            // Read the box dimensions
            std::getline(dumpFile, line);
            std::istringstream(line) >> box.xlo >> box.xhi;
            std::getline(dumpFile, line);
            std::istringstream(line) >> box.ylo >> box.yhi;
            std::getline(dumpFile, line);
            std::istringstream(line) >> box.zlo >> box.zhi;

        } else if (line.find("ITEM: ATOMS") != std::string::npos) {
            readingAtoms = true;  // Start reading atoms in the following lines
        } else if (readingAtoms) {
            // Reading atoms' id, type, and positions
            Atom atom;
            iss >> atom.id >> atom.type >> atom.x >> atom.y >> atom.z;
            atoms.push_back(atom);

            if (atoms.size() == numAtoms) {
                // Finished reading atoms for this frame, write to XYZ
                box.sidex=box.xhi-box.xlo;
                box.sidey=box.yhi-box.ylo;
                box.sidez=box.zhi-box.zlo;
                writeXYZFrame(xyzFile, atoms, box, elementMap);
                readingAtoms = false;
            }
        }
    }

    dumpFile.close();
    xyzFile.close();

    std::cout << "Conversion completed!" << std::endl;
    return 0;
}
