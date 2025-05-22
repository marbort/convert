#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <regex>

int main(int argc, char* argv[]) {
    if (argc < 2) {
        std::cerr << "Usage: " << argv[0] << " <filename>\n";
        return 1;
    }

    std::string filename = argv[1];
    std::ifstream infile(filename);
    if (!infile) {
        std::cerr << "Error opening file for reading: " << filename << "\n";
        return 1;
    }

    std::vector<std::string> lines;
    std::string line;
    while (std::getline(infile, line)) {
        lines.push_back(line + "\n");  // preserve newline
    }
    infile.close();

    std::regex pattern("^0\\n$");
    int begin = 0;

    for (size_t i = 0; i < lines.size(); ++i) {
        if (std::regex_match(lines[i], pattern)) {
            std::cout << "Found\n";
            begin = (i > 0) ? i - 1 : 0;
            break;  // Assuming you want only the first match
        }
    }

    std::ofstream outfile(filename, std::ios::trunc);
    if (!outfile) {
        std::cerr << "Error opening file for writing: " << filename << "\n";
        return 1;
    }

    for (size_t i = begin; i < lines.size(); ++i) {
        outfile << lines[i];
    }

    outfile.close();
    return 0;
}
