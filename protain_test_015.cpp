#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <cmath>
#include <iomanip>
#include <map>

// --- Project-137: Universe OS v9.2 "The Architect" ---
// 修正：リング（輪）とシート（板）の共存。構造の優先順位を定義する。

const double PI = 3.1415926535;
const double CA_DIST = 3.8; 

struct Residue {
    int resSeq;
    char code;
    std::string resName;
    double x, y, z;
};

std::pair<std::string, int> getAAInfo(char code) {
    static std::map<char, std::pair<std::string, int>> m = {
        {'G', {"GLY", 0}}, {'A', {"ALA", 1}}, {'P', {"PRO", 2}}, {'V', {"VAL", 2}},
        {'C', {"CYS", 1}}, {'I', {"ILE", 3}}, {'L', {"LEU", 3}}, {'Y', {"TYR", 7}},
        {'Q', {"GLN", 4}}, {'N', {"ASN", 3}}, {'S', {"SER", 1}}, {'T', {"THR", 2}}
    };
    return m.count(code) ? m[code] : std::make_pair("UNK", 0);
}

int main() {
    std::string inputSeq;
    std::cout << "Project-137 v9.2 [The Architect Hybrid]" << std::endl;
    std::cout << "Input: ";
    std::cin >> inputSeq;

    std::vector<Residue> protein;
    double currX = 0, currY = 0, currZ = 0;
    double zDir = 1.0, xOff = 0.0, angle = 0.0;

    // オキシトシン環（C...C）の事前検知
    int c1 = -1, c6 = -1;
    for(int i=0; i < (int)inputSeq.length() - 5; i++) {
        if(inputSeq[i] == 'C' && inputSeq[i+5] == 'C') { c1 = i; c6 = i + 5; break; }
    }

    for (int i = 0; i < (int)inputSeq.length(); ++i) {
        Residue res;
        res.resSeq = i + 1;
        res.code = inputSeq[i];
        res.resName = getAAInfo(res.code).first;

        double dx, dy, dz;

        if (c1 != -1 && i >= c1 && i <= c6) {
            // --- 【Ring Mode】優先 ---
            double progress = (double)(i - c1) / (c6 - c1);
            double theta = progress * PI;
            dx = std::cos(theta);
            dy = std::sin(theta);
            dz = 0.3; 
        } 
        else if (res.code == 'P' || (i > 0 && res.code == 'G' && inputSeq[i-1] == 'G')) {
            // --- 【Sheet Turn Mode】 ---
            zDir *= -1.0;
            xOff += 4.8;
            dx = (i % 2 == 0 ? 0.5 : -0.5);
            dy = 0.1;
            dz = zDir;
        } 
        else {
            // --- 【Standard Zigzag Mode】 ---
            dx = (i % 2 == 0 ? 0.8 : -0.8);
            dy = 0.2;
            dz = zDir * 1.5;
        }

        double len = std::sqrt(dx*dx + dy*dy + dz*dz);
        currX += (dx/len) * CA_DIST;
        currY += (dy/len) * CA_DIST;
        currZ += (dz/len) * CA_DIST;

        res.x = currX + xOff; res.y = currY; res.z = currZ;
        protein.push_back(res);
    }

    std::ofstream out("Universe_v9_2_Architect.pdb");
    out << "HEADER    PROJECT-137 V9.2 HYBRID ARCHITECT\n";
    for (auto& r : protein) {
        char buf[100];
        sprintf(buf, "ATOM  %5d  CA  %3s A%4d    %8.3f%8.3f%8.3f  1.00  4.00           C\n",
                r.resSeq, r.resName.c_str(), r.resSeq, r.x, r.y, r.z);
        out << buf;
    }
    if(c1 != -1) out << "CONECT " << std::setw(5) << c1 + 1 << std::setw(5) << c6 + 1 << "\n";
    out << "END\n";
    
    std::cout << "[Complete] Architect engine finished." << std::endl;
    return 0;
}
