#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <cmath>
#include <iomanip>
#include <map>

// --- Project-137: Universe OS v9.3 "The Sculptor" ---
// コンセプト: シートの幾何学に「宇宙のねじれ（1/137）」を再統合し、完全な3D構造を復元する

const double PI = 3.1415926535;
const double GOLDEN_ANGLE = PI * (3.0 - std::sqrt(5.0)); // 137.508...度
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
    std::cout << "Project-137 v9.3 [The Sculptor 3D]" << std::endl;
    std::cout << "Input: ";
    std::cin >> inputSeq;

    std::vector<Residue> protein;
    double currX = 0, currY = 0, currZ = 0;
    double zDir = 1.0, xOff = 0.0;

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
        // 宇宙のねじれ角
        double phi = i * GOLDEN_ANGLE;

        if (c1 != -1 && i >= c1 && i <= c6) {
            // --- 【3D Ring Mode】 ---
            double progress = (double)(i - c1) / (c6 - c1);
            double theta = progress * 2.0 * PI;
            dx = std::cos(theta) * 1.5;
            dy = std::sin(theta) * 1.5;
            dz = std::cos(phi) * 0.8; // Z方向にも宇宙のゆらぎを追加
        } 
        else if (res.code == 'P' || (i > 0 && res.code == 'G' && inputSeq[i-1] == 'G')) {
            // --- 【Sheet Turn Mode】 ---
            zDir *= -1.0;
            xOff += 5.0; // 隣の列へ
            dx = std::cos(phi) * 2.0;
            dy = std::sin(phi) * 2.0;
            dz = zDir * 2.0;
        } 
        else {
            // --- 【3D Zigzag Mode】 ---
            // 黄金角による回転を加えながらジグザグに進む
            double zigzag = (i % 2 == 0 ? 1.0 : -1.0);
            dx = std::cos(phi) * zigzag + std::sin(phi) * 0.5;
            dy = std::sin(phi) * zigzag + std::cos(phi) * 0.5;
            dz = zDir * 2.5; 
        }

        double len = std::sqrt(dx*dx + dy*dy + dz*dz);
        currX += (dx/len) * CA_DIST;
        currY += (dy/len) * CA_DIST;
        currZ += (dz/len) * CA_DIST;

        res.x = currX + xOff; res.y = currY; res.z = currZ;
        protein.push_back(res);
    }

    std::ofstream out("Universe_v9_3_Sculptor.pdb");
    out << "HEADER    PROJECT-137 V9.3 3D SCULPTOR\n";
    for (auto& r : protein) {
        char buf[100];
        sprintf(buf, "ATOM  %5d  CA  %3s A%4d    %8.3f%8.3f%8.3f  1.00  4.00           C\n",
                r.resSeq, r.resName.c_str(), r.resSeq, r.x, r.y, r.z);
        out << buf;
    }
    if(c1 != -1) out << "CONECT " << std::setw(5) << c1 + 1 << std::setw(5) << c6 + 1 << "\n";
    out << "END\n";
    
    std::cout << "[Success] 3D Sculpting complete. v9.3 generated." << std::endl;
    return 0;
}

