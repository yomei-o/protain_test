#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <cmath>
#include <iomanip>
#include <map>

// --- Project-137: Universe OS v8.2 "Aqueous Drag" ---
// コンセプト: 「水の抵抗」によって、直進エネルギーを強制的に螺旋へと変換する

const double PI = 3.1415926535;
const double GOLDEN_ANGLE = PI * (3.0 - std::sqrt(5.0)); // 137.508...
const double CA_DIST = 3.8; 

struct Residue {
    int resSeq; char code; std::string resName;
    double x, y, z;
};

void normalize(double &x, double &y, double &z) {
    double len = std::sqrt(x*x + y*y + z*z);
    if (len > 0) { x /= len; y /= len; z /= len; }
}

std::pair<std::string, int> getAAInfo(char code) {
    static std::map<char, std::pair<std::string, int>> m = {
        {'A', {"ALA", 1}}, {'R', {"ARG", 6}}, {'N', {"ASN", 3}}, {'C', {"CYS", 1}},
        {'Q', {"GLN", 4}}, {'G', {"GLY", 0}}, {'L', {"LEU", 3}}, {'P', {"PRO", 2}}
    }; // (省略)
    return m.count(code) ? m[code] : std::make_pair("UNK", 0);
}

int main() {
    std::string inputSeq;
    std::cout << "Project-137 [v8.2 Aqueous Drag Mode]\nInput: ";
    std::cin >> inputSeq;

    std::vector<Residue> protein;
    double currX = 0, currY = 0, currZ = 0;
    
    // 物理パラメータ
    double angleX = 0, angleY = 0;
    double dragCoeff = 0.05; // 水の抵抗係数

    for (int i = 0; i < (int)inputSeq.length(); ++i) {
        Residue res;
        res.resSeq = i + 1;
        res.code = inputSeq[i];
        res.resName = getAAInfo(res.code).first;

        // --- 【Aqueous Drag Engine】 ---
        // 1. 直進しようとするほど、黄金角方向への「回転トルク」が強まる
        // 2. 水の抵抗(dragCoeff)が、直進軸(Z)を少しずつ横へ押し出す
        
        angleX += GOLDEN_ANGLE; 
        angleY += 0.2 + (dragCoeff * i * 0.01); // 進むほど抵抗によるブレが増す

        // 抵抗による「しなり」をベクトルに変換
        double dirX = std::sin(angleX) * std::cos(angleY);
        double dirY = std::sin(angleX) * std::sin(angleY);
        double dirZ = std::cos(angleX) * (1.0 - dragCoeff); // Z軸への伸びを抑制

        normalize(dirX, dirY, dirZ);
        
        currX += dirX * CA_DIST;
        currY += dirY * CA_DIST;
        currZ += dirZ * CA_DIST;

        res.x = currX; res.y = currY; res.z = currZ;
        protein.push_back(res);
    }

    std::ofstream out("Universe_v8_2_Drag.pdb");
    out << "HEADER    PROJECT-137 V8.2 AQUEOUS DRAG\n";
    for (auto& r : protein) {
        char buf[100];
        sprintf(buf, "ATOM  %5d  CA  %3s A%4d    %8.3f%8.3f%8.3f  1.00  4.00           C\n",
                r.resSeq, r.resName.c_str(), r.resSeq, r.x, r.y, r.z);
        out << buf;
    }
    out << "END\n";
    std::cout << "[Success] v8.2 Aqueous Drag generated.\n";
    return 0;
}

