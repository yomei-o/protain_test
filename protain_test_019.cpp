#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <cmath>

// --- Project-137: v8.3 "Ribosome Extrusion" ---
// コンセプト: 合成時に「押し出される」際の、一定方向への微小な曲がりが螺旋を作る

const double PI = 3.1415926535;
const double CA_DIST = 3.8; 

// おじさんの理論：押し出される時にかかる「一定のクセ」
const double BEND_X = 0.15; // 毎回少しずつ「うなだれる」角度
const double TWIST_Y = 0.55; // 毎回少しずつ「回転する」角度

struct Residue {
    int resSeq; std::string resName;
    double x, y, z;
};

int main() {
    std::string inputSeq;
    std::cout << "Project-137 [v8.3 Ribosome Extrusion]\nSequence: ";
    std::cin >> inputSeq;

    std::vector<Residue> protein;
    double currX = 0, currY = 0, currZ = 0;
    
    // 押し出し方向の初期値（最初は真上に押し出される）
    double angleX = 0;
    double angleY = 0;

    for (int i = 0; i < (int)inputSeq.length(); ++i) {
        Residue res;
        res.resSeq = i + 1;
        res.resName = "ALA"; // 簡易化

        // 【おじさんの理論】
        // 押し出すたびに、一定の抵抗を受けて「同じ方向」に少しずつ曲がる
        angleX += BEND_X;
        angleY += TWIST_Y;

        // 向いている方向に 3.8 オングストローム押し出す
        double dirX = std::sin(angleX) * std::cos(angleY);
        double dirY = std::sin(angleX) * std::sin(angleY);
        double dirZ = std::cos(angleX);

        currX += dirX * CA_DIST;
        currY += dirY * CA_DIST;
        currZ += dirZ * CA_DIST;

        res.x = currX; res.y = currY; res.z = currZ;
        protein.push_back(res);
    }

    std::ofstream out("Universe_v8_3_Extrusion.pdb");
    out << "HEADER    PROJECT-137 V8.3 EXTRUSION MODEL\n";
    for (auto& r : protein) {
        char buf[100];
        sprintf(buf, "ATOM  %5d  CA  %3s A%4d    %8.3f%8.3f%8.3f  1.00  4.00           C\n",
                r.resSeq, r.resName.c_str(), r.resSeq, r.x, r.y, r.z);
        out << buf;
    }
    std::cout << "[Success] v8.3 Generated. This is the 'Force of Life'." << std::endl;
    return 0;
}


