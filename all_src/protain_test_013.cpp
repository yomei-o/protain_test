#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <cmath>
#include <iomanip>
#include <map>

// --- Project-137: Universe OS v9.0 "The Weaver" ---
// コンセプト: 螺旋をほどき、ジグザグの反復によって「面（シート）」を構築する

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
        {'A', {"ALA", 1}}, {'G', {"GLY", 0}}, {'V', {"VAL", 2}}, {'S', {"SER", 1}},
        {'C', {"CYS", 1}}, {'I', {"ILE", 3}}, {'L', {"LEU", 3}}, {'Y', {"TYR", 7}},
        {'Q', {"GLN", 4}}, {'N', {"ASN", 3}}, {'P', {"PRO", 2}}, {'M', {"MET", 4}}
    };
    return m.count(code) ? m[code] : std::make_pair("UNK", 0);
}

int main() {
    std::string inputSeq;
    std::cout << "Project-137 v9.0 [The Weaver Mode]" << std::endl;
    std::cout << "Sequence (Try: GAGAGAGAGAGAGAGAGAGA): ";
    std::cin >> inputSeq;

    std::vector<Residue> protein;
    double currX = 0, currY = 0, currZ = 0;
    
    // シート形成のための周期フラグ
    bool isSheetPattern = (inputSeq.find("GAGA") != std::string::npos || inputSeq.find("VTV") != std::string::npos);

    for (int i = 0; i < (int)inputSeq.length(); ++i) {
        Residue res;
        res.resSeq = i + 1;
        res.code = inputSeq[i];
        res.resName = getAAInfo(res.code).first;

        double dirX, dirY, dirZ;

        if (isSheetPattern) {
            // --- 【Weaver Logic: ジグザグ・エンジン】 ---
            // 1歩ごとに左右に 120度 振り、前進し続けることで「伸びた鎖」を作る
            double sideSweep = (i % 2 == 0) ? 1.0 : -1.0;
            dirX = 0.5 * sideSweep; // 左右のジグザグ
            dirY = 0.2;            // わずかな厚み
            dirZ = 1.0;            // メインの前進軸
        } else {
            // 既存の Cosmic Backbone (螺旋モード)
            double phi = i * PI * (3.0 - std::sqrt(5.0));
            dirX = std::cos(phi);
            dirY = std::sin(phi);
            dirZ = 0.5;
        }

        // ベクトルの正規化
        double len = std::sqrt(dirX*dirX + dirY*dirY + dirZ*dirZ);
        dirX /= len; dirY /= len; dirZ /= len;

        currX += dirX * CA_DIST;
        currY += dirY * CA_DIST;
        currZ += dirZ * CA_DIST;

        res.x = currX; res.y = currY; res.z = currZ;
        protein.push_back(res);
    }

    std::ofstream out("Universe_v9_0_Weaver.pdb");
    out << "HEADER    PROJECT-137 V9.0 THE WEAVER\n";
    for (auto& r : protein) {
        char buf[100];
        sprintf(buf, "ATOM  %5d  CA  %3s A%4d    %8.3f%8.3f%8.3f  1.00  4.00           C\n",
                r.resSeq, r.resName.c_str(), r.resSeq, r.x, r.y, r.z);
        out << buf;
    }
    out << "END\n";
    
    std::cout << "[Success] Generated Universe_v9_0_Weaver.pdb" << std::endl;
    if(isSheetPattern) std::cout << "[System] Sheet Pattern Detected. Zigzag engine engaged." << std::endl;

    return 0;
}