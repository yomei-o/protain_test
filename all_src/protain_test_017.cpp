#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <cmath>
#include <iomanip>
#include <map>

// --- Project-137: Universe OS v8.0 [Restored] ---
// 最も「立体的で美しかった」黄金螺旋ロジックへの原点回帰

const double PI = 3.1415926535;
const double GOLDEN_ANGLE = PI * (3.0 - std::sqrt(5.0)); // 137.508...
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
    std::cout << "Project-137 [v8.0 Restore: Cosmic Backbone]" << std::endl;
    std::cout << "Input: ";
    std::cin >> inputSeq;

    std::vector<Residue> protein;
    
    // オキシトシン環（C...C）の検知（CONECT行用）
    int c1 = -1, c6 = -1;
    for(int i=0; i < (int)inputSeq.length() - 5; i++) {
        if(inputSeq[i] == 'C' && inputSeq[i+5] == 'C') { c1 = i; c6 = i + 5; break; }
    }

    for (int i = 0; i < (int)inputSeq.length(); ++i) {
        Residue res;
        res.resSeq = i + 1;
        res.code = inputSeq[i];
        res.resName = getAAInfo(res.code).first;

        // --- 【v8.0: Cosmic Spiral Engine】 ---
        // 黄金角による回転と、一定のZ軸上昇が生む「宇宙の螺旋」
        double phi = i * GOLDEN_ANGLE;
        double radius = 5.0 + std::sin(i * 0.5) * 2.0; // 緩やかな半径の変化
        
        res.x = radius * std::cos(phi);
        res.y = radius * std::sin(phi);
        res.z = i * 2.5; // 適度なピッチで上昇
        
        protein.push_back(res);
    }

    std::ofstream out("Universe_v8_0_Restored.pdb");
    out << "HEADER    PROJECT-137 V8.0 RESTORED\n";
    for (auto& r : protein) {
        char buf[100];
        sprintf(buf, "ATOM  %5d  CA  %3s A%4d    %8.3f%8.3f%8.3f  1.00  4.00           C\n",
                r.resSeq, r.resName.c_str(), r.resSeq, r.x, r.y, r.z);
        out << buf;
    }
    // オキシトシンの「指輪」を繋ぐ
    if(c1 != -1) out << "CONECT " << std::setw(5) << c1 + 1 << std::setw(5) << c6 + 1 << "\n";
    out << "END\n";
    
    std::cout << "[Complete] Reverted to the most beautiful cosmic logic." << std::endl;
    return 0;
}
