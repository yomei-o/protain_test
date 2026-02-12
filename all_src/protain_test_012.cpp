#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <cmath>
#include <iomanip>
#include <map>

// --- Project-137: Universe OS v8.0 "Cosmic Backbone" ---
// コンセプト: 幾何学的な「円」を卒業し、物理的な「しなり」を実装する

const double PI = 3.1415926535;
const double CA_DIST = 3.8; 

struct Residue {
    int resSeq;
    char code;
    std::string resName;
    double x, y, z;
};

// アミノ酸ライブラリ
std::pair<std::string, int> getAAInfo(char code) {
    static std::map<char, std::pair<std::string, int>> m = {
        {'A', {"ALA", 1}}, {'R', {"ARG", 6}}, {'N', {"ASN", 3}}, {'D', {"ASP", 3}},
        {'C', {"CYS", 1}}, {'Q', {"GLN", 4}}, {'E', {"GLU", 4}}, {'G', {"GLY", 0}},
        {'H', {"HIS", 5}}, {'I', {"ILE", 3}}, {'L', {"LEU", 3}}, {'K', {"LYS", 4}},
        {'M', {"MET", 4}}, {'F', {"PHE", 6}}, {'P', {"PRO", 2}}, {'S', {"SER", 1}},
        {'T', {"THR", 2}}, {'W', {"TRP", 9}}, {'Y', {"TYR", 7}}, {'V', {"VAL", 2}}
    };
    return m.count(code) ? m[code] : std::make_pair("UNK", 0);
}

void normalize(double& x, double& y, double& z) {
    double len = std::sqrt(x*x + y*y + z*z);
    if (len > 1e-6) { x /= len; y /= len; z /= len; }
}

int main() {
    std::string inputSeq;
    std::cout << "Project-137 v8.0 [Cosmic Backbone]" << std::endl;
    std::cout << "Input: ";
    std::cin >> inputSeq;

    std::vector<Residue> protein;
    double currX = 0, currY = 0, currZ = 0;
    
    // 慣性ベクトル（滑らかな曲線の維持）
    double dirX = 0, dirY = 0, dirZ = 1.0;
    double angleX = 0, angleY = 0;

    // オキシトシン環の自動検知 (Cys...Cys)
    int c1 = -1, c6 = -1;
    for(int i=0; i < (int)inputSeq.length() - 5; i++) {
        if(inputSeq[i] == 'C' && inputSeq[i+5] == 'C') {
            c1 = i; c6 = i + 5; break; // 最初に見つかったペアを採用
        }
    }

    for (int i = 0; i < (int)inputSeq.length(); ++i) {
        Residue res;
        res.resSeq = i + 1;
        res.code = inputSeq[i];
        res.resName = getAAInfo(res.code).first;

        // --- 宇宙OS 流体動力学エンジン ---
        if (c1 != -1 && i >= c1 && i <= c6) {
            // 【環形成モード】 
            // 6つのアミノ酸で 180度以上曲げて戻ってくる「Uターン」を描く
            double progress = (double)(i - c1) / (c6 - c1);
            angleX += 0.8; // 急カーブ
            angleY += 0.4 * std::sin(progress * PI);
        } 
        else if (i < 20 && inputSeq.length() > 50) {
            // 【シグナルモード】 直進しながら微細な螺旋 (Alpha-Helix)
            angleX += 0.1;
            angleY += 0.5;
        } 
        else {
            // 【フォールディングモード】 黄金角によるランダム・ウォーク
            angleX += (PI * (3.0 - std::sqrt(5.0)));
            angleY += 0.3;
        }

        // 角度からベクトルへ変換（慣性を考慮）
        dirX = std::sin(angleX) * std::cos(angleY);
        dirY = std::sin(angleX) * std::sin(angleY);
        dirZ = std::cos(angleX);

        normalize(dirX, dirY, dirZ);
        currX += dirX * CA_DIST;
        currY += dirY * CA_DIST;
        currZ += dirZ * CA_DIST;

        res.x = currX; res.y = currY; res.z = currZ;
        protein.push_back(res);
    }

    // PDB出力
    std::ofstream out("Universe_v8_0_Final.pdb");
    out << "HEADER    PROJECT-137 V8.0 COSMIC BACKBONE\n";
    for (auto& r : protein) {
        char buf[100];
        sprintf(buf, "ATOM  %5d  CA  %3s A%4d    %8.3f%8.3f%8.3f  1.00  4.00           C\n",
                r.resSeq, r.resName.c_str(), r.resSeq, r.x, r.y, r.z);
        out << buf;
    }
    // 正しい接続情報の書き出し
    if (c1 != -1) {
        out << "CONECT " << std::setw(5) << c1 + 1 << std::setw(5) << c6 + 1 << "\n";
    }
    out << "END\n";
    
    std::cout << "[Complete] 'Universe_v8_0_Final.pdb' generated." << std::endl;
    if(c1 != -1) std::cout << "[System] Disulfide Bridge detected: " << c1+1 << " - " << c6+1 << std::endl;

    return 0;
}

