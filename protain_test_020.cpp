#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <cmath>
#include <iomanip>
#include <map>

// --- Project-137: Universe OS v9.0 "Helix Gravity" ---
// コンセプト: 直線を廃止し、すべての線に「重力によるねじれ」を与える

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
    std::cout << "Project-137 v9.0 [Helix Gravity]" << std::endl;
    std::cout << "Input: ";
    std::cin >> inputSeq;

    std::vector<Residue> protein;
    double currX = 0, currY = 0, currZ = 0;
    
    // 慣性ベクトル
    double dirX = 0, dirY = 0, dirZ = 1.0;
    
    // 角度パラメータ
    double angleX = 0; // 緯度（Z軸との角度）
    double angleY = 0; // 経度（XY平面の回転）

    // オキシトシン環の自動検知 (Cys...Cys)
    int c1 = -1, c6 = -1;
    for(int i=0; i < (int)inputSeq.length() - 5; i++) {
        if(inputSeq[i] == 'C' && inputSeq[i+5] == 'C') {
            c1 = i; c6 = i + 5; break; 
        }
    }

    for (int i = 0; i < (int)inputSeq.length(); ++i) {
        Residue res;
        res.resSeq = i + 1;
        res.code = inputSeq[i];
        res.resName = getAAInfo(res.code).first;

        // --- 宇宙OS 流体動力学エンジン v9 ---
        
        if (c1 != -1 && i >= c1 && i <= c6) {
            // 【環形成モード】 Uターン
            double progress = (double)(i - c1) / (c6 - c1);
            angleX += 0.8; 
            angleY += 0.4 * std::sin(progress * PI);
        } 
        else if (i < 20 && inputSeq.length() > 50) {
            // 【シグナルモード：螺旋重力パッチ】
            // 直線にするな！ねじれろ！
            
            // 1. 緯度を固定する（これ重要！）
            // angleXを増やし続けると「円」になるが、固定すると「直線的な螺旋」になる
            angleX = 0.6; // Z軸から約35度開く（傘のような開き具合）
            
            // 2. 経度を回す
            // 1ステップごとに約60度回す（トイレットペーパーのねじれ）
            angleY += 1.0; 
            
            // ※おまけ：少しだけ「ゆらぎ」を入れて、人工的な完全螺旋を防ぐ
            angleX += 0.05 * std::sin(i * 0.5);
        } 
        else {
            // 【フォールディングモード】 黄金角ランダム
            angleX += (PI * (3.0 - std::sqrt(5.0)));
            angleY += 0.3;
        }

        // 角度からベクトルへ変換
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
    std::ofstream out("Universe_v9_0_Gravity.pdb");
    out << "HEADER    PROJECT-137 V9.0 HELIX GRAVITY\n";
    for (auto& r : protein) {
        char buf[100];
        sprintf(buf, "ATOM  %5d  CA  %3s A%4d    %8.3f%8.3f%8.3f  1.00  4.00           C\n",
                r.resSeq, r.resName.c_str(), r.resSeq, r.x, r.y, r.z);
        out << buf;
    }
    if (c1 != -1) {
        out << "CONECT " << std::setw(5) << c1 + 1 << std::setw(5) << c6 + 1 << "\n";
    }
    out << "END\n";
    
    std::cout << "[Complete] 'Universe_v9_0_Gravity.pdb' generated." << std::endl;
    std::cout << "[Info] Signal peptide (1-20) is now twisted." << std::endl;

    return 0;
}

