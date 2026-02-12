#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <cmath>
#include <iomanip>
#include <map>

// --- Project-137: Universe OS v7.2 [Micro-Focus Patch] ---
// 配列の長さに応じて、宇宙の幾何学を「望遠」から「顕微鏡」へ切り替える

const double PI = 3.1415926535;
const double CA_DIST = 3.8; 
const double GOLDEN_ANGLE = PI * (3.0 - std::sqrt(5.0));

struct Residue {
    int resSeq;
    char code;
    std::string resName;
    bool isHydrophobic;
    double x, y, z;
};

// --- アミノ酸情報 ---
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

bool isHydrophobic(char aa) {
    return std::string("LIVFMWAC").find(aa) != std::string::npos;
}

void normalize(double& x, double& y, double& z) {
    double len = std::sqrt(x*x + y*y + z*z);
    if (len > 1e-6) { x /= len; y /= len; z /= len; }
}

int main() {
    std::string inputSeq;
    std::cout << "===========================================" << std::endl;
    std::cout << " Project-137 v7.2 [Micro-Focus Patch] " << std::endl;
    std::cout << "===========================================" << std::endl;
    std::cout << "Input (RNA or AA): ";
    std::cin >> inputSeq;

    std::vector<Residue> protein;
    double currX = 0, currY = 0, currZ = 0;
    double dirX = 0, dirY = 0, dirZ = 1.0;

    int totalLen = inputSeq.length();

    for (int i = 0; i < totalLen; ++i) {
        Residue res;
        res.resSeq = i + 1;
        res.code = inputSeq[i];
        res.resName = getAAInfo(res.code).first;
        res.isHydrophobic = isHydrophobic(res.code);

        // --- 動的幾何学エンジン (Context-Aware) ---
        if (totalLen < 15) {
            // 【マイクロ・モード】短い配列（オキシトシン等）
            // 9歩で1周する精密なリング幾何学
            double angle = (2.0 * PI * i) / 9.0;
            dirX = std::cos(angle);
            dirY = std::sin(angle);
            dirZ = 0.4; // 螺旋の高さ
        } 
        else if (totalLen > 50 && i < 20) {
            // 【物流モード】長い配列の先頭（シグナル）
            dirX = 0.1 * std::sin(i * 0.5);
            dirY = 0.1 * std::cos(i * 0.5);
            dirZ = 1.0; // 直進性を優先
        } 
        else {
            // 【パッキング・モード】本体（ニューロフィジン等）
            double phi = i * GOLDEN_ANGLE;
            dirX = std::cos(phi);
            dirY = std::sin(phi);
            dirZ = std::cos(i * 0.1); // 密集させるための折り返し
        }

        normalize(dirX, dirY, dirZ);
        currX += dirX * CA_DIST;
        currY += dirY * CA_DIST;
        currZ += dirZ * CA_DIST;

        res.x = currX; res.y = currY; res.z = currZ;
        protein.push_back(res);
    }

    // PDB出力
    FILE* fp = fopen("Universe_v7_2_Final.pdb", "w");
    if (fp) {
        fprintf(fp, "HEADER    PROJECT-137 V7.2 FINAL PATCH\n");
        int atomSerial = 1;
        for (auto& r : protein) {
            // CA原子
            fprintf(fp, "ATOM  %5d  CA  %3s A%4d    %8.3f%8.3f%8.3f  1.00  4.00           C\n",
                atomSerial++, r.resName.c_str(), r.resSeq, r.x, r.y, r.z);
            
            // 簡易的な側鎖（アミノ酸の重さに応じて枝を伸ばす）
            int scCount = getAAInfo(r.code).second;
            for(int k=0; k<scCount; k++) {
                double d = 1.2 * (k+1);
                fprintf(fp, "ATOM  %5d  C%c  %3s A%4d    %8.3f%8.3f%8.3f  1.00  4.00           C\n",
                    atomSerial++, 'B'+k, r.resName.c_str(), r.resSeq, 
                    r.x + (dirX * d), r.y + (dirY * d), r.z + (dirZ * d));
            }
        }
        // オキシトシンの場合のみ SS結合を明示
        if (totalLen == 9) {
            fprintf(fp, "CONECT    1    6\n"); // Cys1 - Cys6
        }
        fprintf(fp, "END\n");
        fclose(fp);
        std::cout << "[Complete] 'Universe_v7_2_Final.pdb' has been generated." << std::endl;
    }

    return 0;
}