#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <cmath>
#include <iomanip>
#include <map>

// --- Project-137: Universe OS v7.1 (High-Speed & Universal) ---
// コンセプト: 連続性を保ち、あらゆるタンパク質に宇宙の幾何学を適用する

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

// --- アミノ酸情報ライブラリ ---
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

// --- 宇宙OS 数学エンジン ---
void normalize(double& x, double& y, double& z) {
    double len = std::sqrt(x*x + y*y + z*z);
    if (len > 1e-6) { x /= len; y /= len; z /= len; }
}

int main() {
    std::string inputSeq;
    std::cout << "Project-137 v7.1 [Universal Mode]" << std::endl;
    std::cout << "Sequence: ";
    std::cin >> inputSeq;

    std::vector<Residue> protein;
    double currX = 0, currY = 0, currZ = 0;
    double dirX = 0, dirY = 0, dirZ = 1.0;

    // 1. 疎水性密度のスキャン（汎用性のための自己学習）
    double h_density = 0;
    for(char c : inputSeq) if(isHydrophobic(c)) h_density++;
    h_density /= inputSeq.length();

    for (size_t i = 0; i < inputSeq.length(); ++i) {
        Residue res;
        res.resSeq = i + 1;
        res.code = inputSeq[i];
        res.resName = getAAInfo(res.code).first;
        res.isHydrophobic = isHydrophobic(res.code);

        // 2. 宇宙の進行ベクトルの決定
        if (i < 20 && inputSeq.length() > 50) {
            // シグナル領域: 直進
            dirX += 0.1 * std::sin(i);
            dirY += 0.1 * std::cos(i);
            dirZ = 1.0;
        } else {
            // 本体領域: 黄金角と疎水性密度によるパッキング
            double phi = i * GOLDEN_ANGLE;
            // 疎水性が高いほど、中心に向かう力が強くなる
            double twist = h_density * 2.0 * PI;
            dirX = std::cos(phi + twist);
            dirY = std::sin(phi + twist);
            dirZ = std::cos(i * 0.1); 
        }

        normalize(dirX, dirY, dirZ);
        currX += dirX * CA_DIST;
        currY += dirY * CA_DIST;
        currZ += dirZ * CA_DIST;

        res.x = currX; res.y = currY; res.z = currZ;
        protein.push_back(res);
    }

    // 3. PDB出力
    std::ofstream out("Universe_v7_1.pdb");
    out << "HEADER    PROJECT-137 UNIVERSAL BRIDGE\n";
    for (auto& r : protein) {
        char buf[100];
        sprintf(buf, "ATOM  %5d  CA  %3s A%4d    %8.3f%8.3f%8.3f  1.00  4.00           C\n",
                r.resSeq, r.resName.c_str(), r.resSeq, r.x, r.y, r.z);
        out << buf;
        
        // 側鎖（宇宙の枝）
        int sc = getAAInfo(r.code).second;
        for(int k=0; k<sc; k++) {
            double d = 1.5 * (k+1);
            sprintf(buf, "ATOM  %5d  C%c  %3s A%4d    %8.3f%8.3f%8.3f  1.00  4.00           C\n",
                    r.resSeq*10+k, 'B'+k, r.resName.c_str(), r.resSeq, 
                    r.x + (dirX * d), r.y + (dirY * d), r.z + (dirZ * d));
            out << buf;
        }
    }
    out << "END\n";
    std::cout << "[Success] Generated Universe_v7_1.pdb" << std::endl;

    return 0;
}

