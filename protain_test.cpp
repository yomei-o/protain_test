#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <cmath>
#include <iomanip>
#include <map>

// --- 宇宙OS 定数定義 ---
const double PI = 3.1415926535;
// 宇宙の「巻き取り」黄金角 (Golden Angle)
// 植物の葉の配置や銀河の螺旋と同じ、最も効率よく球を埋める角度
const double GOLDEN_ANGLE = PI * (3.0 - std::sqrt(5.0)); 

struct AminoAcid {
    char code;
    double x, y, z;
    bool isHydrophobic;
};

// 疎水性判定
bool isHydrophobic(char aa) {
    // 内側に隠れたがるアミノ酸たち
    std::string hydrophobic = "LIVFMWAC"; 
    return hydrophobic.find(aa) != std::string::npos;
}

// コドン変換 (省略せずに書きます)
char codonToAA(const std::string& codon) {
    static std::map<std::string, char> m = {
        {"UUU",'F'}, {"UUC",'F'}, {"UUA",'L'}, {"UUG",'L'},
        {"CUU",'L'}, {"CUC",'L'}, {"CUA",'L'}, {"CUG",'L'},
        {"AUU",'I'}, {"AUC",'I'}, {"AUA",'I'}, {"AUG",'M'},
        {"GUU",'V'}, {"GUC",'V'}, {"GUA",'V'}, {"GUG",'V'},
        {"UCU",'S'}, {"UCC",'S'}, {"UCA",'S'}, {"UCG",'S'},
        {"CCU",'P'}, {"CCC",'P'}, {"CCA",'P'}, {"CCG",'P'},
        {"ACU",'T'}, {"ACC",'T'}, {"ACA",'T'}, {"ACG",'T'},
        {"GCU",'A'}, {"GCC",'A'}, {"GCA",'A'}, {"GCG",'A'},
        {"UAU",'Y'}, {"UAC",'Y'}, {"UAA",'_'}, {"UAG",'_'},
        {"CAU",'H'}, {"CAC",'H'}, {"CAA",'Q'}, {"CAG",'Q'},
        {"AAU",'N'}, {"AAC",'N'}, {"AAA",'K'}, {"AAG",'K'},
        {"GAU",'D'}, {"GAC",'D'}, {"GAA",'E'}, {"GAG",'E'},
        {"UGU",'C'}, {"UGC",'C'}, {"UGA",'_'}, {"UGG",'W'},
        {"CGU",'R'}, {"CGC",'R'}, {"CGA",'R'}, {"CGG",'R'},
        {"AGU",'S'}, {"AGC",'S'}, {"AGA",'R'}, {"AGG",'R'},
        {"GGU",'G'}, {"GGC",'G'}, {"GGA",'G'}, {"GGG",'G'}
    };
    if (m.find(codon) != m.end()) return m[codon];
    return '?';
}

void savePDB(const std::vector<AminoAcid>& protein, const std::string& filename) {
    std::ofstream file(filename);
    int atomSerial = 1;
    for (size_t i = 0; i < protein.size(); ++i) {
        file << "ATOM  " << std::setw(5) << atomSerial << "  CA  GLY A" 
             << std::setw(4) << (i + 1) << "    "
             << std::setw(8) << std::fixed << std::setprecision(3) << protein[i].x
             << std::setw(8) << std::fixed << std::setprecision(3) << protein[i].y
             << std::setw(8) << std::fixed << std::setprecision(3) << protein[i].z
             << "  1.00  4.00           C  " << std::endl;
        atomSerial++;
        
        // 視認性を良くするために、原子間を結ぶCONECTレコードを追加（疑似的な結合）
        if (i > 0) {
            file << "CONECT" 
                 << std::setw(5) << (atomSerial - 1) 
                 << std::setw(5) << (atomSerial - 2) << std::endl;
        }
    }
    file << "END" << std::endl;
    std::cout << "Generated: " << filename << std::endl;
}

int main() {
    std::string rnaSeq;
    std::cout << "Enter RNA: ";
    std::cin >> rnaSeq;

    std::vector<AminoAcid> protein;
    
    // 1. 翻訳
    for (size_t i = 0; i < rnaSeq.length(); i += 3) {
        if (i + 3 > rnaSeq.length()) break;
        char code = codonToAA(rnaSeq.substr(i, 3));
        if (code == '_') break;
        AminoAcid aa;
        aa.code = code;
        aa.isHydrophobic = isHydrophobic(code);
        protein.push_back(aa);
    }

    // 2. 宇宙OS「球体パッキング」ロジック (Spherical Folding)
    // 紐をダラダラ伸ばすのではなく、最初から「球」として巻いていく
    
    double radiusBase = 5.0 + (protein.size() * 0.15); // タンパク質の大きさに応じた基本半径
    double y = 0;
    double dy = 2.0 / (double)protein.size(); // 高さの刻み幅

    for (size_t i = 0; i < protein.size(); ++i) {
        // 球面上の螺旋配置 (Fibonacci Sphere Algorithmの応用)
        double myY = 1.0 - (i * dy); // 1 to -1
        double r_plane = std::sqrt(1.0 - myY * myY); // 平面上の半径
        
        double theta = i * GOLDEN_ANGLE; // 黄金角で回転

        // --- ここが宇宙OSのハック ---
        // 疎水性なら「内側（核）」へ、親水性なら「外側（殻）」へ半径をズラす
        double layer = radiusBase;
        if (protein[i].isHydrophobic) {
            layer *= 0.6; // 内側に潜る (Core)
        } else {
            layer *= 1.2; // 外側に出る (Surface)
        }

        // 座標変換 (極座標 -> 直交座標)
        protein[i].x = std::cos(theta) * r_plane * layer;
        protein[i].z = std::sin(theta) * r_plane * layer;
        protein[i].y = myY * radiusBase; // 高さはそのまま
    }

    savePDB(protein, "universe_os_sphere.pdb");
    return 0;
}


