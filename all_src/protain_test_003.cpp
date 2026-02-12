#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <cmath>
#include <iomanip>
#include <map>

// --- 宇宙OS 定数 ---
const double PI = 3.1415926535;
const double GOLDEN_ANGLE = PI * (3.0 - std::sqrt(5.0));
const double RAD_TO_DEG = 180.0 / PI;

// --- 原子構造体 ---
struct Atom {
    int serial;
    std::string name; // "CA", "CB", "N", "O" etc.
    std::string element; // "C", "N", "O", "S"
    double x, y, z;
};

// --- アミノ酸（残基）構造体 ---
struct Residue {
    int resSeq;
    char code; // 'A', 'G', etc.
    std::string resName; // "ALA", "GLY"
    bool isHydrophobic;
    std::vector<Atom> atoms; // この残基に含まれる全原子
    double centerX, centerY, centerZ; // CAの位置（基準点）
};

// --- 3文字コード変換 & 側鎖サイズ定義 ---
// 戻り値: {3文字名, 側鎖の原子数(簡易)}
std::pair<std::string, int> getAAInfo(char code) {
    static std::map<char, std::pair<std::string, int>> m = {
        {'A', {"ALA", 1}}, {'R', {"ARG", 6}}, {'N', {"ASN", 3}}, {'D', {"ASP", 3}},
        {'C', {"CYS", 1}}, {'Q', {"GLN", 4}}, {'E', {"GLU", 4}}, {'G', {"GLY", 0}},
        {'H', {"HIS", 5}}, {'I', {"ILE", 3}}, {'L', {"LEU", 3}}, {'K', {"LYS", 4}},
        {'M', {"MET", 4}}, {'F', {"PHE", 6}}, {'P', {"PRO", 2}}, {'S', {"SER", 1}},
        {'T', {"THR", 2}}, {'W', {"TRP", 9}}, {'Y', {"TYR", 7}}, {'V', {"VAL", 2}}
    };
    if (m.count(code)) return m[code];
    return {"UNK", 0};
}

// 疎水性判定
bool isHydrophobic(char aa) {
    return std::string("LIVFMWAC").find(aa) != std::string::npos;
}

// コドン変換
char codonToAA(const std::string& codon) {
    // (短縮のため主要なもののみ記述、実際は全テーブルが必要)
    // ここではテスト用配列が動くように完全なマップを想定
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
    if (m.count(codon)) return m[codon];
    return '?';
}

// --- ベクトル計算用 ---
void normalize(double& x, double& y, double& z) {
    double len = std::sqrt(x*x + y*y + z*z);
    if (len > 0) { x /= len; y /= len; z /= len; }
}

int main() {
    std::string rnaSeq;
    std::cout << "Enter RNA: ";
    std::cin >> rnaSeq;

    std::vector<Residue> protein;
    int atomSerialCounter = 1;
    
    // 1. 基本骨格（Sphere Logic）の計算
    // まずはCA（アルファ炭素）の位置だけを決める
    std::vector<char> seqCodes;
    for (size_t i = 0; i < rnaSeq.length(); i += 3) {
        if (i+3 > rnaSeq.length()) break;
        char c = codonToAA(rnaSeq.substr(i, 3));
        if (c == '_' || c == '?') break;
        seqCodes.push_back(c);
    }

    double radiusBase = 6.0 + (seqCodes.size() * 0.12);
    double dy = 2.0 / (double)seqCodes.size();

    for (size_t i = 0; i < seqCodes.size(); ++i) {
        Residue res;
        res.resSeq = i + 1;
        res.code = seqCodes[i];
        res.isHydrophobic = isHydrophobic(res.code);
        auto info = getAAInfo(res.code);
        res.resName = info.first;

        // 球面配置
        double myY = 1.0 - (i * dy);
        double r_plane = std::sqrt(1.0 - myY * myY);
        double theta = i * GOLDEN_ANGLE;

        // 疎水性レイヤー調整
        double layer = radiusBase;
        if (res.isHydrophobic) layer *= 0.6; 
        else layer *= 1.2;

        // CA座標の確定
        res.centerX = std::cos(theta) * r_plane * layer;
        res.centerZ = std::sin(theta) * r_plane * layer;
        res.centerY = myY * radiusBase;

        // --- 詳細原子（Atom）の生成 (Procedural Generation) ---
        
        // 1. Backbone (N, CA, C, O)
        // 宇宙OSの手抜き: 厳密な結合角ではなく、CAの周囲に簡易配置する
        // CA
        res.atoms.push_back({atomSerialCounter++, "CA", "C", res.centerX, res.centerY, res.centerZ});
        
        // N (少し上に)
        res.atoms.push_back({atomSerialCounter++, "N", "N", res.centerX, res.centerY + 1.4, res.centerZ});
        
        // C (少し下に)
        res.atoms.push_back({atomSerialCounter++, "C", "C", res.centerX, res.centerY - 1.5, res.centerZ});
        
        // O (Cの横に)
        res.atoms.push_back({atomSerialCounter++, "O", "O", res.centerX + 1.2, res.centerY - 1.5, res.centerZ});

        // 2. Side Chain (側鎖)
        // 中心から外向き（あるいは内向き）のベクトルを計算
        double vecX = res.centerX;
        double vecY = res.centerY;
        double vecZ = res.centerZ;
        normalize(vecX, vecY, vecZ);

        // 疎水性なら内側へ、親水性なら外側へ原子を伸ばす
        double direction = res.isHydrophobic ? -1.0 : 1.0;
        int sideChainCount = info.second;

        for (int k = 0; k < sideChainCount; ++k) {
            std::string atomName = "C" + std::string(1, 'B' + k); // CB, CC, CD... (簡易命名)
            std::string elem = "C";
            // 末端原子の種類の簡易判定 (NH2とかOHとか)
            if (k == sideChainCount - 1) {
                if (res.code == 'K' || res.code == 'R') elem = "N"; // Lys, Argは窒素
                if (res.code == 'D' || res.code == 'E' || res.code == 'S') elem = "O"; // 酸素
            }

            double dist = 1.5 * (k + 1); // 結合距離
            double atomX = res.centerX + (vecX * dist * direction);
            double atomY = res.centerY + (vecY * dist * direction); // Y方向にも少し広げる
            double atomZ = res.centerZ + (vecZ * dist * direction);

            // 乱数的なゆらぎ（137由来）を加えて、直線にならないようにする
            double wobble = std::sin((i + k) * 137.0) * 0.5;
            atomX += wobble;

            res.atoms.push_back({atomSerialCounter++, atomName, elem, atomX, atomY, atomZ});
        }

        protein.push_back(res);
    }

    // --- PDBファイル出力 ---
    std::ofstream file("universe_os_full_render.pdb");
    for (const auto& res : protein) {
        for (const auto& atom : res.atoms) {
            // PDBフォーマット整形 (桁合わせ重要)
            file << "ATOM  " 
                 << std::setw(5) << atom.serial << " "
                 << std::setw(4) << std::left << atom.name << " " // 左寄せ
                 << std::setw(3) << res.resName << " A"
                 << std::setw(4) << std::right << res.resSeq << "    " // 右寄せ
                 << std::setw(8) << std::fixed << std::setprecision(3) << atom.x
                 << std::setw(8) << std::fixed << std::setprecision(3) << atom.y
                 << std::setw(8) << std::fixed << std::setprecision(3) << atom.z
                 << "  1.00  4.00           "
                 << atom.element << "  " << std::endl;
        }
    }
    file << "END" << std::endl;
    std::cout << "Full atomic model generated: universe_os_full_render.pdb" << std::endl;

    return 0;
}