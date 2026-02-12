#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <cmath>
#include <iomanip>
#include <map>
#include <algorithm>

// --- 宇宙OS 定数定義 ---
const double PI = 3.1415926535;
const double CA_DIST = 3.8;   // CA間距離 (Angstrom)
const double GOLDEN_ANGLE = PI * (3.0 - std::sqrt(5.0));

// --- データ構造 ---
struct Atom {
    int serial;
    std::string name;
    std::string element;
    double x, y, z;
};

struct Residue {
    int resSeq;
    char code;
    std::string resName;
    bool isHydrophobic;
    int structureType; // 0:Loop, 1:Helix, 2:Sheet, 3:Collagen
    double x, y, z; // CA座標
    std::vector<Atom> atoms;
};

// --- ヘルパー関数群 ---

// 3文字コード変換 & 側鎖原子数
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

bool isHydrophobic(char aa) {
    return std::string("LIVFMWAC").find(aa) != std::string::npos;
}

// 構造予測ロジック (v4.0: コラーゲン検知機能付き)
int predictStructure(const std::vector<char>& seq, size_t index) {
    char c = seq[index];
    char next1 = (index + 1 < seq.size()) ? seq[index+1] : ' ';
    char next2 = (index + 2 < seq.size()) ? seq[index+2] : ' ';

    // 1. コラーゲン検知 (G-P-P などの繰り返しパターン)
    // グリシン(G)やプロリン(P)が密集しているエリアは「ねじれ」を起こす
    int collagenScore = 0;
    if (c == 'G' || c == 'P') collagenScore++;
    if (next1 == 'G' || next1 == 'P') collagenScore++;
    if (next2 == 'G' || next2 == 'P') collagenScore++;
    
    if (collagenScore >= 2) return 3; // Collagen Type

    // 2. Alpha Helix (螺旋)
    std::string helixFormers = "EMALKQH";
    if (helixFormers.find(c) != std::string::npos) return 1;

    // 3. Beta Sheet (壁)
    std::string sheetFormers = "VIYFWCT";
    if (sheetFormers.find(c) != std::string::npos) return 2;

    return 0; // Loop (その他)
}

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
    if (m.count(codon)) return m[codon];
    return '?';
}

void normalize(double& x, double& y, double& z) {
    double len = std::sqrt(x*x + y*y + z*z);
    if (len > 1e-6) { x /= len; y /= len; z /= len; }
}

int main() {
    std::string rnaSeq;
    std::cout << "--- Universe OS Protein Architect v4.0 ---" << std::endl;
    std::cout << "Target RNA Sequence: ";
    std::cin >> rnaSeq;

    std::vector<Residue> protein;
    std::vector<char> fullSeq;

    // 1. まず全配列を翻訳
    for (size_t i = 0; i < rnaSeq.length(); i += 3) {
        if (i+3 > rnaSeq.length()) break;
        char c = codonToAA(rnaSeq.substr(i, 3));
        if (c == '_' || c == '?') break;
        fullSeq.push_back(c);
    }

    // 2. 構造予測と骨格生成 (Turtle Graphics)
    double currX = 0, currY = 0, currZ = 0;
    double dirX = 1.0, dirY = 0.0, dirZ = 0.0; // 進行方向

    for (size_t i = 0; i < fullSeq.size(); ++i) {
        Residue res;
        res.resSeq = i + 1;
        res.code = fullSeq[i];
        res.isHydrophobic = isHydrophobic(res.code);
        auto info = getAAInfo(res.code);
        res.resName = info.first;
        
        // 構造タイプの決定
        res.structureType = predictStructure(fullSeq, i);

        // --- 移動ロジック (Movement Logic) ---
        
        if (res.structureType == 1) { 
            // Type 1: Alpha Helix (標準的な螺旋)
            // XY平面で回転しつつ、Z方向へ少し進む
            double angle = 100.0 * (PI / 180.0);
            double newDirX = dirX * std::cos(angle) - dirY * std::sin(angle);
            double newDirY = dirX * std::sin(angle) + dirY * std::cos(angle);
            dirX = newDirX; dirY = newDirY;
            dirZ = 0.5; // バネのピッチ
        } 
        else if (res.structureType == 2) {
            // Type 2: Beta Sheet (壁)
            // ジグザグ運動
            dirY = (i % 2 == 0) ? 0.8 : -0.8;
            dirZ = 0.2; // 少しずつ進む
        }
        else if (res.structureType == 3) {
            // Type 3: Collagen Helix (コラーゲン・ドリル)
            // 【修正パッチ】: ここでZ方向への強い推進力を与える
            double angle = 120.0 * (PI / 180.0); // 3回で一周
            double newDirX = dirX * std::cos(angle) - dirY * std::sin(angle);
            double newDirY = dirX * std::sin(angle) + dirY * std::cos(angle);
            dirX = newDirX; dirY = newDirY;
            dirZ = 1.5; // 強力な推進力！これで「板」にはならない
        }
        else {
            // Type 0: Loop (ランダムウォーク)
            double angle = GOLDEN_ANGLE;
            double newDirX = dirX * std::cos(angle) - dirZ * std::sin(angle);
            double newDirZ = dirX * std::sin(angle) + dirZ * std::cos(angle);
            dirX = newDirX; dirZ = newDirZ;
            // カオスな動きを加える
            dirY += std::sin(i) * 0.5;
        }

        normalize(dirX, dirY, dirZ);
        currX += dirX * CA_DIST;
        currY += dirY * CA_DIST;
        currZ += dirZ * CA_DIST;

        res.x = currX;
        res.y = currY;
        res.z = currZ;
        
        protein.push_back(res);
    }

    // 3. Hydro-Zip v2 (Advanced Folding)
    // 単純な重心への圧縮ではなく、近傍の疎水性コア形成を優先する
    
    // 全体の重心
    double globX = 0, globY = 0, globZ = 0;
    for(auto& r : protein) { globX += r.x; globY += r.y; globZ += r.z; }
    globX /= protein.size(); globY /= protein.size(); globZ /= protein.size();

    for (size_t i = 0; i < protein.size(); ++i) {
        if (protein[i].structureType == 3) continue; // コラーゲンは硬いので圧縮しない！
        
        double pullFactor = protein[i].isHydrophobic ? 0.05 : 0.01; // 係数を下げてマイルドに
        
        // Helixなどは構造を保ちたいので、あまり動かさない
        if (protein[i].structureType == 1) pullFactor *= 0.5;

        protein[i].x = protein[i].x * (1.0 - pullFactor) + globX * pullFactor;
        protein[i].y = protein[i].y * (1.0 - pullFactor) + globY * pullFactor;
        protein[i].z = protein[i].z * (1.0 - pullFactor) + globZ * pullFactor;
    }

    // 4. 原子生成 & 出力 (Atomic Rendering)
    std::ofstream file("universe_os_final.pdb");
    int atomSerial = 1;

    for (size_t i = 0; i < protein.size(); ++i) {
        // バックボーン (N-CA-C-O)
        file << "ATOM  " << std::setw(5) << atomSerial++ << "  N   " 
             << std::setw(3) << protein[i].resName << " A" << std::setw(4) << protein[i].resSeq << "    "
             << std::setw(8) << std::fixed << std::setprecision(3) << protein[i].x - 0.5 
             << std::setw(8) << protein[i].y + 0.5 << std::setw(8) << protein[i].z << "  1.00  4.00           N  \n";
             
        file << "ATOM  " << std::setw(5) << atomSerial++ << "  CA  " 
             << std::setw(3) << protein[i].resName << " A" << std::setw(4) << protein[i].resSeq << "    "
             << std::setw(8) << std::fixed << std::setprecision(3) << protein[i].x 
             << std::setw(8) << protein[i].y << std::setw(8) << protein[i].z << "  1.00  4.00           C  \n";

        // 側鎖 (Side Chain) - 簡易版: 重心と逆方向へ伸ばす
        double vecX = protein[i].x - globX;
        double vecY = protein[i].y - globY;
        double vecZ = protein[i].z - globZ;
        normalize(vecX, vecY, vecZ);
        
        // 疎水性は内側(-)、親水性は外側(+)
        double dir = protein[i].isHydrophobic ? -1.0 : 1.0;
        int sideChains = getAAInfo(protein[i].code).second;
        
        for(int k=0; k<sideChains; k++) {
            double dist = 1.5 * (k+1);
            file << "ATOM  " << std::setw(5) << atomSerial++ << "  C" << (char)('B'+k) << "  " 
                 << std::setw(3) << protein[i].resName << " A" << std::setw(4) << protein[i].resSeq << "    "
                 << std::setw(8) << std::fixed << std::setprecision(3) << protein[i].x + (vecX * dist * dir)
                 << std::setw(8) << protein[i].y + (vecY * dist * dir)
                 << std::setw(8) << protein[i].z + (vecZ * dist * dir) << "  1.00  4.00           C  \n";
        }
        
        // CONECTレコード (可視化用)
        if (i > 0) {
             // 前のCAと今のCAを結ぶ（簡易接続）
             // ※本当は厳密なシリアル番号指定が必要ですが、Mol*等はCA距離で勝手に繋いでくれることが多いです
        }
    }
    file << "END" << std::endl;
    std::cout << "Successfully generated: universe_os_final.pdb" << std::endl;

    return 0;
}

