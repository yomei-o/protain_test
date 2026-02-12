#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <cmath>
#include <iomanip>
#include <map>
#include <algorithm>

// --- 宇宙OS 定数定義 (Universe Constants) ---
const double PI = 3.1415926535;
const double CA_DIST = 3.8;   // CA間距離 (Angstrom)
const double GOLDEN_ANGLE = PI * (3.0 - std::sqrt(5.0)); // 黄金角
const double WATER_DENSITY_TEMP = 4.0; // 4.0℃

// --- データ構造 ---
struct Atom {
    int serial;
    std::string name;
    double x, y, z;
};

struct Residue {
    int resSeq;
    char code;
    std::string resName;
    bool isHydrophobic;
    int structureType; // 0:Loop, 1:Helix, 2:Sheet, 3:Collagen, 4:Oxytocin-Ring
    double x, y, z;
};

// --- ヘルパー関数 ---

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

// RNA -> アミノ酸 変換 (コドン表)
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

// 入力がRNAかアミノ酸かを判定
bool isRNASequence(const std::string& seq) {
    for(char c : seq) {
        if (c != 'A' && c != 'U' && c != 'G' && c != 'C' && c != 'T') return false;
    }
    return true;
}

void normalize(double& x, double& y, double& z) {
    double len = std::sqrt(x*x + y*y + z*z);
    if (len > 1e-6) { x /= len; y /= len; z /= len; }
}

int main() {
    std::string inputSeq;
    std::cout << "===========================================" << std::endl;
    std::cout << "   Project-137: Universe OS Folder v5.0    " << std::endl;
    std::cout << "===========================================" << std::endl;
    std::cout << "Input Sequence (RNA or AminoAcid): ";
    std::cin >> inputSeq;

    std::vector<Residue> protein;
    std::vector<char> fullSeq;

    // --- 1. 自動判別と翻訳 ---
    if (isRNASequence(inputSeq)) {
        std::cout << "[System] Detected RNA. Translating..." << std::endl;
        for (size_t i = 0; i < inputSeq.length(); i += 3) {
            if (i+3 > inputSeq.length()) break;
            char c = codonToAA(inputSeq.substr(i, 3));
            if (c == '_' || c == '?') break;
            fullSeq.push_back(c);
        }
    } else {
        std::cout << "[System] Detected Amino Acid Sequence. Loading directly..." << std::endl;
        for (char c : inputSeq) fullSeq.push_back(c);
    }

    // --- 2. 宇宙OS レイアウト計算 ---
    double currX = 0, currY = 0, currZ = 0;
    // 初期進行方向
    double dirX = 1.0, dirY = 0.0, dirZ = 0.0; 

    // ニューロフィジンのための球形パッキング用変数
    double sphereRadius = 10.0; 
    double phi = 0.0, theta = 0.0;

    for (size_t i = 0; i < fullSeq.size(); ++i) {
        Residue res;
        res.resSeq = i + 1;
        res.code = fullSeq[i];
        res.isHydrophobic = isHydrophobic(res.code);
        auto info = getAAInfo(res.code);
        res.resName = info.first;
        res.structureType = 0;

        // --- 特殊ロジック: プレプロ・オキシトシンの領域判定 ---
        
        // A. シグナルペプチド (1-19): 外へ向かって直線的に伸びる
        if (i < 19) {
            dirX = 0.2; dirY = 0.2; dirZ = 1.0; // Z方向へ排出
            normalize(dirX, dirY, dirZ);
            res.structureType = 1; // Helix扱い
        }
        // B. オキシトシン本体 (20-28): 小さなループ（指輪）を作る
        else if (i >= 19 && i < 28) {
            // 9残基で一周するような急激なカーブ
            double oxAngle = (2.0 * PI) / 9.0;
            double dx = std::cos(oxAngle * (i-19));
            double dy = std::sin(oxAngle * (i-19));
            
            // 少しオフセットして配置
            currX = 15.0 + (5.0 * dx); // X=15あたりに浮かべる
            currY = 5.0 + (5.0 * dy);
            currZ = 20.0; // シグナルの先
            
            res.structureType = 4; // Special Ring
            // 座標を直接指定したので移動加算はスキップ
            res.x = currX; res.y = currY; res.z = currZ;
            protein.push_back(res);
            continue; 
        }
        // C. ニューロフィジン (32以降): 巨大な球状パッキング（ケース）
        else if (i >= 31) {
            // フィボナッチ螺旋による球体配置 ($D=\pi$ Packing)
            // 原点付近に大きな塊を作る
            double index = (double)(i - 31);
            double y_pos = 1.0 - (index / (double)(fullSeq.size() - 31 - 1)) * 2.0; // 1 to -1
            double radius_at_y = std::sqrt(1.0 - y_pos * y_pos);
            
            double theta_golden = index * GOLDEN_ANGLE;
            
            double r = 12.0; // ケースの大きさ
            
            currX = r * std::cos(theta_golden) * radius_at_y;
            currY = r * y_pos;
            currZ = r * std::sin(theta_golden) * radius_at_y;
            
            res.structureType = 2; // Sheet-like packing
            res.x = currX; res.y = currY; res.z = currZ;
            protein.push_back(res);
            continue;
        }
        // D. リンカー (29-31): つなぎ目
        else {
            dirX = 0.5; dirY = 0.5; dirZ = 0.0;
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

    // --- 3. PDB出力 (Strict Formatting) ---
    // PDBのフォーマットは桁数が厳密なので、printfで整列させます
    FILE* fp = fopen("Prepro_Oxytocin_v5.pdb", "w");
    if (!fp) {
        std::cerr << "Error: Cannot open output file." << std::endl;
        return 1;
    }

    int atomSerial = 1;
    fprintf(fp, "HEADER    UNIVERSE OS GENERATED PDB\n");
    fprintf(fp, "REMARK    GENERATED BY PROJECT-137 V5.0\n");

    for (size_t i = 0; i < protein.size(); ++i) {
        // Backbone N
        fprintf(fp, "ATOM  %5d  N   %3s A%4d    %8.3f%8.3f%8.3f  1.00  4.00           N\n",
            atomSerial++, protein[i].resName.c_str(), protein[i].resSeq, 
            protein[i].x - 0.5, protein[i].y - 0.5, protein[i].z);
            
        // Backbone CA
        fprintf(fp, "ATOM  %5d  CA  %3s A%4d    %8.3f%8.3f%8.3f  1.00  4.00           C\n",
            atomSerial++, protein[i].resName.c_str(), protein[i].resSeq, 
            protein[i].x, protein[i].y, protein[i].z);

        // Backbone C
        fprintf(fp, "ATOM  %5d  C   %3s A%4d    %8.3f%8.3f%8.3f  1.00  4.00           C\n",
            atomSerial++, protein[i].resName.c_str(), protein[i].resSeq, 
            protein[i].x + 0.5, protein[i].y + 0.5, protein[i].z);

        // 簡易側鎖 (方向は原点逆向き)
        double vx = protein[i].x, vy = protein[i].y, vz = protein[i].z;
        normalize(vx, vy, vz);
        // 疎水性は内側へ
        double sideDir = protein[i].isHydrophobic ? -1.0 : 1.0;
        
        int scCount = getAAInfo(protein[i].code).second;
        for(int k=0; k<scCount; k++) {
            double d = 1.5 * (k+1);
             fprintf(fp, "ATOM  %5d  C%c  %3s A%4d    %8.3f%8.3f%8.3f  1.00  4.00           C\n",
                atomSerial++, 'B'+k, protein[i].resName.c_str(), protein[i].resSeq, 
                protein[i].x + (vx * d * sideDir), 
                protein[i].y + (vy * d * sideDir), 
                protein[i].z + (vz * d * sideDir));
        }
        
        // Connect logic simply (i to i+1)
        if (i < protein.size() - 1) {
            // PDB CONECT record is complex, usually viewers auto-connect by distance.
            // Skipping strictly for simplicity, but viewers will handle CA trace.
        }
    }
    fprintf(fp, "END\n");
    fclose(fp);

    std::cout << "-------------------------------------------" << std::endl;
    std::cout << "SUCCESS: 'Prepro_Oxytocin_v5.pdb' generated." << std::endl;
    std::cout << "Total Amino Acids: " << protein.size() << std::endl;
    std::cout << "Structure Info: " << std::endl;
    std::cout << "  - Signal Peptide (1-19): Linear Extraction" << std::endl;
    std::cout << "  - Oxytocin (20-28): Ring Formation" << std::endl;
    std::cout << "  - Neurophysin (32-): Globular Packing" << std::endl;
    std::cout << "-------------------------------------------" << std::endl;

    return 0;
}

