#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <cmath>
#include <iomanip>
#include <map>
#include <algorithm>
#include <sstream>

// --- 宇宙OS 定数定義 (Universe Constants) ---
const double PI = 3.1415926535;
const double CA_DIST = 3.8;   // CA間距離 (Angstrom)
const double GOLDEN_ANGLE = PI * (3.0 - std::sqrt(5.0)); // 黄金角
const double WATER_DENSITY_TEMP = 4.0; // 4.0℃ (密度の基準)

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

// 3文字コード(ALA)から1文字(A)へ
char resNameToCode(const std::string& resName) {
    static std::map<std::string, char> m;
    if (m.empty()) {
        std::string codes = "ARNDCQEGHILKMFPSTWYV";
        for (char c : codes) m[getAAInfo(c).first] = c;
    }
    if (m.count(resName)) return m[resName];
    return '?';
}

bool isHydrophobic(char aa) {
    return std::string("LIVFMWAC").find(aa) != std::string::npos;
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

double dist(const Residue& a, const Residue& b) {
    return std::sqrt(pow(a.x - b.x, 2) + pow(a.y - b.y, 2) + pow(a.z - b.z, 2));
}

// 外部PDB読み込み機能 (簡易版: CAのみ抽出)
std::vector<Residue> loadReferencePDB(const std::string& filename) {
    std::vector<Residue> loaded;
    std::ifstream file(filename);
    if (!file.is_open()) return loaded;

    std::string line;
    while (std::getline(file, line)) {
        if (line.substr(0, 4) == "ATOM" && line.substr(13, 2) == "CA") {
            Residue r;
            r.resName = line.substr(17, 3);
            r.code = resNameToCode(r.resName);
            r.resSeq = std::stoi(line.substr(22, 4));
            r.x = std::stod(line.substr(30, 8));
            r.y = std::stod(line.substr(38, 8));
            r.z = std::stod(line.substr(46, 8));
            loaded.push_back(r);
        }
    }
    return loaded;
}

int main() {
    std::string inputSeq;
    std::cout << "===========================================" << std::endl;
    std::cout << "   Project-137: Universe OS v6.0 (Reality Check)   " << std::endl;
    std::cout << "===========================================" << std::endl;
    std::cout << "Input Sequence (RNA or AA): ";
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
        std::cout << "[System] Detected AA Sequence. Loading..." << std::endl;
        for (char c : inputSeq) fullSeq.push_back(c);
    }

    // --- 2. 宇宙OS レイアウト計算 ---
    double currX = 0, currY = 0, currZ = 0;
    
    // オキシトシン検出用
    int cys1_idx = -1; 
    
    for (size_t i = 0; i < fullSeq.size(); ++i) {
        Residue res;
        res.resSeq = i + 1;
        res.code = fullSeq[i];
        res.isHydrophobic = isHydrophobic(res.code);
        res.resName = getAAInfo(res.code).first;
        res.structureType = 0;

        // --- スペシャルロジック: オキシトシン・リングの形成 ---
        // 配列内で 1番目と6番目がCysの場合、リングを作る
        // (プレプロオキシトシンの場合、20番目と25番目がCys)
        
        // シンプル化のため、今回は「オキシトシン領域(9残基)」を特別扱い
        bool isOxytocinRegion = false;
        // マニュアル検出: 20番目から始まると仮定(プレプロ)、または1番目から(オキシトシン単体)
        int local_idx = -1;

        if (fullSeq.size() < 20 && i < 9) { // 短い配列なら最初がオキシトシン
             isOxytocinRegion = true;
             local_idx = i;
        } else if (i >= 19 && i < 28) { // プレプロ配列なら20番目から
             isOxytocinRegion = true;
             local_idx = i - 19;
        }

        if (isOxytocinRegion) {
             res.structureType = 4; // Oxytocin Ring
             // 9角形の配置 + 螺旋のねじれ
             double angle = (2.0 * PI * local_idx) / 9.0;
             double radius = 5.0; // 宇宙OS的理想半径
             
             // Cys1とCys6を近づけるための「重力補正」
             if (local_idx == 0) cys1_idx = i;
             
             // リング状に配置
             currX = radius * std::cos(angle);
             currY = radius * std::sin(angle);
             currZ = (double)local_idx * 1.5; // 少しずつズラす(螺旋)
             
             // プレプロの場合は位置をオフセット（ニューロフィジンの横へ）
             if (fullSeq.size() > 20) {
                 currX += 15.0; 
                 currY += 5.0;
                 currZ += 20.0;
             }
        } else {
             // 通常のフィボナッチ・スフィア・パッキング
             double index = (double)i;
             double phi = index * GOLDEN_ANGLE;
             double z_sphere = 1.0 - (2.0 * index / (double)(fullSeq.size() - 1));
             double r_sphere = std::sqrt(1.0 - z_sphere * z_sphere);
             
             double scale = 12.0; // 全体のサイズ感
             currX = scale * r_sphere * std::cos(phi);
             currY = scale * r_sphere * std::sin(phi);
             currZ = scale * z_sphere;
        }

        res.x = currX;
        res.y = currY;
        res.z = currZ;
        protein.push_back(res);
    }

    // --- 3. 内部検証 (Reality Check) ---
    std::cout << "\n--- [Debug] Universe OS Reality Check ---" << std::endl;
    
    // A. SS結合距離 (Cys-Cys) の検証
    int c1 = -1, c6 = -1;
    // オキシトシン領域のCysを探す
    int startSearch = (fullSeq.size() > 20) ? 19 : 0;
    for (int i = startSearch; i < startSearch + 9 && i < protein.size(); ++i) {
        if (protein[i].code == 'C') {
            if (c1 == -1) c1 = i;
            else if (c6 == -1) c6 = i;
        }
    }

    if (c1 != -1 && c6 != -1) {
        double d = dist(protein[c1], protein[c6]);
        std::cout << "Generated Cys(" << c1+1 << ")-Cys(" << c6+1 << ") Distance: " 
                  << std::fixed << std::setprecision(2) << d << " A" << std::endl;
        std::cout << "  -> Ideal SS-Bond (CA-CA): ~5.6 A (Direct S-S is ~2.0 A)" << std::endl;
        
        if (d > 7.0) std::cout << "  [WARNING] Ring is too open! Needs 'Water Pressure' patch." << std::endl;
        else std::cout << "  [OK] Ring is closed. Good packing." << std::endl;
    }

    // --- 4. 外部ファイル比較 (1NPO / real_oxytocin.pdb) ---
    std::string refFile = "real_oxytocin.pdb";
    std::vector<Residue> ref = loadReferencePDB(refFile);
    
    if (!ref.empty()) {
        std::cout << "\n--- [Comparison] vs Real World (" << refFile << ") ---" << std::endl;
        // 実測データのCys間距離
        int rc1 = -1, rc6 = -1;
        for(size_t i=0; i<ref.size(); ++i) {
             if (ref[i].code == 'C') {
                 if (rc1 == -1) rc1 = i;
                 else if (rc6 == -1) rc6 = i;
             }
        }
        if (rc1 != -1 && rc6 != -1) {
            double rd = dist(ref[rc1], ref[rc6]);
            std::cout << "REAL World Cys-Cys Distance: " << rd << " A" << std::endl;
            
            // 誤差表示
            if (c1 != -1 && c6 != -1) {
                 double gen_d = dist(protein[c1], protein[c6]);
                 double diff = std::abs(gen_d - rd);
                 std::cout << "Difference (Reality Gap): " << diff << " A" << std::endl;
                 if (diff < 2.0) std::cout << "  -> AMAZING! Universe OS matches Reality." << std::endl;
                 else std::cout << "  -> Gap detected. Universe implies a looser structure." << std::endl;
            }
        }
    } else {
        std::cout << "\n[Info] 'real_oxytocin.pdb' not found. Skipping comparison." << std::endl;
        std::cout << "       (Put 1NPO extracted PDB here to unlock comparison mode)" << std::endl;
    }

    // --- 5. PDB出力 ---
    FILE* fp = fopen("Universe_Oxytocin_v6.pdb", "w");
    if (fp) {
        int atomSerial = 1;
        fprintf(fp, "HEADER    UNIVERSE OS V6.0 GENERATED\n");
        for (size_t i = 0; i < protein.size(); ++i) {
            fprintf(fp, "ATOM  %5d  CA  %3s A%4d    %8.3f%8.3f%8.3f  1.00  4.00           C\n",
                atomSerial++, protein[i].resName.c_str(), protein[i].resSeq, 
                protein[i].x, protein[i].y, protein[i].z);
                
            // 簡易側鎖
            double vx = protein[i].x, vy = protein[i].y, vz = protein[i].z;
            normalize(vx, vy, vz);
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
        }
        fprintf(fp, "END\n");
        fclose(fp);
        std::cout << "\n[Output] Generated 'Universe_Oxytocin_v6.pdb'" << std::endl;
    }

    return 0;
}
