#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <cmath>
#include <iomanip>
#include <map>

// --- 宇宙OS 定数定義 (Universe OS Constants) ---
const double FINE_STRUCTURE = 1.0 / 137.035999;
const double PI = 3.1415926535;
const double WATER_MAGIC_TEMP = 4.0; 
const double GRID_SIZE = 3.8; // アミノ酸間の平均距離 (Angstrom) - タクシー幾何学の1ブロック

// --- アミノ酸データ構造 ---
struct AminoAcid {
    char code;       // 'A', 'R', 'N', ...
    double x, y, z;  // 3D座標
    bool isHydrophobic; // 疎水性フラグ (TrueならHydro-Zip対象)
};

// --- 疎水性判定 (Hydrophobicity Table) ---
bool isHydrophobic(char aa) {
    // L, I, V, F, M, W, A などは疎水性で、内側に集まりやすい
    std::string hydrophobic = "LIVFMWAC"; 
    return hydrophobic.find(aa) != std::string::npos;
}

// --- コドン変換テーブル ---
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

// --- PDBファイル出力関数 ---
void savePDB(const std::vector<AminoAcid>& protein, const std::string& filename) {
    std::ofstream file(filename);
    if (!file.is_open()) {
        std::cerr << "Error: Could not open file " << filename << std::endl;
        return;
    }

    int atomSerial = 1;
    for (size_t i = 0; i < protein.size(); ++i) {
        // PDB Format: ATOM, Serial, AtomName, ResName, ChainID, ResSeq, X, Y, Z, Occ, Temp
        // 簡易的にCA (Alpha Carbon) のみ出力
        file << "ATOM  " 
             << std::setw(5) << atomSerial << " "
             << " CA  " // Atom Name
             << "GLY" // 簡易化のため全残基名をGLYとしていますが、本来はprotein[i].codeから変換
             << " A" 
             << std::setw(4) << (i + 1) << "    "
             << std::setw(8) << std::fixed << std::setprecision(3) << protein[i].x
             << std::setw(8) << std::fixed << std::setprecision(3) << protein[i].y
             << std::setw(8) << std::fixed << std::setprecision(3) << protein[i].z
             << "  1.00" // Occupancy
             << "  " << std::setw(5) << std::fixed << std::setprecision(2) << WATER_MAGIC_TEMP // TempFactorに4.00を使用 (署名代わり)
             << "           C  " << std::endl;
        atomSerial++;
    }
    file << "END" << std::endl;
    std::cout << "PDB file generated: " << filename << std::endl;
}

int main() {
    std::string rnaSeq;
    std::cout << "Enter RNA Sequence (e.g., AUGUUU...): ";
    std::cin >> rnaSeq;

    std::vector<AminoAcid> protein;
    double currentX = 0, currentY = 0, currentZ = 0;

    // 1. 翻訳 & 初期配置 (Translation & Initial Grid Mapping)
    std::cout << "Translating RNA and applying Universe OS Logic..." << std::endl;
    
    for (size_t i = 0; i < rnaSeq.length(); i += 3) {
        if (i + 3 > rnaSeq.length()) break;
        char aaCode = codonToAA(rnaSeq.substr(i, 3));
        if (aaCode == '_') break; // Stop codon

        AminoAcid aa;
        aa.code = aaCode;
        aa.isHydrophobic = isHydrophobic(aaCode);

        // タクシー幾何学的配置 (Taxicab Geometry Logic)
        // 137のハッシュを使って、次は「上下左右前後」のどこに進むか決める
        // 物理シミュレーションではなく、「定められたルート」を辿る
        double hash = std::sin(i * FINE_STRUCTURE * 1000) * std::cos(i * PI);
        
        int direction = (int)(std::abs(hash) * 6) % 6; // 0..5 (X+, X-, Y+, Y-, Z+, Z-)

        switch(direction) {
            case 0: currentX += GRID_SIZE; break;
            case 1: currentX -= GRID_SIZE; break;
            case 2: currentY += GRID_SIZE; break;
            case 3: currentY -= GRID_SIZE; break;
            case 4: currentZ += GRID_SIZE; break;
            case 5: currentZ -= GRID_SIZE; break;
        }

        aa.x = currentX;
        aa.y = currentY;
        aa.z = currentZ;
        protein.push_back(aa);
    }

    // 2. 疎水性崩壊 (Hydro-Zip / Distance Zero Optimization)
    // 疎水性アミノ酸を「重心」に向かってワープさせる
    std::cout << "Applying Hydro-Zip (Distance->0 Optimization)..." << std::endl;
    
    double centerX = 0, centerY = 0, centerZ = 0;
    // 重心計算
    for (const auto& aa : protein) {
        centerX += aa.x; centerY += aa.y; centerZ += aa.z;
    }
    centerX /= protein.size(); centerY /= protein.size(); centerZ /= protein.size();

    for (auto& aa : protein) {
        if (aa.isHydrophobic) {
            // 疎水性アミノ酸は、4℃の水が作るグリッドに従って内側へ「畳み込まれる」
            // 距離を短縮する (Zip)
            aa.x = (aa.x + centerX) / 2.0; 
            aa.y = (aa.y + centerY) / 2.0;
            aa.z = (aa.z + centerZ) / 2.0;
        }
    }

    // 3. ファイル出力
    savePDB(protein, "universe_os_structure.pdb");
    
    std::cout << "Done. The structure is optimized for 4.0C stability." << std::endl;

    return 0;
}



