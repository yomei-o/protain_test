#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <cmath>
#include <iomanip>
#include <map>

// --- 宇宙OS 定数 ---
const double PI = 3.1415926535;
const double CA_DIST = 3.8; // アミノ酸間の基本距離 (Angstrom)

// --- 構造体 ---
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
    int structureType; // 0:Loop, 1:Helix, 2:Sheet
    double x, y, z; // CA座標
    std::vector<Atom> atoms;
};

// 3文字コード変換
std::string getResName(char code) {
    static std::map<char, std::string> m = {
        {'A', "ALA"}, {'R', "ARG"}, {'N', "ASN"}, {'D', "ASP"}, {'C', "CYS"},
        {'Q', "GLN"}, {'E', "GLU"}, {'G', "GLY"}, {'H', "HIS"}, {'I', "ILE"},
        {'L', "LEU"}, {'K', "LYS"}, {'M', "MET"}, {'F', "PHE"}, {'P', "PRO"},
        {'S', "SER"}, {'T', "THR"}, {'W', "TRP"}, {'Y', "TYR"}, {'V', "VAL"}
    };
    return m.count(code) ? m[code] : "UNK";
}

// 疎水性判定 (Hydro-Zip対象)
bool isHydrophobic(char aa) {
    return std::string("LIVFMWAC").find(aa) != std::string::npos;
}

// 簡易的な二次構造予測 (Chou-Fasmanの宇宙OS版近似)
// OSがアミノ酸の「性格」を見て、バネにするか壁にするか決める
int predictStructure(char code) {
    std::string helixFormers = "EMALKQH"; // 螺旋になりやすい
    std::string sheetFormers = "VIYFWCT"; // 壁になりやすい
    
    if (helixFormers.find(code) != std::string::npos) return 1; // Helix
    if (sheetFormers.find(code) != std::string::npos) return 2; // Sheet
    return 0; // Loop (それ以外は曲がり角)
}

// コドン変換
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

int main() {
    std::string rnaSeq;
    std::cout << "Enter RNA (Ubiquitin recommended): ";
    std::cin >> rnaSeq;

    std::vector<Residue> protein;
    
    // 1. 翻訳 & 構造タグ付け
    for (size_t i = 0; i < rnaSeq.length(); i += 3) {
        if (i+3 > rnaSeq.length()) break;
        char c = codonToAA(rnaSeq.substr(i, 3));
        if (c == '_' || c == '?') break;
        
        Residue res;
        res.resSeq = i/3 + 1;
        res.code = c;
        res.resName = getResName(c);
        res.isHydrophobic = isHydrophobic(c);
        res.structureType = predictStructure(c);
        protein.push_back(res);
    }

    // 2. 「生きた形」の生成 (Procedural Folding)
    // カメ（Turtle Graphics）のように、先頭から順に空間を泳がせる
    
    double currX = 0, currY = 0, currZ = 0;
    // 進行方向ベクトル
    double dirX = 1.0, dirY = 0.0, dirZ = 0.0;
    
    for (size_t i = 0; i < protein.size(); ++i) {
        
        // --- 局所ルール (Local Rule) ---
        // 構造タイプによって動き方を変える
        
        if (protein[i].structureType == 1) { 
            // Alpha Helix: 螺旋状に進む
            // 進行方向に対して少しずつ回転を加える
            double angle = 100.0 * (PI / 180.0); // 1残基あたり100度回転
            double newDirX = dirX * std::cos(angle) - dirZ * std::sin(angle);
            double newDirZ = dirX * std::sin(angle) + dirZ * std::cos(angle);
            dirX = newDirX; dirZ = newDirZ;
            // 縦方向（Y）にも進む（バネの伸び）
            dirY += 0.2; 
        } 
        else if (protein[i].structureType == 2) {
            // Beta Sheet: ジグザグに進む
            // Y軸方向へ交互に振れる
            dirY = (i % 2 == 0) ? 0.8 : -0.8;
            // XZ平面では少し曲がる（シートのねじれ）
            double angle = 15.0 * (PI / 180.0);
            double newDirX = dirX * std::cos(angle) - dirZ * std::sin(angle);
            double newDirZ = dirX * std::sin(angle) + dirZ * std::cos(angle);
            dirX = newDirX; dirZ = newDirZ;
        }
        else {
            // Loop: ランダムな方向転換（グリシンやプロリンによるターン）
            // ここで大きく形が崩れる＝「不格好」になる要素
            double angle = 137.5 * (PI / 180.0); // 黄金角で大きく曲がる
            double newDirX = dirX * std::cos(angle) - dirY * std::sin(angle);
            double newDirY = dirX * std::sin(angle) + dirY * std::cos(angle);
            dirX = newDirX; dirY = newDirY;
            
            // ループ部分では疎水性の引力が強く働き、中心へ向かおうとするベクトルを足す
            dirX -= currX * 0.1;
            dirY -= currY * 0.1;
            dirZ -= currZ * 0.1;
        }

        // ベクトルの正規化
        double len = std::sqrt(dirX*dirX + dirY*dirY + dirZ*dirZ);
        dirX /= len; dirY /= len; dirZ /= len;

        // 座標更新 (CAの位置決定)
        currX += dirX * CA_DIST;
        currY += dirY * CA_DIST;
        currZ += dirZ * CA_DIST;

        protein[i].x = currX;
        protein[i].y = currY;
        protein[i].z = currZ;
    }

    // 3. Hydro-Zip (疎水性崩壊) の適用
    // 生成された「ヒモ」はまだ伸びているので、疎水性コアを中心に「グシャッ」と潰す
    
    double centerX = 0, centerY = 0, centerZ = 0;
    int hydrophobicCount = 0;
    for (const auto& r : protein) {
        if (r.isHydrophobic) {
            centerX += r.x; centerY += r.y; centerZ += r.z;
            hydrophobicCount++;
        }
    }
    if (hydrophobicCount > 0) {
        centerX /= hydrophobicCount; centerY /= hydrophobicCount; centerZ /= hydrophobicCount;
    }

    // 重心に向かって引き寄せる（親水性はあまり動かさない）
    for (auto& r : protein) {
        double pullStrength = r.isHydrophobic ? 0.7 : 0.3; // 疎水性は強く引く
        r.x = r.x * (1.0 - pullStrength) + centerX * pullStrength;
        r.y = r.y * (1.0 - pullStrength) + centerY * pullStrength;
        r.z = r.z * (1.0 - pullStrength) + centerZ * pullStrength;
        
        // 衝突回避（簡易）：近づきすぎたら少し離す（パウリの排他律）
        double distToCenter = std::sqrt((r.x-centerX)*(r.x-centerX) + (r.y-centerY)*(r.y-centerY) + (r.z-centerZ)*(r.z-centerZ));
        if (distToCenter < 5.0) {
            r.x *= 1.2; r.y *= 1.2; r.z *= 1.2;
        }
    }

    // --- PDB出力 (CAのみ。肉付けはMol* Viewerに任せる) ---
    std::ofstream file("universe_os_living_shape.pdb");
    int atomSerial = 1;
    for (size_t i = 0; i < protein.size(); ++i) {
        file << "ATOM  " 
             << std::setw(5) << atomSerial++ << "  CA  " 
             << std::setw(3) << protein[i].resName << " A" 
             << std::setw(4) << protein[i].resSeq << "    "
             << std::setw(8) << std::fixed << std::setprecision(3) << protein[i].x
             << std::setw(8) << std::fixed << std::setprecision(3) << protein[i].y
             << std::setw(8) << std::fixed << std::setprecision(3) << protein[i].z
             << "  1.00  4.00           C  " << std::endl;
             
        if (i > 0) {
             file << "CONECT" 
                 << std::setw(5) << (atomSerial - 1) 
                 << std::setw(5) << (atomSerial - 2) << std::endl;
        }
    }
    file << "END" << std::endl;
    std::cout << "Living shape generated: universe_os_living_shape.pdb" << std::endl;
    return 0;
}


