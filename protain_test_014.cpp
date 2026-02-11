#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <cmath>
#include <iomanip>
#include <map>

// --- Project-137: Universe OS v9.1 "Pattern-Folder" ---
// コンセプト: 10歩という固定概念を捨て、アミノ酸の「性格（PやG）」で折り返しを決定する

const double PI = 3.1415926535;
const double CA_DIST = 3.8; 

struct Residue {
    int resSeq;
    char code;
    std::string resName;
    double x, y, z;
};

std::pair<std::string, int> getAAInfo(char code) {
    static std::map<char, std::pair<std::string, int>> m = {
        {'G', {"GLY", 0}}, {'A', {"ALA", 1}}, {'P', {"PRO", 2}}, {'V', {"VAL", 2}},
        {'S', {"SER", 1}}, {'T', {"THR", 2}}, {'L', {"LEU", 3}}, {'I', {"ILE", 3}}
    };
    return m.count(code) ? m[code] : std::make_pair("UNK", 0);
}

int main() {
    std::string inputSeq;
    std::cout << "Project-137 v9.1 [Pattern-Folder Mode]" << std::endl;
    // テスト推奨配列: GAGAGAPGAGAGA (Pで折り返します)
    std::cout << "Input Sequence: ";
    std::cin >> inputSeq;

    std::vector<Residue> protein;
    double currX = 0, currY = 0, currZ = 0;
    
    // 状態管理変数
    double zDirection = 1.0; // 1.0 = 往路, -1.0 = 復路
    double xOffset = 0.0;    // 折り返すたびに横にずらす
    bool turning = false;

    for (int i = 0; i < (int)inputSeq.length(); ++i) {
        Residue res;
        res.resSeq = i + 1;
        res.code = inputSeq[i];
        res.resName = getAAInfo(res.code).first;

        // --- 【Pattern-Logic: 転回シグナルの検知】 ---
        // プロリン(P)を見つけるか、グリシン(G)が連続して遊びが生まれた時に曲がる
        if (!turning && (res.code == 'P' || (i > 0 && res.code == 'G' && inputSeq[i-1] == 'G'))) {
            turning = true;
            zDirection *= -1.0; // 進行方向を反転
            xOffset += 4.5;     // 隣の列に移動 (水素結合が作れる距離)
            std::cout << "[System] Turn Signal detected at residue " << i+1 << " (" << res.code << ")" << std::endl;
        } else {
            turning = false;
        }

        // --- 【Weaver Engine 2.0】 ---
        double dirX, dirY, dirZ;
        
        // 基本はジグザグに進む
        double zigzag = (i % 2 == 0) ? 0.7 : -0.7;
        
        dirX = zigzag; 
        dirY = 0.1; // わずかな歪み
        dirZ = zDirection * 1.5; // 往路か復路か

        // 正規化
        double len = std::sqrt(dirX*dirX + dirY*dirY + dirZ*dirZ);
        dirX /= len; dirY /= len; dirZ /= len;

        currX += dirX * CA_DIST;
        currY += dirY * CA_DIST;
        currZ += dirZ * CA_DIST;

        // 計算した座標にオフセットを適用
        res.x = currX + xOffset; 
        res.y = currY; 
        res.z = currZ;
        
        protein.push_back(res);
    }

    std::ofstream out("Universe_v9_1_Folder.pdb");
    out << "HEADER    PROJECT-137 V9.1 PATTERN FOLDER\n";
    for (auto& r : protein) {
        char buf[100];
        sprintf(buf, "ATOM  %5d  CA  %3s A%4d    %8.3f%8.3f%8.3f  1.00  4.00           C\n",
                r.resSeq, r.resName.c_str(), r.resSeq, r.x, r.y, r.z);
        out << buf;
    }
    out << "END\n";
    
    std::cout << "[Success] Generated Universe_v9_1_Folder.pdb" << std::endl;
    return 0;
}


