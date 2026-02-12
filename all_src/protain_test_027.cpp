#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <cmath>
#include <iomanip>
#include <map>
#include <algorithm>

// --- Project-137: Universe OS v16.1 "Collagen Wire Patch" ---
// コンセプト: プロリンの「剛性（硬さ）」を物理エンジンに実装。
// プロリンが密集する領域（コラーゲン）は、丸まらずに強靭なワイヤーになる。

const double PI = 3.1415926535;
const double CA_DIST = 3.8;
const double MIN_DIST = 4.0;

struct Vec3 {
    double x, y, z;
    Vec3 operator-(const Vec3& o) const { return {x - o.x, y - o.y, z - o.z}; }
    Vec3 operator+(const Vec3& o) const { return {x + o.x, y + o.y, z + o.z}; }
    Vec3 operator*(double s) const { return {x * s, y * s, z * s}; }
};

Vec3 cross(Vec3 a, Vec3 b) { return {a.y*b.z - a.z*b.y, a.z*b.x - a.x*b.z, a.x*b.y - a.y*b.x}; }
double length(Vec3 v) { return std::sqrt(v.x*v.x + v.y*v.y + v.z*v.z); }
Vec3 normalize(Vec3 v) { double l = length(v); return (l>1e-6)? Vec3{v.x/l, v.y/l, v.z/l} : Vec3{0,0,0}; }

struct Residue {
    int resSeq; char code; std::string resName;
    double hydrophobicity; double charge; double p_alpha;
    bool isHelix; bool isWire; // 【新フラグ】ワイヤー(剛体)かどうか
    Vec3 ca, n, c, o;
};

struct AAProps { std::string name; double hydro; double charge; double p_alpha; };
AAProps getAAProps(char code) {
    static std::map<char, AAProps> m = {
        {'E',{"GLU",-3.5,-1, 1.51}}, {'M',{"MET", 1.9, 0, 1.45}}, {'A',{"ALA", 1.8, 0, 1.42}},
        {'L',{"LEU", 3.8, 0, 1.21}}, {'K',{"LYS",-3.9, 1, 1.16}}, {'F',{"PHE", 2.8, 0, 1.13}},
        {'Q',{"GLN",-3.5, 0, 1.11}}, {'W',{"TRP",-0.9, 0, 1.08}}, {'I',{"ILE", 4.5, 0, 1.08}},
        {'V',{"VAL", 4.2, 0, 1.06}}, {'D',{"ASP",-3.5,-1, 1.01}}, {'H',{"HIS",-3.2, 1, 1.00}},
        {'R',{"ARG",-4.5, 1, 0.98}}, {'T',{"THR",-0.7, 0, 0.83}}, {'S',{"SER",-0.8, 0, 0.77}},
        {'C',{"CYS", 2.5, 0, 0.70}}, {'Y',{"TYR",-1.3, 0, 0.69}}, {'N',{"ASN",-3.5, 0, 0.67}},
        {'P',{"PRO",-1.6, 0, 0.57}}, {'G',{"GLY",-0.4, 0, 0.57}}
    };
    if(m.count(code)) return m[code];
    return {"UNK",0,0, 0.5};
}

// 翻訳エンジン省略（アミノ酸直接入力）

int main() {
    std::string rawInput;
    std::cout << "--- Project-137 v16.1 [Collagen Wire Engine] ---" << std::endl;
    std::cout << "Input Amino Acid Sequence: ";
    std::cin >> rawInput;
    std::transform(rawInput.begin(), rawInput.end(), rawInput.begin(), ::toupper);

    std::vector<Residue> protein;
    for(int i=0; i<(int)rawInput.length(); ++i) {
        Residue res; res.resSeq = i+1; res.code = rawInput[i];
        AAProps p = getAAProps(res.code);
        res.resName = p.name; res.hydrophobicity = p.hydro; 
        res.charge = p.charge; res.p_alpha = p.p_alpha;
        res.isHelix = false; res.isWire = false;
        protein.push_back(res);
    }

    double currX=0, currY=0, currZ=0;
    double angleX=0.8, angleY=0.0;

    for(int i=0; i<(int)protein.size(); ++i) {
        double avgProp = 0; 
        double proCount = 0; // プロリンの数を数える
        int count=0;
        
        for(int k=i-2; k<=i+2; k++) {
            if(k>=0 && k<(int)protein.size()) { 
                avgProp += protein[k].p_alpha; 
                if (protein[k].code == 'P') proCount++; // プロリン発見！
                count++; 
            }
        }
        avgProp /= (count>0?count:1);
        double proRatio = proCount / count; // プロリンの密集度

        double localChargeSum = 0;
        for(int k=i-2; k<=i+2; k++) {
            if(k>=0 && k<(int)protein.size()) localChargeSum += std::abs(protein[k].charge);
        }

        // --- 究極の物理法則分岐 ---

        // 【新法則: コラーゲンモード】プロリンが密集(40%以上)しているなら「剛性ワイヤー」になる！
        if (proRatio >= 0.4) {
            protein[i].isWire = true;
            angleX = 0.4; // 螺旋よりも大きくZ軸方向へ真っ直ぐ伸びる
            angleY += 2.094; // コラーゲン特有の「約120度の緩やかな旋回」
        } 
        // 【アルファヘリックスモード】油っぽくて、かつPとGではない
        else if (avgProp > 1.05 && protein[i].code != 'P' && protein[i].code != 'G') {
            protein[i].isHelix = true;
            angleX = 0.9;
            angleY += 1.74533; // 完璧な螺旋(100度)
        } 
        // 【ランダムコイルモード】水っぽい部分は曲がる
        else {
            double repulsion = 1.0 + (localChargeSum * 0.2);
            angleX += (PI*(3.0-std::sqrt(5.0))) * 0.5 * repulsion;
            angleY += 0.5 * repulsion;
        }

        Vec3 dir = {std::sin(angleX)*std::cos(angleY), std::sin(angleX)*std::sin(angleY), std::cos(angleX)};
        dir = normalize(dir);
        currX += dir.x*CA_DIST; currY += dir.y*CA_DIST; currZ += dir.z*CA_DIST;
        protein[i].ca = {currX, currY, currZ};
    }

    // 剛体ロック付きリラクゼーション
    for(int iter=0; iter<200; ++iter) {
        for(int i=0; i<(int)protein.size(); ++i) {
            for(int j=i+1; j<(int)protein.size(); ++j) {
                // アルファヘリックス、または「ワイヤー」の部分は形を崩さない！
                bool rigidI = (protein[i].isHelix || protein[i].isWire);
                bool rigidJ = (protein[j].isHelix || protein[j].isWire);
                if(rigidI && rigidJ && std::abs(i - j) < 6) continue;

                Vec3 d = protein[j].ca - protein[i].ca;
                double dist = length(d);
                if(dist < MIN_DIST && dist > 0.1) {
                    double push = (MIN_DIST - dist) * 0.5 * 0.1;
                    Vec3 move = normalize(d) * push;
                    if(!rigidI) protein[i].ca = protein[i].ca - move;
                    if(!rigidJ) protein[j].ca = protein[j].ca + move;
                }
            }
        }
    }

    // フルアトム計算(ペプチド平面)
    for(int i=0; i<(int)protein.size(); ++i) {
        Vec3 ca_prev = (i > 0) ? protein[i-1].ca : protein[i].ca - Vec3{1,0,0};
        Vec3 ca_next = (i < (int)protein.size()-1) ? protein[i+1].ca : protein[i].ca + Vec3{1,0,0};

        Vec3 v_prev = normalize(protein[i].ca - ca_prev);
        Vec3 v_next = normalize(ca_next - protein[i].ca);

        protein[i].n = protein[i].ca - v_prev * 1.45;
        protein[i].c = protein[i].ca + v_next * 1.52;

        Vec3 ortho = cross(v_prev, v_next);
        if(length(ortho) < 0.1) ortho = {0,1,0};
        ortho = normalize(ortho);
        protein[i].o = protein[i].c + ortho * 1.23;
    }

    // PDB出力
    std::ofstream out("Universe_v16_1_Collagen.pdb");
    out << "HEADER    PROJECT-137 V16.1 COLLAGEN WIRE\n";

    // HELIXレコードの出力 (ビューワーにワイヤーをリボンとして描かせる)
    int helixId = 1;
    for(int i=0; i<(int)protein.size(); i++) {
        if(protein[i].isHelix || protein[i].isWire) { // ワイヤーもHELIXとして登録
            int start = i;
            while(i < (int)protein.size() && (protein[i].isHelix || protein[i].isWire)) i++;
            int end = i - 1;
            if(end - start >= 4) {
                char buf[100];
                sprintf(buf, "HELIX  %3d %3d %3s A %4d  %3s A %4d  1\n", 
                        helixId, helixId, protein[start].resName.c_str(), protein[start].resSeq,
                        protein[end].resName.c_str(), protein[end].resSeq);
                out << buf;
                helixId++;
            }
        }
    }

    int atomId = 1;
    for (auto& r : protein) {
        char buf[100];
        sprintf(buf, "ATOM  %5d  N   %3s A%4d    %8.3f%8.3f%8.3f  1.00 20.00           N\n", atomId++, r.resName.c_str(), r.resSeq, r.n.x, r.n.y, r.n.z); out << buf;
        sprintf(buf, "ATOM  %5d  CA  %3s A%4d    %8.3f%8.3f%8.3f  1.00 20.00           C\n", atomId++, r.resName.c_str(), r.resSeq, r.ca.x, r.ca.y, r.ca.z); out << buf;
        sprintf(buf, "ATOM  %5d  C   %3s A%4d    %8.3f%8.3f%8.3f  1.00 20.00           C\n", atomId++, r.resName.c_str(), r.resSeq, r.c.x, r.c.y, r.c.z); out << buf;
        sprintf(buf, "ATOM  %5d  O   %3s A%4d    %8.3f%8.3f%8.3f  1.00 20.00           O\n", atomId++, r.resName.c_str(), r.resSeq, r.o.x, r.o.y, r.o.z); out << buf;
    }

    out << "END\n";
    std::cout << "[Complete] Generated Collagen Wire!" << std::endl;
    return 0;
}

