#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <cmath>
#include <iomanip>
#include <map>
#include <algorithm>

// --- Project-137: Universe OS v18.0 "Local Frame Engine" ---
// コンセプト: 絶対座標(angleX)を廃止。向いている方向(Axis)を基準に螺旋を作る。
// 疎水性引力（中心に集まる力）を実装し、細長いソーセージ化を防ぐ。

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
    bool isHelix; bool isWire;
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

char codonToAA(const std::string& codon) {
     static std::map<std::string, char> m = {
        {"UUU",'F'}, {"UUC",'F'}, {"UUA",'L'}, {"UUG",'L'}, {"CUU",'L'}, {"CUC",'L'}, {"CUA",'L'}, {"CUG",'L'},
        {"AUU",'I'}, {"AUC",'I'}, {"AUA",'I'}, {"AUG",'M'}, {"GUU",'V'}, {"GUC",'V'}, {"GUA",'V'}, {"GUG",'V'},
        {"UCU",'S'}, {"UCC",'S'}, {"UCA",'S'}, {"UCG",'S'}, {"CCU",'P'}, {"CCC",'P'}, {"CCA",'P'}, {"CCG",'P'},
        {"ACU",'T'}, {"ACC",'T'}, {"ACA",'T'}, {"ACG",'T'}, {"GCU",'A'}, {"GCC",'A'}, {"GCA",'A'}, {"GCG",'A'},
        {"UAU",'Y'}, {"UAC",'Y'}, {"UAA",'_'}, {"UAG",'_'}, {"CAU",'H'}, {"CAC",'H'}, {"CAA",'Q'}, {"CAG",'Q'},
        {"AAU",'N'}, {"AAC",'N'}, {"AAA",'K'}, {"AAG",'K'}, {"GAU",'D'}, {"GAC",'D'}, {"GAA",'E'}, {"GAG",'E'},
        {"UGU",'C'}, {"UGC",'C'}, {"UGA",'_'}, {"UGG",'W'}, {"CGU",'R'}, {"CGC",'R'}, {"CGA",'R'}, {"CGG",'R'},
        {"AGU",'S'}, {"AGC",'S'}, {"AGA",'R'}, {"AGG",'R'}, {"GGU",'G'}, {"GGC",'G'}, {"GGA",'G'}, {"GGG",'G'}
    };
    if (m.count(codon)) return m[codon];
    return '?';
}

int main() {
    std::string rawInput;
    std::cout << "--- Project-137 v18.0 [Local Frame Engine] ---" << std::endl;
    std::cout << "Input Sequence (RNA/DNA OR Amino Acid): ";
    std::cin >> rawInput;
    std::transform(rawInput.begin(), rawInput.end(), rawInput.begin(), ::toupper);

    bool isNucleotide = true;
    for (char c : rawInput) {
        if (c != 'A' && c != 'T' && c != 'G' && c != 'C' && c != 'U') {
            isNucleotide = false; break;
        }
    }

    std::string aaSeq = "";
    if (isNucleotide && rawInput.length() >= 3) {
        std::cout << "[System] Genetic code detected. Booting Ribosome..." << std::endl;
        for(char& c : rawInput) { if(c == 'T') c = 'U'; }
        for(size_t i = 0; i + 2 < rawInput.length(); i += 3) {
            std::string codon = rawInput.substr(i, 3);
            char aa = codonToAA(codon);
            if(aa == '_' || aa == '?') break;
            aaSeq += aa;
        }
    } else {
        std::cout << "[System] Amino Acid sequence detected. Bypassing Ribosome..." << std::endl;
        aaSeq = rawInput;
    }

    if (aaSeq.empty()) return 1;

    std::vector<Residue> protein;
    for(int i=0; i<(int)aaSeq.length(); ++i) {
        Residue res; res.resSeq = i+1; res.code = aaSeq[i];
        AAProps p = getAAProps(res.code);
        res.resName = p.name; res.hydrophobicity = p.hydro; 
        res.charge = p.charge; res.p_alpha = p.p_alpha;
        res.isHelix = false; res.isWire = false;
        protein.push_back(res);
    }

    // --- ローカル座標系（Frenet-Serret frame）による立体の構築 ---
    Vec3 pos = {0, 0, 0};       // 現在の中心位置
    Vec3 axis = {0, 0, 1};      // 進行方向のベクトル
    Vec3 u = {1, 0, 0};         // 横方向1
    Vec3 v = {0, 1, 0};         // 横方向2
    double theta = 0;           // 螺旋の回転角

    for(int i=0; i<(int)protein.size(); ++i) {
        double avgProp = 0; double proCount = 0; int count = 0;
        for(int k=i-2; k<=i+2; k++) {
            if(k>=0 && k<(int)protein.size()) { 
                avgProp += protein[k].p_alpha; 
                if (protein[k].code == 'P') proCount++; 
                count++; 
            }
        }
        avgProp /= (count>0?count:1);
        double proRatio = proCount / count;

        if (proRatio >= 0.4) {
            // 【コラーゲン】直進するワイヤー。Axisは変更しない！
            protein[i].isWire = true;
            pos = pos + axis * 2.9; // 1歩が長い
            theta += 2.094;         // 120度回転
            protein[i].ca = pos + (u * std::cos(theta) + v * std::sin(theta)) * 1.4; // 半径1.4
        } else if (avgProp > 1.05 && protein[i].code != 'P' && protein[i].code != 'G') {
            // 【螺旋】直進するバネ。Axisは変更しない！
            protein[i].isHelix = true;
            pos = pos + axis * 1.5; // 1歩が短い
            theta += 1.745;         // 100度回転
            protein[i].ca = pos + (u * std::cos(theta) + v * std::sin(theta)) * 2.3; // 半径2.3
        } else {
            // 【関節】ここでAxis（進行方向）を曲げる！
            
            // 疎水性引力：今までの重心（中心）を計算し、そこへ向かう力を作る
            Vec3 com = {0,0,0};
            for(int j=0; j<i; j++) com = com + protein[j].ca;
            if(i > 0) com = com * (1.0 / i);
            
            Vec3 toCenter = com - pos;
            if (length(toCenter) > 0) toCenter = normalize(toCenter);

            double hydroPull = (protein[i].hydrophobicity > 0) ? 1.5 : -0.5; // 油なら中心へ！
            double scatter = 1.0 + std::abs(protein[i].charge);              // 電荷でランダムに暴れる

            axis.x += (std::sin(i * 1.3) * scatter) + toCenter.x * hydroPull;
            axis.y += (std::cos(i * 1.7) * scatter) + toCenter.y * hydroPull;
            axis.z += (std::sin(i * 2.3) * scatter) + toCenter.z * hydroPull;
            axis = normalize(axis); // 新しい進行方向が決定！

            // 横方向ベクトル(u,v)を再計算
            u = cross(axis, {0,1,0});
            if (length(u) < 0.1) u = cross(axis, {1,0,0});
            u = normalize(u);
            v = normalize(cross(axis, u));

            pos = pos + axis * 3.8; // 関節は大きく進む
            protein[i].ca = pos;
            theta = 0; // 次の螺旋のためにリセット
        }
    }

    // リラクゼーション（剛体を守りつつ、関節の距離を整える）
    for(int iter=0; iter<300; ++iter) {
        for(int i=0; i<(int)protein.size(); ++i) {
            for(int j=i+1; j<(int)protein.size(); ++j) {
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
        for(int i=0; i<(int)protein.size()-1; ++i) {
             Vec3 d = protein[i+1].ca - protein[i].ca;
             double dist = length(d);
             if(dist > 0.1) {
                 double correction = (dist - CA_DIST) * 0.5;
                 Vec3 move = normalize(d) * correction;
                 protein[i].ca = protein[i].ca + move;
                 protein[i+1].ca = protein[i+1].ca - move;
             }
        }
    }

    // フルアトム計算
    for(int i=0; i<(int)protein.size(); ++i) {
        Vec3 ca_prev = (i > 0) ? protein[i-1].ca : protein[i].ca - Vec3{1,0,0};
        Vec3 ca_next = (i < (int)protein.size()-1) ? protein[i+1].ca : protein[i].ca + Vec3{1,0,0};
        Vec3 v_prev = normalize(protein[i].ca - ca_prev);
        Vec3 v_next = normalize(ca_next - protein[i].ca);

        protein[i].n = protein[i].ca - v_prev * 1.45;
        protein[i].c = protein[i].ca + v_next * 1.52;
        Vec3 ortho = cross(v_prev, v_next);
        if(length(ortho) < 0.1) ortho = {0,1,0};
        protein[i].o = protein[i].c + normalize(ortho) * 1.23;
    }

    std::ofstream out("Universe_v18_LocalFrame.pdb");
    out << "HEADER    PROJECT-137 V18.0 LOCAL FRAME ENGINE\n";

    int helixId = 1;
    for(int i=0; i<(int)protein.size(); i++) {
        if(protein[i].isHelix || protein[i].isWire) {
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
    std::cout << "[Complete] Check 'Universe_v18_LocalFrame.pdb'!" << std::endl;
    return 0;
}

