#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <cmath>
#include <iomanip>
#include <map>
#include <algorithm>

// --- Project-137: Universe OS v29.0 "Water Shield & 137 Resonance" ---
// コンセプト: 
// 1. 1Dの魔法の数字（j - i > 45）を廃止。オキシトシンの「数珠化バグ」を解消。
// 2. 初期の空間距離（3D）に水分子が埋まることによる「通信遮断（Water Shield）」を実装。
//    これにより、カルモジュリンは2つのダンベルに分かれ、オキシトシンは1つの塊に美しく丸まる。
// 3. 螺旋の割り込み（P, G）と、距離ゼロワープ（SF公式）の完全統合。

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

// 宇宙OS「ハッシュ・パラメータ」
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
    std::cout << "--- Project-137 v29.0 [Water Shield & 137 Resonance] ---" << std::endl;
    std::cout << "Input Sequence: ";
    std::cin >> rawInput;
    std::transform(rawInput.begin(), rawInput.end(), rawInput.begin(), ::toupper);

    bool isNucleotide = true;
    for (char c : rawInput) {
        if (c != 'A' && c != 'T' && c != 'G' && c != 'C' && c != 'U') { isNucleotide = false; break; }
    }

    std::string aaSeq = "";
    if (isNucleotide && rawInput.length() >= 3) {
        for(char& c : rawInput) { if(c == 'T') c = 'U'; }
        for(size_t i = 0; i + 2 < rawInput.length(); i += 3) {
            std::string codon = rawInput.substr(i, 3);
            char aa = codonToAA(codon);
            if(aa == '_' || aa == '?') break;
            aaSeq += aa;
        }
    } else {
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

    // --- 1. OSプレレンダリング (デジタル・フラグによる判定) ---
    // P(曲げ)やG(自由)によるハードウェア割り込み
    int continuous_helix = 0;
    for(int i=0; i<(int)protein.size(); ++i) {
        if (protein[i].code == 'P') {
            protein[i].isWire = true; 
            protein[i].isHelix = false;
            continuous_helix = 0;
        } else if (protein[i].code == 'G') {
            protein[i].isWire = false;
            protein[i].isHelix = false; 
            continuous_helix = 0;
        } else if (protein[i].p_alpha >= 1.05 || protein[i].hydrophobicity > 1.5) {
            protein[i].isHelix = true;
            protein[i].isWire = false;
            continuous_helix++;
        } else {
            protein[i].isHelix = false;
            protein[i].isWire = false;
            continuous_helix = 0;
        }
    }

    Vec3 pos = {0, 0, 0};
    Vec3 axis = {0, 0, 1};
    Vec3 u = {1, 0, 0};
    Vec3 v = {0, 1, 0};
    double theta = 0;

    // 初期骨格の配置
    for(int i=0; i<(int)protein.size(); ++i) {
        if (protein[i].isWire) {
            pos = pos + axis * 2.9;
            theta += 2.094;
            protein[i].ca = pos + (u * std::cos(theta) + v * std::sin(theta)) * 1.4;
        } else if (protein[i].isHelix) {
            pos = pos + axis * 1.5;
            theta += 1.745;
            protein[i].ca = pos + (u * std::cos(theta) + v * std::sin(theta)) * 2.3;
        } else {
            double pseudo1 = std::abs(std::sin(i * 12.9898) * 43758.5453);
            double pseudo2 = std::abs(std::cos(i * 4.1414) * 43758.5453);
            double rx = std::sqrt(std::max(0.0, 1.0 - std::pow(2.0*pseudo2-1.0, 2))) * std::cos(2.0*PI*(pseudo1-std::floor(pseudo1)));
            double ry = std::sqrt(std::max(0.0, 1.0 - std::pow(2.0*pseudo2-1.0, 2))) * std::sin(2.0*PI*(pseudo1-std::floor(pseudo1)));
            double rz = 2.0*(pseudo2-std::floor(pseudo2)) - 1.0;

            axis.x += rx * 1.5;
            axis.y += ry * 1.5;
            axis.z += rz * 1.5;
            axis = normalize(axis);

            u = cross(axis, {0,1,0});
            if (length(u) < 0.1) u = cross(axis, {1,0,0});
            u = normalize(u); v = normalize(cross(axis, u));

            pos = pos + axis * 3.8;
            protein[i].ca = pos;
            theta = 0;
        }
    }

    // 【パッチ1: 水分子シールド（空間通信限界）の計測】
    // 配列番号の差 (j - i > 45) という人間的なインチキを廃止。
    // 初期の空間距離を記憶し、間に水が詰まってパケットが届かない限界値（38.0Å）をOSの仕様とする。
    std::vector<Vec3> initial_pos(protein.size());
    for(size_t i = 0; i < protein.size(); ++i) {
        initial_pos[i] = protein[i].ca;
    }

    // --- 2. 究極のリラクゼーション ---
    int total_iters = 600;
    for(int iter=0; iter<total_iters; ++iter) { 
        double anneal_factor = (double)iter / (double)total_iters;
        
        // 【螺旋のテンション(張力)デーモン】
        for(int i=0; i<(int)protein.size() - 3; ++i) {
            if (protein[i].isHelix && protein[i+1].isHelix && protein[i+2].isHelix && protein[i+3].isHelix) {
                Vec3 d3 = protein[i+3].ca - protein[i].ca;
                double dist3 = length(d3);
                if(dist3 > 0.1) {
                    double force3 = (dist3 - 5.4) * 0.05; // 5.4Aを保とうとする弱いバネ
                    Vec3 move3 = normalize(d3) * force3;
                    protein[i].ca = protein[i].ca + move3;
                    protein[i+3].ca = protein[i+3].ca - move3;
                }
            }
        }

        // 引力・斥力の計算
        for(int i=0; i<(int)protein.size(); ++i) {
            for(int j=i+1; j<(int)protein.size(); ++j) {
                
                // 【水和殻パッチ】初期空間距離によるOSのサボり判定
                // 初期の時点で 38.0A 以上離れている場合、間に水分子がビッシリ詰まって通信パケット（引力）が遮断される。
                // これにより、カルモジュリンは2つに分かれ、オキシトシンは全体が一気に丸まる。
                double init_dist = length(initial_pos[j] - initial_pos[i]);
                if (init_dist > 38.0) continue;

                double massI = 1.0;
                double massJ = 1.0;
                double totalMass = 2.0;

                Vec3 d = protein[j].ca - protein[i].ca;
                double euclid_dist = length(d); 
                double n_taxi = std::abs(d.x) + std::abs(d.y) + std::abs(d.z); // タクシー幾何学
                
                if (euclid_dist < 0.01) { d = {0.01, 0.01, 0.0}; euclid_dist = 0.01; n_taxi = 0.02; }
                
                // 【Soft Lock】めり込みを多少許容し、2Dペシャンコ化を防ぐ
                if(euclid_dist < MIN_DIST) {
                    double push = (MIN_DIST - euclid_dist) * 0.5;
                    if (euclid_dist < 3.0) push += (3.0 / euclid_dist); 
                    if (push > 5.0) push = 5.0; // マイルドなリミッター

                    Vec3 move = normalize(d) * push;
                    protein[i].ca = protein[i].ca - move * (massI / totalMass);
                    protein[j].ca = protein[j].ca + move * (massJ / totalMass);
                } 
                // 【引力 (SF オーバーライド)】
                else if (n_taxi < 15.0) { 
                    if (protein[i].hydrophobicity > 0.0 && protein[j].hydrophobicity > 0.0) {
                        double D = 2.0; 
                        
                        bool isCoreI = (protein[i].code == 'L' || protein[i].code == 'I' || protein[i].code == 'V');
                        bool isCoreJ = (protein[j].code == 'L' || protein[j].code == 'I' || protein[j].code == 'V');
                        
                        if (isCoreI && isCoreJ) D = 3.0; // ワープ引力

                        // 1.37 共鳴を組み込んだSF公式
                        double OS_CONST = 1.37;
                        double SF = 1.0 / std::pow((n_taxi / OS_CONST), D);
                        if (SF > 1.5) SF = 1.5; 
                        
                        double pull = SF * 3.0 * (protein[i].hydrophobicity + protein[j].hydrophobicity) * anneal_factor; 
                        if (pull > 0.6) pull = 0.6; // 吸着の抑制

                        Vec3 move = normalize(d) * pull;
                        protein[i].ca = protein[i].ca + move * (massI / totalMass);
                        protein[j].ca = protein[j].ca - move * (massJ / totalMass);
                    }
                }
            }
        }
        // ペプチド結合（鎖）の絶対維持
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

    // --- 3. フルアトム生成 ---
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

    std::ofstream out("Universe_v29_WaterShield.pdb");
    out << "HEADER    PROJECT-137 V29.0 WATER SHIELD & 137 RESONANCE\n";

    int helixId = 1;
    for(int i=0; i<(int)protein.size(); i++) {
        if(protein[i].isHelix || protein[i].isWire) {
            int start = i;
            while(i < (int)protein.size() && (protein[i].isHelix || protein[i].isWire)) i++;
            int end = i - 1;
            if(end - start >= 3) {
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
    std::cout << "[Complete] Check 'Universe_v29_WaterShield.pdb'!" << std::endl;
    return 0;
}

