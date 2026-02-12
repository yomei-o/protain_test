#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <cmath>
#include <map>
#include <algorithm>

// --- Project-137: Universe OS v40.0 "Elastic Coil Spring (弾性コイルばねモデル)" ---
// おじさんの発案により、距離のバネ（2Dフラット化の元凶）を廃止し、
// 3D空間の「ねじれ・曲がり」を直接計算する「コイルばねモデル」を導入！
// 3つのアミノ酸が作る局所座標系に対して「4つ目の理想的な3D位置」を計算し、そこへ誘導する。

const double CA_DIST = 3.8;   
const double MIN_DIST = 4.0;  
const double DT = 0.01;       
const double DAMPING = 0.85;  
const double MAX_VEL = 2.5;   

struct Vec3 {
    double x, y, z;
    Vec3 operator-(const Vec3& o) const { return {x - o.x, y - o.y, z - o.z}; }
    Vec3 operator+(const Vec3& o) const { return {x + o.x, y + o.y, z + o.z}; }
    Vec3 operator*(double s) const { return {x * s, y * s, z * s}; }
    Vec3& operator+=(const Vec3& o) { x += o.x; y += o.y; z += o.z; return *this; }
    Vec3& operator-=(const Vec3& o) { x -= o.x; y -= o.y; z -= o.z; return *this; }
};

Vec3 cross(Vec3 a, Vec3 b) { return {a.y*b.z - a.z*b.y, a.z*b.x - a.x*b.z, a.x*b.y - a.y*b.x}; }
double length(Vec3 v) { return std::sqrt(v.x*v.x + v.y*v.y + v.z*v.z); } 
Vec3 normalize(Vec3 v) { double l = length(v); return (l>1e-6)? v * (1.0/l) : Vec3{0,0,0}; }

struct AAProps { std::string name; double hydro; double p_alpha; };
AAProps getAAProps(char code) {
    static std::map<char, AAProps> m = {
        {'E',{"GLU",-3.5, 1.51}}, {'M',{"MET", 1.9, 1.45}}, {'A',{"ALA", 1.8, 1.42}},
        {'L',{"LEU", 3.8, 1.21}}, {'K',{"LYS",-3.9, 1.16}}, {'F',{"PHE", 2.8, 1.13}},
        {'Q',{"GLN",-3.5, 1.11}}, {'W',{"TRP",-0.9, 1.08}}, {'I',{"ILE", 4.5, 1.08}},
        {'V',{"VAL", 4.2, 1.06}}, {'D',{"ASP",-3.5, 1.01}}, {'H',{"HIS",-3.2, 1.00}},
        {'R',{"ARG",-4.5, 0.98}}, {'T',{"THR",-0.7, 0.83}}, {'S',{"SER",-0.8, 0.77}},
        {'C',{"CYS", 2.5, 0.70}}, {'Y',{"TYR",-1.3, 0.69}}, {'N',{"ASN",-3.5, 0.67}},
        {'P',{"PRO",-1.6, 0.57}}, {'G',{"GLY",-0.4, 0.57}}
    };
    if(m.count(code)) return m[code];
    return {"UNK",0, 0.5};
}

struct Particle {
    int id; char code; std::string name;
    double hydro, p_alpha;
    Vec3 pos, vel, force;
};

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
    std::cout << "--- Project-137 v40.0 [Elastic Coil Spring (3Dコイルばね)] ---" << std::endl;
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

    std::vector<Particle> p(aaSeq.length());
    
    // 【初期配置】対称性を破るノイズ
    for(int i=0; i<(int)aaSeq.length(); ++i) {
        AAProps props = getAAProps(aaSeq[i]);
        double noise_y = std::sin(i * 13.5) * 0.5; 
        double noise_z = std::cos(i * 7.2) * 0.5;  
        p[i] = {i+1, aaSeq[i], props.name, props.hydro, props.p_alpha, 
                {i * CA_DIST, noise_y, noise_z}, {0,0,0}, {0,0,0}};
    }

    int N = p.size();
    int TOTAL_STEPS = 20000; 

    for (int step = 0; step < TOTAL_STEPS; step++) {
        for (auto& particle : p) particle.force = {0,0,0};

        // 1. 【基本骨格】結合距離 (i と i+1 を 3.8A に保つ)
        for (int i = 0; i < N - 1; i++) {
            Vec3 dir = p[i+1].pos - p[i].pos;
            double d = length(dir);
            double force_mag = 150.0 * (d - CA_DIST); 
            Vec3 f = normalize(dir) * force_mag;
            p[i].force += f;
            p[i+1].force -= f;
        }

        // 2. 【究極パッチ：3D形状記憶コイルばね】(i, i+1, i+2 から i+3 を導く)
        for (int i = 0; i < N - 3; i++) {
            // ローカル座標系（直交基底）の作成
            Vec3 u = normalize(p[i+1].pos - p[i].pos);
            Vec3 v = normalize(p[i+2].pos - p[i+1].pos);
            Vec3 w = normalize(cross(u, v)); // 面の法線（ねじれ方向）
            Vec3 n = cross(w, v);            // 内側へ向かうベクトル

            if (length(w) < 0.1) continue; // 一直線になって外積が死んだらスキップ

            // プロリン(P)の含有量で「右巻きコイル」か「左巻きワイヤー」かを自動判定
            int p_count = (p[i].code=='P') + (p[i+1].code=='P') + (p[i+2].code=='P') + (p[i+3].code=='P');
            double target_a, target_b, target_c;

            if (p_count >= 1 || p[i].code == 'G') {
                // 【左巻き・コラーゲンモード (Polyproline II)】
                target_a = 1.39;  // 進行方向
                target_b = -0.98; // 外側
                target_c = -3.80; // 下へねじれる
            } else {
                // 【右巻き・アルファヘリックスモード】
                target_a = 0.02;  // 進行方向
                target_b = 2.46;  // 内側へ丸まる
                target_c = 2.94;  // 上へねじれる
            }

            // 関節の柔らかさ（p_alphaによる判定）
            double helix_strength = 15.0 * (p[i].p_alpha + p[i+3].p_alpha - 2.0);
            if (p_count >= 1 || p[i].code == 'G') helix_strength = 15.0; // コラーゲンは常に強固
            if (helix_strength < 0) helix_strength = 0; // p_alphaが低ければ「ヒンジ(関節)」になる

            if (helix_strength > 0) {
                // 理想の i+3 の位置を算出（絶対的な3D座標）
                Vec3 ideal_p3 = p[i+2].pos + (v * target_a) + (n * target_b) + (w * target_c);
                
                // 現在の i+3 の位置から、理想の位置へ強力なゴム紐で引っ張る！
                Vec3 correction = ideal_p3 - p[i+3].pos;
                p[i+3].force += correction * helix_strength;
                
                // 反作用（ニュートンの法則）で、他の3点も逆向きに調整
                p[i+2].force -= correction * (helix_strength * 0.33);
                p[i+1].force -= correction * (helix_strength * 0.33);
                p[i].force   -= correction * (helix_strength * 0.33);
            }
        }

        // 3. 【力D】疎水性引力、斥力、シャペロンワープ
        for (int i = 0; i < N; i++) {
            for (int j = i + 5; j < N; j++) {
                Vec3 dir = p[j].pos - p[i].pos;
                double d = length(dir);
                if (d < 0.1) d = 0.1;

                double force_mag = 0.0;
                
                // 斥力（4.0A以内のめり込み防止）
                if (d < MIN_DIST) {
                    force_mag = -100.0 * std::pow((MIN_DIST - d), 2); 
                } 
                // 疎水引力（EFハンドなどを曲げる力）
                else if (d < 25.0 && p[i].hydro > 0 && p[j].hydro > 0) {
                    force_mag = 12.0 * (p[i].hydro + p[j].hydro) / d; 
                }

                // シャペロンワープ（オキシトシンの結び目を作る絶対引力）
                if (p[i].code == 'C' && p[j].code == 'C' && d > 4.5) {
                    force_mag = 30.0; 
                }

                Vec3 f = normalize(dir) * force_mag;
                p[i].force += f;
                p[j].force -= f;
            }
        }

        // 4. ニュートンの運動方程式
        for (int i = 0; i < N; i++) {
            p[i].vel += p[i].force * DT;
            
            double v_len = length(p[i].vel);
            if (v_len > MAX_VEL) p[i].vel = normalize(p[i].vel) * MAX_VEL;

            p[i].vel = p[i].vel * DAMPING; 
            p[i].pos += p[i].vel * DT;
        }
    }

    std::ofstream out("Universe_v40_CoilSpring.pdb");
    out << "HEADER    PROJECT-137 V40.0 ELASTIC COIL SPRING ENGINE\n";
    int atomId = 1;
    for (int i=0; i<N; i++) {
        char buf[100];
        sprintf(buf, "ATOM  %5d  CA  %3s A%4d    %8.3f%8.3f%8.3f  1.00 20.00           C\n", 
                atomId++, p[i].name.c_str(), p[i].id, p[i].pos.x, p[i].pos.y, p[i].pos.z); 
        out << buf;
    }
    out << "END\n";
    std::cout << "[Complete] Check 'Universe_v40_CoilSpring.pdb'!" << std::endl;
    return 0;
}
