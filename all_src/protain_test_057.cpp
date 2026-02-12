#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <cmath>
#include <map>
#include <algorithm>
#include <random>

// --- Project-137 v45.5 "Hyper-Viscosity & Inertia Control" ---
// 「竿の先に重りがついて変な方向に曲がる」という慣性の問題を解決するため、
// 細胞質のような「高粘性（ドロドロ）」の物理環境を実装。
// 1. 【粘性抵抗 (Viscous Drag)】: 速度に比例した強力なブレーキを常にかける。
// 2. 【慣性制御 (Inertia Damping)】: 鎖が長くなるほど（spawned_countが増えるほど）、
//    全体の抵抗値を上げて「暴れ」を物理的に封じ込める。
// これにより、カルモジュリンのリンカーや、プレプロオキシトシンの長い鎖も、
// 重力や慣性に負けず、リボソームの出口で「ピタッ」と静止したまま折り畳まれる。

const double CA_DIST = 3.8;   
const double MIN_DIST = 4.0;  
const double DT = 0.01;       
const double MAX_VEL = 2.0;   // 暴走を防ぐ上限速度
const double TEMPERATURE = 0.5;
const int STEPS_PER_AA = 200; 
const int TUNNEL_LENGTH = 45; 

// ★細胞質の粘り気（水の数十倍）★
const double BASE_VISCOSITY = 5.0; 

struct Vec3 {
    double x, y, z;
    Vec3 operator-(const Vec3& o) const { return {x - o.x, y - o.y, z - o.z}; }
    Vec3 operator+(const Vec3& o) const { return {x + o.x, y + o.y, z + o.z}; }
    Vec3 operator*(double s) const { return {x * s, y * s, z * s}; }
    Vec3& operator+=(const Vec3& o) { x += o.x; y += o.y; z += o.z; return *this; }
    Vec3& operator-=(const Vec3& o) { x -= o.x; y -= o.y; z -= o.z; return *this; }
};

Vec3 cross(Vec3 a, Vec3 b) { return {a.y*b.z - a.z*b.y, a.z*b.x - a.x*b.z, a.x*b.y - a.y*b.x}; }
double dot(Vec3 a, Vec3 b) { return a.x*b.x + a.y*b.y + a.z*b.z; }
double length(Vec3 v) { return std::sqrt(v.x*v.x + v.y*v.y + v.z*v.z); } 
Vec3 normalize(Vec3 v) { double l = length(v); return (l>1e-6)? v * (1.0/l) : Vec3{0,0,0}; }

double scalar_triple(Vec3 p0, Vec3 p1, Vec3 p2, Vec3 p3) {
    return dot(p1 - p0, cross(p2 - p1, p3 - p2));
}

struct AAProps { std::string name; double hydro; double p_alpha; double charge; };
AAProps getAAProps(char code) {
    static std::map<char, AAProps> m = {
        {'E',{"GLU",-3.5, 1.51, -1}}, {'M',{"MET", 1.9, 1.45, 0}}, {'A',{"ALA", 1.8, 1.42, 0}},
        {'L',{"LEU", 3.8, 1.21, 0}}, {'K',{"LYS",-3.9, 1.16, 1}}, {'F',{"PHE", 2.8, 1.13, 0}},
        {'Q',{"GLN",-3.5, 1.11, 0}}, {'W',{"TRP",-0.9, 1.08, 0}}, {'I',{"ILE", 4.5, 1.08, 0}},
        {'V',{"VAL", 4.2, 1.06, 0}}, {'D',{"ASP",-3.5, 1.01, -1}}, {'H',{"HIS",-3.2, 1.00, 0.5}}, 
        {'R',{"ARG",-4.5, 0.98, 1}}, {'T',{"THR",-0.7, 0.83, 0}}, {'S',{"SER",-0.8, 0.77, 0}},
        {'C',{"CYS", 2.5, 0.70, 0}}, {'Y',{"TYR",-1.3, 0.69, 0}}, {'N',{"ASN",-3.5, 0.67, 0}},
        {'P',{"PRO",-1.6, 0.57, 0}}, {'G',{"GLY",-0.4, 0.57, 0}}
    };
    if(m.count(code)) return m[code];
    return {"UNK",0, 0.5, 0};
}

struct Particle {
    int id; char code; std::string name;
    double hydro, p_alpha, charge;
    Vec3 pos, vel, force;
    bool active; 
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
    std::cout << "--- Project-137 v45.5 [Hyper-Viscosity] ---" << std::endl;
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
    for(int i=0; i<(int)aaSeq.length(); ++i) {
        AAProps props = getAAProps(aaSeq[i]);
        p[i] = {i+1, aaSeq[i], props.name, props.hydro, props.p_alpha, props.charge,
                {0,0,0}, {0,0,0}, {0,0,0}, false};
    }

    std::mt19937 gen(137);
    std::uniform_real_distribution<> dis(-1.0, 1.0);

    int N = p.size();
    int TOTAL_STEPS = (N * STEPS_PER_AA) + 15000; 
    
    int spawned_count = 0; 

    for (int step = 0; step < TOTAL_STEPS; step++) {
        
        // --- リボソーム射出 ---
        if (step % STEPS_PER_AA == 0 && spawned_count < N) {
            p[spawned_count].active = true;
            if (spawned_count > 0) {
                Vec3 prev_pos = p[spawned_count-1].pos;
                // まっすぐ静かに送り出す（暴れさせない）
                p[spawned_count].pos = prev_pos + Vec3{3.5, 0, 0}; 
            } else {
                p[spawned_count].pos = {0,0,0};
            }
            spawned_count++;
        }

        // --- 物理演算 ---
        for (int i = 0; i < spawned_count; i++) p[i].force = {0,0,0};

        // 1. 【熱力学】
        for (int i = 0; i < spawned_count; i++) {
            Vec3 noise = {dis(gen), dis(gen), dis(gen)};
            p[i].force += noise * TEMPERATURE;
        }

        // 2. 【結合バネ】(3.8A)
        for (int i = 0; i < spawned_count - 1; i++) {
            Vec3 dir = p[i+1].pos - p[i].pos;
            double d = length(dir);
            double force_mag = 300.0 * (d - CA_DIST); 
            Vec3 f = normalize(dir) * force_mag;
            p[i].force += f;
            p[i+1].force -= f;
        }

        // 3. 【骨格定規】(i, i+2)
        for (int i = 0; i < spawned_count - 2; i++) {
            Vec3 dir = p[i+2].pos - p[i].pos;
            double d = length(dir);
            double target = 5.4; 
            if (p[i+1].code == 'P' || p[i+1].code == 'G') target = 6.0; 
            double force_mag = 60.0 * (d - target); 
            Vec3 f = normalize(dir) * force_mag;
            p[i].force += f;
            p[i+2].force -= f;
        }

        // 4. 【カイラリティ強制】
        for (int i = 0; i < spawned_count - 3; i++) {
            double vol = scalar_triple(p[i].pos, p[i+1].pos, p[i+2].pos, p[i+3].pos);
            double target_vol = 12.0; 
            double strength = 12.0;
            if (p[i].code=='P' || p[i+1].code=='P' || p[i+2].code=='P' || p[i+3].code=='P') {
                target_vol = -15.0; 
                strength = 10.0;
            }
            if (strength > 0) {
                double diff = target_vol - vol;
                Vec3 v1 = p[i+1].pos - p[i].pos;
                Vec3 v2 = p[i+2].pos - p[i+1].pos;
                Vec3 normal = normalize(cross(v1, v2));
                Vec3 force = normal * (diff * strength * 0.01);
                p[i+1].force += force;
                p[i+2].force -= force; 
            }
        }

        // 5. 【水素結合】(i, i+4)
        for (int i = 0; i < spawned_count - 4; i++) {
            if (p[i].code == 'P' || p[i+4].code == 'P' || p[i].code == 'G') continue; 
            Vec3 dir = p[i+4].pos - p[i].pos;
            double d = length(dir);
            double helix_strength = 6.0 * (p[i].p_alpha + p[i+4].p_alpha - 1.5); 
            if (helix_strength < 0) helix_strength = 0; 
            if (helix_strength > 0) {
                double force_mag = helix_strength * (d - 5.5); 
                Vec3 f = normalize(dir) * force_mag;
                p[i].force += f;
                p[i+4].force -= f;
            }
        }

        // 6. 【相互作用の地平線 (Interaction Horizon)】
        int tunnel_exit_index = spawned_count - TUNNEL_LENGTH;
        if (tunnel_exit_index < 0) tunnel_exit_index = 0;

        for (int i = 0; i < tunnel_exit_index; i++) {
            for (int j = i + 5; j < tunnel_exit_index; j++) {
                Vec3 dir = p[j].pos - p[i].pos;
                double d = length(dir);
                if (d < 0.1) d = 0.1;

                double force_mag = 0.0;
                
                // (A) 斥力
                if (d < MIN_DIST) {
                    force_mag = -150.0 * std::pow((MIN_DIST - d), 2); 
                } 
                // (B) 疎水引力（近距離優先）
                else if (d < 30.0 && p[i].hydro > 0 && p[j].hydro > 0) {
                    int seq_dist = std::abs(i - j);
                    double local_bias = 1.0;
                    if (seq_dist < 25) local_bias = 3.0; else local_bias = 0.5;
                    force_mag = 20.0 * (p[i].hydro + p[j].hydro) * local_bias / d; 
                }
                // (C) 静電斥力
                if (d < 25.0 && std::abs(p[i].charge) > 0.1 && std::abs(p[j].charge) > 0.1) {
                    double electric = (p[i].charge * p[j].charge) * 25.0 / (d * d); 
                    force_mag -= electric; 
                }
                // (D) シャペロンワープ
                if (p[i].code == 'C' && p[j].code == 'C' && d > 4.5) {
                    force_mag = 50.0; 
                }

                Vec3 f = normalize(dir) * force_mag;
                p[i].force += f;
                p[j].force -= f;
            }
        }

        // 7. 【ニュートンの運動方程式（ハイパー・ダンピング）】
        // ★ここが今回の核心★
        // 鎖が長くなるほど（spawned_count増）、全体の粘性を高めて暴れを防ぐ。
        double current_viscosity = BASE_VISCOSITY * (1.0 + (double)spawned_count / 50.0);

        for (int i = 0; i < spawned_count; i++) {
            // 粘性抵抗力（速度に比例して逆向き）
            Vec3 drag_force = p[i].vel * (-current_viscosity);
            p[i].force += drag_force;

            p[i].vel += p[i].force * DT;
            
            // 速度リミッター
            double v_len = length(p[i].vel);
            if (v_len > MAX_VEL) p[i].vel = normalize(p[i].vel) * MAX_VEL;

            // 位置更新
            p[i].pos += p[i].vel * DT;
        }
    }

    std::ofstream out("Universe_v45_5_HyperViscosity.pdb");
    out << "HEADER    PROJECT-137 V45.5 HYPER-VISCOSITY ENGINE\n";
    int atomId = 1;
    for (int i=0; i<N; i++) {
        char buf[100];
        sprintf(buf, "ATOM  %5d  CA  %3s A%4d    %8.3f%8.3f%8.3f  1.00 20.00           C\n", 
                atomId++, p[i].name.c_str(), p[i].id, p[i].pos.x, p[i].pos.y, p[i].pos.z); 
        out << buf;
    }
    out << "END\n";
    std::cout << "[Complete] Check 'Universe_v45_5_HyperViscosity.pdb'!" << std::endl;
    return 0;
}

