#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <cmath>
#include <map>
#include <algorithm>

// --- Project-137: Universe OS v39.0 "Taxi-Geometry & 137 Grid Snapping" ---
// 連続空間による「知恵の輪（局所最適解）」のバグを突破するため、案2を実装。
// 1. 空間を滑らかな連続ではなく「OS_GRID = 1.37 A」のマス目（デジタル空間）と定義する。
// 2. 常に最も近い1.37の倍数の座標へスナップ（量子化）しようとする力を加える。
// 3. 遠距離の引力（ワープ）は、ユークリッド距離ではなく「タクシー幾何学距離」で判定する。

const double CA_DIST = 3.8;   
const double MIN_DIST = 4.0;  
const double DT = 0.01;       
const double DAMPING = 0.85;  
const double MAX_VEL = 2.5;   

// ★宇宙OSのデジタルグリッド定数★
const double OS_GRID = 1.37;

struct Vec3 {
    double x, y, z;
    Vec3 operator-(const Vec3& o) const { return {x - o.x, y - o.y, z - o.z}; }
    Vec3 operator+(const Vec3& o) const { return {x + o.x, y + o.y, z + o.z}; }
    Vec3 operator*(double s) const { return {x * s, y * s, z * s}; }
    Vec3& operator+=(const Vec3& o) { x += o.x; y += o.y; z += o.z; return *this; }
    Vec3& operator-=(const Vec3& o) { x -= o.x; y -= o.y; z -= o.z; return *this; }
};

double length(Vec3 v) { return std::sqrt(v.x*v.x + v.y*v.y + v.z*v.z); } 
Vec3 normalize(Vec3 v) { double l = length(v); return (l>1e-6)? v * (1.0/l) : Vec3{0,0,0}; }

// ★タクシー幾何学による距離計算★
double taxi_length(Vec3 v) { return std::abs(v.x) + std::abs(v.y) + std::abs(v.z); }

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
    std::cout << "--- Project-137 v39.0 [137 Grid & Taxi Geometry] ---" << std::endl;
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
    
    // 【初期配置】
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

        // 1. 【力A】結合のバネ (距離 3.8A) - ここはユークリッド距離（鎖が切れないため）
        for (int i = 0; i < N - 1; i++) {
            Vec3 dir = p[i+1].pos - p[i].pos;
            double d = length(dir);
            double force_mag = 150.0 * (d - CA_DIST); 
            Vec3 f = normalize(dir) * force_mag;
            p[i].force += f;
            p[i+1].force -= f;
        }

        // 2. 【力B】第一関節バネ (i と i+2) 
        for (int i = 0; i < N - 2; i++) {
            Vec3 dir = p[i+2].pos - p[i].pos;
            double d = length(dir);
            double target_dist = 5.8; 
            if (p[i+1].code == 'P') target_dist = 5.5; 
            double force_mag = 50.0 * (d - target_dist); 
            Vec3 f = normalize(dir) * force_mag;
            p[i].force += f;
            p[i+2].force -= f;
        }

        // 3. 【力C】水素結合のバネ (i と i+4)
        for (int i = 0; i < N - 4; i++) {
            if (p[i].code == 'P' || p[i+4].code == 'P' || p[i].code == 'G') continue; 
            Vec3 dir = p[i+4].pos - p[i].pos;
            double d = length(dir);
            double helix_strength = 5.0 * (p[i].p_alpha + p[i+4].p_alpha - 1.5); 
            if (helix_strength > 0) {
                double force_mag = helix_strength * (d - 5.5); 
                Vec3 f = normalize(dir) * force_mag;
                p[i].force += f;
                p[i+4].force -= f;
            }
        }

        // 4. 【力D】疎水性引力、斥力、シャペロンワープ（タクシー幾何学による計算）
        for (int i = 0; i < N; i++) {
            for (int j = i + 5; j < N; j++) {
                Vec3 dir = p[j].pos - p[i].pos;
                
                // ★ユークリッド距離（直線距離）ではなく、タクシー距離（ブロックの歩数）で評価★
                double d_taxi = taxi_length(dir);
                double d_euclid = length(dir);
                if (d_euclid < 0.1) { d_euclid = 0.1; d_taxi = 0.1; }

                double force_mag = 0.0;
                
                // 斥力（4.0A以内のめり込み防止）
                if (d_euclid < MIN_DIST) {
                    force_mag = -100.0 * std::pow((MIN_DIST - d_euclid), 2); 
                } 
                // 疎水引力（タクシー距離で評価）
                else if (d_taxi < 35.0 && p[i].hydro > 0 && p[j].hydro > 0) {
                    force_mag = 15.0 * (p[i].hydro + p[j].hydro) / d_taxi; 
                }

                // シャペロンワープ（システイン同士）
                if (p[i].code == 'C' && p[j].code == 'C' && d_euclid > 4.5) {
                    force_mag = 30.0; 
                }

                Vec3 f = normalize(dir) * force_mag;
                p[i].force += f;
                p[j].force -= f;
            }
        }

        // 5. 【新パッチ：137グリッド・スナップ（量子化）】
        // 中途半端な位置（局所最適解）に留まることを許さず、
        // 常に「1.37 A」の倍数の座標にカチッ、カチッと吸着させる。
        for (int i = 0; i < N; i++) {
            Vec3 nearest_grid = {
                std::round(p[i].pos.x / OS_GRID) * OS_GRID,
                std::round(p[i].pos.y / OS_GRID) * OS_GRID,
                std::round(p[i].pos.z / OS_GRID) * OS_GRID
            };
            
            // グリッドの中心に向かって引っ張る力
            Vec3 grid_force = (nearest_grid - p[i].pos) * 5.0; 
            p[i].force += grid_force;
        }

        // 6. ニュートンの運動方程式
        for (int i = 0; i < N; i++) {
            p[i].vel += p[i].force * DT;
            
            double v_len = length(p[i].vel);
            if (v_len > MAX_VEL) p[i].vel = normalize(p[i].vel) * MAX_VEL;

            p[i].vel = p[i].vel * DAMPING; 
            p[i].pos += p[i].vel * DT;
        }
    }

    std::ofstream out("Universe_v39_GridSnap.pdb");
    out << "HEADER    PROJECT-137 V39.0 TAXI GEOMETRY 137 GRID SNAP\n";
    int atomId = 1;
    for (int i=0; i<N; i++) {
        char buf[100];
        sprintf(buf, "ATOM  %5d  CA  %3s A%4d    %8.3f%8.3f%8.3f  1.00 20.00           C\n", 
                atomId++, p[i].name.c_str(), p[i].id, p[i].pos.x, p[i].pos.y, p[i].pos.z); 
        out << buf;
    }
    out << "END\n";
    std::cout << "[Complete] Check 'Universe_v39_GridSnap.pdb'!" << std::endl;
    return 0;
}
