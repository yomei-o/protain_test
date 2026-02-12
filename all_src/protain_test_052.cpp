#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <cmath>
#include <map>
#include <algorithm>

// --- Project-137: Universe OS v41.0 "Chirality Force (カイラリティ強制)" ---
// 座標計算のバグによる「螺旋の崩壊・逆転」を完全に防ぐため、
// 「スカラー三重積（カイラル体積）」を導入し、右巻き・左巻きを物理的に強制する。
// 1. 結合バネ (3.8A)
// 2. 角度バネ (1-3距離): Helix=5.4A, Collagen=6.0A
// 3. カイラリティ力 (スカラー三重積): Helix=+Target, Collagen=-Target
// 4. 水素結合 (1-4距離): Helix=5.5A

const double CA_DIST = 3.8;   
const double MIN_DIST = 4.0;  
const double DT = 0.01;       
const double DAMPING = 0.80;  // 振動を抑えるため減衰を強めに
const double MAX_VEL = 2.0;   

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

// ★スカラー三重積の計算（カイラリティ判定）★
// (p1-p0) . ((p2-p1) x (p3-p2))
double scalar_triple(Vec3 p0, Vec3 p1, Vec3 p2, Vec3 p3) {
    Vec3 v1 = p1 - p0;
    Vec3 v2 = p2 - p1;
    Vec3 v3 = p3 - p2;
    return dot(v1, cross(v2, v3));
}

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
    std::cout << "--- Project-137 v41.0 [Chirality Force Engine] ---" << std::endl;
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
    
    // 【初期配置】緩やかなジグザグで配置（一直線によるロックを防ぐ）
    for(int i=0; i<(int)aaSeq.length(); ++i) {
        AAProps props = getAAProps(aaSeq[i]);
        p[i] = {i+1, aaSeq[i], props.name, props.hydro, props.p_alpha, 
                {i * 3.0, std::sin(i*1.0), std::cos(i*1.0)}, {0,0,0}, {0,0,0}};
    }

    int N = p.size();
    int TOTAL_STEPS = 25000; 

    for (int step = 0; step < TOTAL_STEPS; step++) {
        for (auto& particle : p) particle.force = {0,0,0};

        // 1. 【基本】結合バネ (3.8A)
        for (int i = 0; i < N - 1; i++) {
            Vec3 dir = p[i+1].pos - p[i].pos;
            double d = length(dir);
            double force_mag = 200.0 * (d - CA_DIST); 
            Vec3 f = normalize(dir) * force_mag;
            p[i].force += f;
            p[i+1].force -= f;
        }

        // 2. 【骨格】角度バネ (i と i+2)
        // 5.5A (約90度) 〜 6.5A (約120度) を維持
        for (int i = 0; i < N - 2; i++) {
            Vec3 dir = p[i+2].pos - p[i].pos;
            double d = length(dir);
            double target = 5.8;
            if (p[i+1].code == 'P' || p[i+1].code == 'G') target = 6.0; // コラーゲンは伸び気味
            
            double force_mag = 50.0 * (d - target); 
            Vec3 f = normalize(dir) * force_mag;
            p[i].force += f;
            p[i+2].force -= f;
        }

        // 3. 【究極パッチ：カイラリティ強制 (スカラー三重積)】
        // 4点 (i, i+1, i+2, i+3) の体積を計算し、正しい符号（巻き方向）へ誘導する
        for (int i = 0; i < N - 3; i++) {
            // 現在のカイラリティ（体積）
            double vol = scalar_triple(p[i].pos, p[i+1].pos, p[i+2].pos, p[i+3].pos);
            
            double target_vol = 0.0;
            double strength = 0.0;

            bool is_collagen = (p[i].code=='P' || p[i+1].code=='P' || p[i+2].code=='P' || p[i+3].code=='P');

            if (is_collagen) {
                // 【左巻き強制】コラーゲンはマイナスの体積を持つべき
                target_vol = -15.0; // 理想的な左巻きの強さ
                strength = 5.0;
            } else if (p[i].code != 'G' && p[i+1].code != 'G') { // Gは自由
                // 【右巻き強制】アルファヘリックスはプラスの体積を持つべき
                target_vol = 12.0; // 理想的な右巻きの強さ
                strength = 8.0;
            }

            if (strength > 0) {
                // 簡易的な復元力：i+1 と i+2 をずらして体積を修正する（厳密な勾配ではないが近似的に効く）
                double diff = target_vol - vol;
                Vec3 v1 = p[i+1].pos - p[i].pos;
                Vec3 v2 = p[i+2].pos - p[i+1].pos;
                Vec3 normal = normalize(cross(v1, v2)); // ねじれ軸方向
                
                // 軸方向に力を加えてねじれを修正
                Vec3 force = normal * (diff * strength * 0.01);
                p[i+1].force += force;
                p[i+2].force -= force; 
            }
        }

        // 4. 【形状】水素結合バネ (i と i+4)
        for (int i = 0; i < N - 4; i++) {
            if (p[i].code == 'P' || p[i+4].code == 'P' || p[i].code == 'G') continue; 
            
            Vec3 dir = p[i+4].pos - p[i].pos;
            double d = length(dir);
            
            // 関節の創発（p_alphaチェック）
            double helix_strength = 5.0 * (p[i].p_alpha + p[i+4].p_alpha - 1.5); 
            if (helix_strength > 0) {
                double force_mag = helix_strength * (d - 5.5); // 5.5Aに近づける
                Vec3 f = normalize(dir) * force_mag;
                p[i].force += f;
                p[i+4].force -= f;
            }
        }

        // 5. 【全体】疎水性引力・斥力・ワープ
        for (int i = 0; i < N; i++) {
            for (int j = i + 5; j < N; j++) {
                Vec3 dir = p[j].pos - p[i].pos;
                double d = length(dir);
                if (d < 0.1) d = 0.1;

                double force_mag = 0.0;
                
                // 斥力
                if (d < MIN_DIST) {
                    force_mag = -100.0 * std::pow((MIN_DIST - d), 2); 
                } 
                // 疎水引力
                else if (d < 25.0 && p[i].hydro > 0 && p[j].hydro > 0) {
                    force_mag = 10.0 * (p[i].hydro + p[j].hydro) / d; 
                }

                // システインワープ
                if (p[i].code == 'C' && p[j].code == 'C' && d > 4.5) {
                    force_mag = 30.0; 
                }

                Vec3 f = normalize(dir) * force_mag;
                p[i].force += f;
                p[j].force -= f;
            }
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

    std::ofstream out("Universe_v41_Chirality.pdb");
    out << "HEADER    PROJECT-137 V41.0 CHIRALITY FORCE ENGINE\n";
    int atomId = 1;
    for (int i=0; i<N; i++) {
        char buf[100];
        sprintf(buf, "ATOM  %5d  CA  %3s A%4d    %8.3f%8.3f%8.3f  1.00 20.00           C\n", 
                atomId++, p[i].name.c_str(), p[i].id, p[i].pos.x, p[i].pos.y, p[i].pos.z); 
        out << buf;
    }
    out << "END\n";
    std::cout << "[Complete] Check 'Universe_v41_Chirality.pdb'!" << std::endl;
    return 0;
}

