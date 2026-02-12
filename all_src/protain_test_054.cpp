#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <cmath>
#include <map>
#include <algorithm>
#include <random>

// --- Project-137: Universe OS v43.0 "Ribosome Emulator" ---
// 生命の「時間」を導入。アミノ酸を一斉に配置するのではなく、
// リボソームから1個ずつ射出（スポーン）させながら、リアルタイムで物理演算を行う。
// これにより「全部つながった長いバネ（鉄パイプ）」の形成を防ぎ、
// 局所的な折り畳み（ドメイン形成）を完了させてから次の部分を作る「コ・トランスレーショナル」を実現する。

const double CA_DIST = 3.8;   
const double MIN_DIST = 4.0;  
const double DT = 0.01;       
const double DAMPING = 0.85;  
const double MAX_VEL = 2.5;   
const double TEMPERATURE = 0.5; // 生命の熱揺らぎ

// リボソームの合成速度（物理ステップ数 / 1アミノ酸）
// ゆっくり作るほど、途中で綺麗に畳まる時間ができる。
const int STEPS_PER_AA = 200; 

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
    bool active; // ★リボソームから射出されたかどうかのフラグ
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
    std::cout << "--- Project-137 v43.0 [Ribosome Emulator] ---" << std::endl;
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
    // 初期化（最初は全員非アクティブ）
    for(int i=0; i<(int)aaSeq.length(); ++i) {
        AAProps props = getAAProps(aaSeq[i]);
        // 最初は全部原点に隠しておく
        p[i] = {i+1, aaSeq[i], props.name, props.hydro, props.p_alpha, 
                {0,0,0}, {0,0,0}, {0,0,0}, false};
    }

    // 乱数生成器（生命の熱ノイズ）
    std::mt19937 gen(137);
    std::uniform_real_distribution<> dis(-1.0, 1.0);

    int N = p.size();
    // 合成にかかる時間 ＋ 折り畳みの猶予時間
    int TOTAL_STEPS = (N * STEPS_PER_AA) + 5000; 
    
    int spawned_count = 0; // 現在リボソームから出ている数

    for (int step = 0; step < TOTAL_STEPS; step++) {
        
        // --- リボソーム機能：一定時間ごとにアミノ酸を射出 ---
        if (step % STEPS_PER_AA == 0 && spawned_count < N) {
            p[spawned_count].active = true;
            
            // 新しいアミノ酸は、直前のアミノ酸の少し先にスポーンさせる
            // (最初の1個目は原点)
            if (spawned_count > 0) {
                Vec3 prev_pos = p[spawned_count-1].pos;
                // ランダムな方向に少しずらして配置（衝突を避ける）
                p[spawned_count].pos = prev_pos + Vec3{dis(gen)*0.1 + 3.8, dis(gen)*0.1, dis(gen)*0.1}; 
            } else {
                p[spawned_count].pos = {0,0,0};
            }
            spawned_count++;
        }

        // --- 物理演算ループ（アクティブな粒子のみ計算） ---
        for (int i = 0; i < spawned_count; i++) p[i].force = {0,0,0};

        // 1. 【熱力学】温度（ブラウン運動）
        for (int i = 0; i < spawned_count; i++) {
            Vec3 noise = {dis(gen), dis(gen), dis(gen)};
            p[i].force += noise * TEMPERATURE;
        }

        // 2. 【基本】結合バネ (3.8A)
        for (int i = 0; i < spawned_count - 1; i++) {
            Vec3 dir = p[i+1].pos - p[i].pos;
            double d = length(dir);
            double force_mag = 250.0 * (d - CA_DIST); 
            Vec3 f = normalize(dir) * force_mag;
            p[i].force += f;
            p[i+1].force -= f;
        }

        // 3. 【骨格】正しい定規（角度バネ i, i+2）
        for (int i = 0; i < spawned_count - 2; i++) {
            Vec3 dir = p[i+2].pos - p[i].pos;
            double d = length(dir);
            
            double target = 5.4; // Helix
            double strength = 50.0;

            bool is_collagen_region = (p[i+1].code == 'P' || p[i+1].code == 'G');
            if (is_collagen_region) {
                target = 6.0; // Collagen
            }

            double force_mag = strength * (d - target); 
            Vec3 f = normalize(dir) * force_mag;
            p[i].force += f;
            p[i+2].force -= f;
        }

        // 4. 【形状】カイラリティ強制 (スカラー三重積)
        // 4つ揃って初めて効力が発揮される -> リボソームから出た瞬間のねじれを決定する
        for (int i = 0; i < spawned_count - 3; i++) {
            double vol = scalar_triple(p[i].pos, p[i+1].pos, p[i+2].pos, p[i+3].pos);
            
            double target_vol = 12.0; // 右巻き
            double strength = 10.0;

            if (p[i].code=='P' || p[i+1].code=='P' || p[i+2].code=='P' || p[i+3].code=='P') {
                target_vol = -15.0; // 左巻き
                strength = 8.0;
            }

            // 末端保護: 端っこは強く固定
            if (i < 3 || i > N - 6) strength *= 2.0; 

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

        // 5. 【水素結合】Helix形成 (i, i+4)
        for (int i = 0; i < spawned_count - 4; i++) {
            if (p[i].code == 'P' || p[i+4].code == 'P' || p[i].code == 'G') continue; 
            
            Vec3 dir = p[i+4].pos - p[i].pos;
            double d = length(dir);
            
            double helix_strength = 5.0 * (p[i].p_alpha + p[i+4].p_alpha - 1.5); 
            if (helix_strength < 0) helix_strength = 0; // ヒンジ創発

            if (helix_strength > 0) {
                double force_mag = helix_strength * (d - 5.5); 
                Vec3 f = normalize(dir) * force_mag;
                p[i].force += f;
                p[i+4].force -= f;
            }
        }

        // 6. 【全体】疎水性引力・斥力・ワープ
        // アクティブな粒子同士でのみ相互作用する
        for (int i = 0; i < spawned_count; i++) {
            for (int j = i + 5; j < spawned_count; j++) {
                Vec3 dir = p[j].pos - p[i].pos;
                double d = length(dir);
                if (d < 0.1) d = 0.1;

                double force_mag = 0.0;
                
                // 斥力
                if (d < MIN_DIST) {
                    force_mag = -120.0 * std::pow((MIN_DIST - d), 2); 
                } 
                // 疎水引力
                else if (d < 30.0 && p[i].hydro > 0 && p[j].hydro > 0) {
                    force_mag = 15.0 * (p[i].hydro + p[j].hydro) / d; 
                }

                // シャペロンワープ
                if (p[i].code == 'C' && p[j].code == 'C' && d > 4.5) {
                    force_mag = 40.0; 
                }

                Vec3 f = normalize(dir) * force_mag;
                p[i].force += f;
                p[j].force -= f;
            }
        }

        // 7. ニュートンの運動方程式
        for (int i = 0; i < spawned_count; i++) {
            p[i].vel += p[i].force * DT;
            double v_len = length(p[i].vel);
            if (v_len > MAX_VEL) p[i].vel = normalize(p[i].vel) * MAX_VEL;
            p[i].vel = p[i].vel * DAMPING; 
            p[i].pos += p[i].vel * DT;
        }
    }

    std::ofstream out("Universe_v43_Ribosome.pdb");
    out << "HEADER    PROJECT-137 V43.0 RIBOSOME EMULATOR\n";
    int atomId = 1;
    for (int i=0; i<N; i++) {
        char buf[100];
        sprintf(buf, "ATOM  %5d  CA  %3s A%4d    %8.3f%8.3f%8.3f  1.00 20.00           C\n", 
                atomId++, p[i].name.c_str(), p[i].id, p[i].pos.x, p[i].pos.y, p[i].pos.z); 
        out << buf;
    }
    out << "END\n";
    std::cout << "[Complete] Check 'Universe_v43_Ribosome.pdb'!" << std::endl;
    return 0;
}
