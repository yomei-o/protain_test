#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <cmath>
#include <map>
#include <algorithm>

// --- Project-137: Universe OS v32.0 "Emergence Physics Engine" ---
// if文による強制的なHelix設定などを完全廃止。
// 4つの「力（フォース）」の釣り合いだけで、自発的に螺旋やドメインを形成させる。

const double CA_DIST = 3.8;   // 隣り合うアミノ酸の距離
const double MIN_DIST = 4.0;  // 斥力（弾き返す力）が働く距離
const double DT = 0.01;       // 時間の刻み幅
const double DAMPING = 0.95;  // 水の摩擦（速度を殺して安定させる）

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

// コドン変換は省略（直接アミノ酸配列を入力）
int main() {
    std::string aaSeq;
    std::cout << "--- Project-137 v32.0 [Emergence Physics Engine] ---\nInput Amino Acid Sequence: ";
    std::cin >> aaSeq;
    std::transform(aaSeq.begin(), aaSeq.end(), aaSeq.begin(), ::toupper);

    std::vector<Particle> p(aaSeq.length());
    
    // 初期配置（紐を一直線に伸ばしておく）
    for(int i=0; i<(int)aaSeq.length(); ++i) {
        AAProps props = getAAProps(aaSeq[i]);
        p[i] = {i+1, aaSeq[i], props.name, props.hydro, props.p_alpha, 
                {i * CA_DIST, 0.0, 0.0}, {0,0,0}, {0,0,0}};
    }

    int N = p.size();
    int TOTAL_STEPS = 10000; // 時間をかけてゆっくり折り畳む

    for (int step = 0; step < TOTAL_STEPS; step++) {
        // 1. 力のリセット
        for (auto& particle : p) particle.force = {0,0,0};

        // 2. 結合のバネ (隣同士を 3.8A に保つ)
        for (int i = 0; i < N - 1; i++) {
            Vec3 dir = p[i+1].pos - p[i].pos;
            double d = length(dir);
            double force_mag = 50.0 * (d - CA_DIST); // フックの法則
            Vec3 f = normalize(dir) * force_mag;
            p[i].force += f;
            p[i+1].force -= f;
        }

        // 3. 水素結合のバネ (i と i+4 が近づこうとして、勝手に螺旋になる)
        for (int i = 0; i < N - 4; i++) {
            // プロリン(P)やグリシン(G)がいれば螺旋の力を弱める
            if (p[i].code == 'P' || p[i+4].code == 'P' || p[i].code == 'G') continue;
            
            Vec3 dir = p[i+4].pos - p[i].pos;
            double d = length(dir);
            // 螺旋になりやすいアミノ酸ほど引力が強い
            double helix_strength = 2.0 * (p[i].p_alpha + p[i+4].p_alpha); 
            double force_mag = helix_strength * (d - 5.5); // 1ターンのピッチ(5.5A)に近づける
            Vec3 f = normalize(dir) * force_mag;
            p[i].force += f;
            p[i+4].force -= f;
        }

        // 4. 疎水性引力（水から逃げる力）と 斥力（弾き返す力）
        for (int i = 0; i < N; i++) {
            for (int j = i + 5; j < N; j++) {
                Vec3 dir = p[j].pos - p[i].pos;
                double d = length(dir);
                if (d < 0.1) d = 0.1;

                double force_mag = 0.0;
                
                // 斥力：4.0A以内に近づいたら強烈に弾き返す
                if (d < MIN_DIST) {
                    force_mag = -50.0 * std::pow((MIN_DIST - d), 2); 
                } 
                // 疎水性引力：遠くにあれば、水から逃げて集まろうとする
                else if (d < 20.0 && p[i].hydro > 0 && p[j].hydro > 0) {
                    // 疎水性が高いほど強く引く
                    force_mag = 5.0 * (p[i].hydro + p[j].hydro) / (d * d); 
                }

                // システイン(C)のシャペロンワープ（例外的な絶対引力）
                if (p[i].code == 'C' && p[j].code == 'C' && d > 4.5) {
                    force_mag = 10.0; // どんなに離れていても一定の力で引き寄せる
                }

                Vec3 f = normalize(dir) * force_mag;
                p[i].force += f;
                p[j].force -= f;
            }
        }

        // 5. ニュートンの運動方程式（座標の更新と冷却）
        for (int i = 0; i < N; i++) {
            p[i].vel += p[i].force * DT;
            p[i].vel = p[i].vel * DAMPING; // 水の摩擦で動きを落ち着かせる
            p[i].pos += p[i].vel * DT;
        }
    }

    // --- 擬似フルアトム出力（主鎖のみ） ---
    std::ofstream out("Universe_v32_PhysicsEngine.pdb");
    out << "HEADER    PROJECT-137 V32.0 EMERGENCE PHYSICS ENGINE\n";
    int atomId = 1;
    for (int i=0; i<N; i++) {
        char buf[100];
        sprintf(buf, "ATOM  %5d  CA  %3s A%4d    %8.3f%8.3f%8.3f  1.00 20.00           C\n", 
                atomId++, p[i].name.c_str(), p[i].id, p[i].pos.x, p[i].pos.y, p[i].pos.z); 
        out << buf;
    }
    out << "END\n";
    std::cout << "[Complete] Check 'Universe_v32_PhysicsEngine.pdb'!\n";
    return 0;
}
