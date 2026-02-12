#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <cmath>
#include <map>
#include <algorithm>
#include <random>

// --- Project-137 v48.0 [Kinetic Spawner] ---
// 最終安全装置:
// 1. 「ベクトル延長スポーン」で、自分の身体の中に原子が出る自爆事故を回避。
// 2. 物理演算の力が暴走しないよう、フォース・クランプ(力の上限)を導入。
// 3. 側鎖の計算を、端っこの原子でもバグらないように補正。

const double CA_DIST = 3.8;   
const double DT = 0.01;       
const double DAMPING = 0.85;  
const double MAX_VEL = 3.0;   
const double BASE_TEMP = 0.3; 
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

struct AAProps { std::string name; double hydro; double sc_len; double sc_radius; };

AAProps getAAProps(char code) {
    static std::map<char, AAProps> m = {
        {'G',{"GLY",-0.4, 0.00, 1.7}}, 
        {'A',{"ALA", 1.8, 1.50, 1.8}}, 
        {'V',{"VAL", 4.2, 2.00, 2.2}}, 
        {'L',{"LEU", 3.8, 2.50, 2.5}}, 
        {'I',{"ILE", 4.5, 2.50, 2.5}}, 
        {'F',{"PHE", 2.8, 3.00, 2.8}}, 
        {'Y',{"TYR",-1.3, 3.20, 2.9}}, 
        {'W',{"TRP",-0.9, 3.50, 3.2}}, 
        {'S',{"SER",-0.8, 1.80, 1.9}}, 
        {'T',{"THR",-0.7, 2.00, 2.1}}, 
        {'C',{"CYS", 2.5, 1.80, 1.9}}, 
        {'M',{"MET", 1.9, 3.00, 2.7}}, 
        {'P',{"PRO",-1.6, 1.50, 1.6}}, 
        {'H',{"HIS",-3.2, 2.80, 2.4}}, 
        {'D',{"ASP",-3.5, 2.20, 2.1}}, 
        {'E',{"GLU",-3.5, 2.80, 2.4}}, 
        {'N',{"ASN",-3.5, 2.20, 2.1}}, 
        {'Q',{"GLN",-3.5, 2.80, 2.4}}, 
        {'K',{"LYS",-3.9, 3.50, 2.6}}, 
        {'R',{"ARG",-4.5, 4.00, 2.7}}
    };
    if(m.count(code)) return m[code];
    return {"UNK",0, 1.5, 2.0};
}

struct Particle {
    int id; char code; std::string name;
    double hydro, sc_len, sc_radius;
    Vec3 pos, vel, force;
    Vec3 sc_pos; 
    bool active;
};

// 安全な側鎖方向計算
Vec3 computeSideChainDir(const Particle& p_prev, const Particle& p_curr, const Particle& p_next) {
    Vec3 v1 = normalize(p_curr.pos - p_prev.pos);
    Vec3 v2 = normalize(p_next.pos - p_curr.pos);
    Vec3 bisector = normalize(v1 - v2); 
    if (length(bisector) < 0.1) {
        Vec3 arbitrary = {0,1,0};
        if (std::abs(v1.y) > 0.9) arbitrary = {1,0,0};
        bisector = normalize(cross(v1, arbitrary));
    }
    return bisector;
}

void applyForce(Particle& p1, Particle& p2, Vec3 pos1, Vec3 pos2, double radius1, double radius2, bool useHydro) {
    Vec3 dir = pos2 - pos1;
    double d = length(dir);
    if (d < 0.1) d = 0.1; // ゼロ除算防止

    double force_mag = 0.0;
    double contact_dist = radius1 + radius2;

    if (d < contact_dist) {
        force_mag = -250.0 * (contact_dist - d); // 斥力
    }
    else if (useHydro && d < 12.0) {
        if (p1.hydro > 0 && p2.hydro > 0) {
            force_mag = 80.0 * (p1.hydro + p2.hydro) / d; // 引力
        } else if (p1.hydro < -1.0 && p2.hydro < -1.0) {
            force_mag = -20.0 / (d * d); // 電荷反発
        }
    }

    // ★安全弁: 力が強すぎたらカットする
    if (force_mag > 500.0) force_mag = 500.0;
    if (force_mag < -500.0) force_mag = -500.0;

    Vec3 f = normalize(dir) * force_mag;
    p1.force += f;
    p2.force -= f;
}

int main() {
    std::string rawInput;
    std::cout << "--- Project-137 v48.0 [Kinetic Spawner] ---" << std::endl;
    std::cout << "Sequence: ";
    std::cin >> rawInput;
    std::transform(rawInput.begin(), rawInput.end(), rawInput.begin(), ::toupper);

    std::string aaSeq = rawInput; 
    std::vector<Particle> p(aaSeq.length());
    for(int i=0; i<(int)aaSeq.length(); ++i) {
        AAProps props = getAAProps(aaSeq[i]);
        p[i] = {i+1, aaSeq[i], props.name, props.hydro, props.sc_len, props.sc_radius,
                {0,0,0}, {0,0,0}, {0,0,0}, {0,0,0}, false};
    }

    std::mt19937 gen(137);
    std::uniform_real_distribution<> dis(-1.0, 1.0);

    int N = p.size();
    int TOTAL_STEPS = (N * STEPS_PER_AA) + 8000; 
    int spawned_count = 0;

    for (int step = 0; step < TOTAL_STEPS; step++) {
        
        // --- ★改良版: 運動学的スポーン ---
        if (step % STEPS_PER_AA == 0 && spawned_count < N) {
            p[spawned_count].active = true;
            
            Vec3 spawn_pos = {0,0,0};
            if (spawned_count == 0) {
                spawn_pos = {0,0,0};
            } else if (spawned_count == 1) {
                spawn_pos = p[0].pos + Vec3{3.8, 0, 0};
            } else {
                // 1個前と2個前の原子を結んだ線を延長する
                Vec3 prev = p[spawned_count-1].pos;
                Vec3 prev2 = p[spawned_count-2].pos;
                Vec3 direction = normalize(prev - prev2);
                // 少しランダムに振って、一直線になりすぎるのを防ぐ
                direction.x += dis(gen) * 0.2;
                direction.y += dis(gen) * 0.2;
                direction.z += dis(gen) * 0.2;
                spawn_pos = prev + normalize(direction) * 3.8;
            }
            p[spawned_count].pos = spawn_pos;
            
            // 初速も進行方向に与えて、逆戻りを防ぐ
            if (spawned_count > 0) {
                p[spawned_count].vel = normalize(p[spawned_count].pos - p[spawned_count-1].pos) * 0.5;
            }
            
            spawned_count++;
        }

        // --- 側鎖位置計算 ---
        for (int i = 1; i < spawned_count - 1; i++) {
            Vec3 dir = computeSideChainDir(p[i-1], p[i], p[i+1]);
            p[i].sc_pos = p[i].pos + dir * p[i].sc_len;
        }
        // 端っこの処理 (バグ防止)
        if (spawned_count > 1) {
            // 先頭: 2番目の向きを逆にする
            p[0].sc_pos = p[0].pos + normalize(p[0].pos - p[1].pos) * p[0].sc_len;
            // 末尾: 1つ前の進行方向をそのまま使う
            int last = spawned_count - 1;
            Vec3 dir = normalize(p[last].pos - p[last-1].pos);
            // 進行方向と垂直な適当なベクトルを作る
            Vec3 up = {0,1,0}; if(std::abs(dir.y)>0.9) up={1,0,0};
            Vec3 sc_dir = normalize(cross(dir, up));
            p[last].sc_pos = p[last].pos + sc_dir * p[last].sc_len;
        }

        for (int i = 0; i < spawned_count; i++) p[i].force = {0,0,0};

        // 1. 熱揺らぎ
        for (int i = 0; i < spawned_count; i++) {
            p[i].force += Vec3{dis(gen), dis(gen), dis(gen)} * BASE_TEMP;
        }

        // 2. 結合バネ
        for (int i = 0; i < spawned_count - 1; i++) {
            Vec3 dir = p[i+1].pos - p[i].pos;
            double d = length(dir);
            double force_mag = 250.0 * (d - CA_DIST); 
            Vec3 f = normalize(dir) * force_mag;
            p[i].force += f;
            p[i+1].force -= f;
        }

        // 3. 角度維持
        for (int i = 0; i < spawned_count - 2; i++) {
            Vec3 dir = p[i+2].pos - p[i].pos;
            double d = length(dir);
            double target = (p[i+1].code == 'P' || p[i+1].code == 'G') ? 6.0 : 5.2;
            double force_mag = 50.0 * (d - target); 
            Vec3 f = normalize(dir) * force_mag;
            p[i].force += f;
            p[i+2].force -= f;
        }

        // 4. 二重装甲 + 安全弁付き相互作用
        for (int i = 0; i < spawned_count; i++) {
            for (int j = i + 3; j < spawned_count; j++) {
                
                // (A) CA vs CA
                applyForce(p[i], p[j], p[i].pos, p[j].pos, 1.9, 1.9, false); // CA半径1.9

                // (B) SC vs SC
                if (p[i].sc_len > 0.1 && p[j].sc_len > 0.1) {
                    applyForce(p[i], p[j], p[i].sc_pos, p[j].sc_pos, p[i].sc_radius, p[j].sc_radius, true);
                }

                // (C) SC vs CA
                if (p[i].sc_len > 0.1) {
                    applyForce(p[i], p[j], p[i].sc_pos, p[j].pos, p[i].sc_radius, 1.9, false);
                }
                if (p[j].sc_len > 0.1) {
                    applyForce(p[i], p[j], p[i].pos, p[j].sc_pos, 1.9, p[j].sc_radius, false);
                }
            }
        }

        // 5. 積分
        for (int i = 0; i < spawned_count; i++) {
            p[i].vel += p[i].force * DT;
            if (length(p[i].vel) > MAX_VEL) p[i].vel = normalize(p[i].vel) * MAX_VEL;
            p[i].vel = p[i].vel * DAMPING; 
            p[i].pos += p[i].vel * DT;
        }
    }

    std::ofstream out("Universe_v48_KineticSpawner.pdb");
    int atomId = 1;
    for (int i=0; i<N; i++) {
        char buf[100];
        sprintf(buf, "ATOM  %5d  CA  %3s A%4d    %8.3f%8.3f%8.3f  1.00 20.00           C\n", 
                atomId++, p[i].name.c_str(), p[i].id, p[i].pos.x, p[i].pos.y, p[i].pos.z); 
        out << buf;

        if (p[i].sc_len > 0.1) {
            sprintf(buf, "ATOM  %5d  CB  %3s A%4d    %8.3f%8.3f%8.3f  1.00 20.00           C\n", 
                    atomId++, p[i].name.c_str(), p[i].id, p[i].sc_pos.x, p[i].sc_pos.y, p[i].sc_pos.z); 
            out << buf;
        }
    }
    out << "END\n";
    std::cout << "[Complete] Universe_v48_KineticSpawner.pdb generated." << std::endl;
    return 0;
}
