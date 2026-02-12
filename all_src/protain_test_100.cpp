#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <cmath>
#include <map>
#include <algorithm>
#include <random>

// --- Project-137 v46.0 [Double Armor] ---
// 修正点:
// 1. グリシンの半径を物理的に正しい値(1.7A)に修正し、無限圧縮を回避。
// 2. 「側鎖同士」だけでなく「主鎖(CA)同士」の衝突判定を追加し、すり抜けを防止。
// 3. 側鎖の向き計算を安定化。

const double CA_DIST = 3.8;   
const double DT = 0.01;       
const double DAMPING = 0.80;  
const double MAX_VEL = 2.5;   
const double BASE_TEMP = 0.5;
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

// sc_len: CA-CB距離, sc_radius: 側鎖の衝突半径
struct AAProps { std::string name; double hydro; double sc_len; double sc_radius; };

AAProps getAAProps(char code) {
    static std::map<char, AAProps> m = {
        // ★修正: Gの半径を0.5->1.7へ。これで潰れない。
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
        {'P',{"PRO",-1.6, 1.50, 2.1}}, 
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

// 側鎖の向き計算（安定版）
Vec3 computeSideChainDir(const Particle& p_prev, const Particle& p_curr, const Particle& p_next) {
    Vec3 v1 = normalize(p_curr.pos - p_prev.pos);
    Vec3 v2 = normalize(p_next.pos - p_curr.pos);
    Vec3 bisector = normalize(v1 - v2); // 外側へ向くベクトル
    
    // 一直線に近い場合(外積で上方向を決める)
    if (length(bisector) < 0.1) {
        Vec3 arbitrary = {0,1,0};
        if (std::abs(v1.y) > 0.9) arbitrary = {1,0,0};
        bisector = normalize(cross(v1, arbitrary));
    }
    return bisector;
}

int main() {
    std::string rawInput;
    std::cout << "--- Project-137 v46.0 [Double Armor] ---" << std::endl;
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
    int TOTAL_STEPS = (N * STEPS_PER_AA) + 5000; 
    int spawned_count = 0;

    for (int step = 0; step < TOTAL_STEPS; step++) {
        
        // リボソーム射出
        if (step % STEPS_PER_AA == 0 && spawned_count < N) {
            p[spawned_count].active = true;
            if (spawned_count > 0) {
                Vec3 prev_pos = p[spawned_count-1].pos;
                p[spawned_count].pos = prev_pos + Vec3{3.5, dis(gen)*0.5, dis(gen)*0.5}; 
            }
            spawned_count++;
        }

        // 側鎖位置更新
        for (int i = 1; i < spawned_count - 1; i++) {
            Vec3 dir = computeSideChainDir(p[i-1], p[i], p[i+1]);
            p[i].sc_pos = p[i].pos + dir * p[i].sc_len;
        }
        if (spawned_count > 0) p[0].sc_pos = p[0].pos + Vec3{0, 1, 0} * p[0].sc_len;
        if (spawned_count > 1) p[spawned_count-1].sc_pos = p[spawned_count-1].pos + Vec3{0, 1, 0} * p[spawned_count-1].sc_len;

        for (int i = 0; i < spawned_count; i++) p[i].force = {0,0,0};

        // 1. 熱揺らぎ
        for (int i = 0; i < spawned_count; i++) {
            p[i].force += Vec3{dis(gen), dis(gen), dis(gen)} * BASE_TEMP;
        }

        // 2. 結合バネ (CA-CA)
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

        // 4. ★二重装甲相互作用★
        for (int i = 0; i < spawned_count; i++) {
            for (int j = i + 3; j < spawned_count; j++) {
                
                // (A) 主鎖(CA)同士の衝突判定 [NEW!]
                // これがないと背骨同士がすり抜ける
                Vec3 caDir = p[j].pos - p[i].pos;
                double caDist = length(caDir);
                if (caDist < 3.5) { // CAのファンデルワールス直径和に近い値
                    Vec3 f = normalize(caDir) * (-200.0 * (3.5 - caDist));
                    p[i].force += f;
                    p[j].force -= f;
                }

                // (B) 側鎖(SC)同士の判定
                Vec3 posI = (p[i].sc_len > 0.1) ? p[i].sc_pos : p[i].pos;
                Vec3 posJ = (p[j].sc_len > 0.1) ? p[j].sc_pos : p[j].pos;
                Vec3 scDir = posJ - posI;
                double scDist = length(scDir);
                if (scDist < 0.1) scDist = 0.1;

                double force_mag = 0.0;
                double contact_dist = p[i].sc_radius + p[j].sc_radius;

                // 斥力
                if (scDist < contact_dist) {
                    force_mag = -200.0 * (contact_dist - scDist);
                }
                // 疎水性引力
                else if (scDist < 12.0 && p[i].hydro > 0 && p[j].hydro > 0) {
                    force_mag = 40.0 * (p[i].hydro + p[j].hydro) / scDist;
                }
                // 親水性反発
                else if (scDist < 15.0 && p[i].hydro < -1.0 && p[j].hydro < -1.0) {
                    force_mag = -20.0 / (scDist * scDist);
                }

                Vec3 f = normalize(scDir) * force_mag;
                p[i].force += f;
                p[j].force -= f;
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

    std::ofstream out("Universe_v46_DoubleArmor.pdb");
    int atomId = 1;
    for (int i=0; i<N; i++) {
        // CA出力
        char buf[100];
        sprintf(buf, "ATOM  %5d  CA  %3s A%4d    %8.3f%8.3f%8.3f  1.00 20.00           C\n", 
                atomId++, p[i].name.c_str(), p[i].id, p[i].pos.x, p[i].pos.y, p[i].pos.z); 
        out << buf;

        // CB出力
        if (p[i].sc_len > 0.1) {
            sprintf(buf, "ATOM  %5d  CB  %3s A%4d    %8.3f%8.3f%8.3f  1.00 20.00           C\n", 
                    atomId++, p[i].name.c_str(), p[i].id, p[i].sc_pos.x, p[i].sc_pos.y, p[i].sc_pos.z); 
            out << buf;
        }
    }
    out << "END\n";
    std::cout << "[Complete] Universe_v46_DoubleArmor.pdb generated." << std::endl;
    return 0;
}

