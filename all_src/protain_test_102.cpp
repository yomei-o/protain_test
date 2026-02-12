#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <cmath>
#include <map>
#include <algorithm>
#include <random>

// --- Project-137 v50.1 [Final Golden Master] ---
// 最終調整: Alpha-Helix形成力の暴走を防ぐ「ソフトランディング(力の上限)」を追加。
// これにより、カルモジュリンの優雅な巻きと、コラーゲンの剛直な伸びを両立。
// 生命シミュレーションの集大成。

const double CA_DIST = 3.8;   
const double DT = 0.01;       
const double DAMPING = 0.85;  
const double MAX_VEL = 2.5;   
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

struct AAProps { std::string name; double hydro; double sc_len; double sc_radius; bool is_helix_breaker; };

AAProps getAAProps(char code) {
    static std::map<char, AAProps> m = {
        {'G',{"GLY",-0.4, 0.00, 1.7, true}}, // Helix Breaker
        {'P',{"PRO",-1.6, 1.50, 1.6, true}}, // Helix Breaker
        {'A',{"ALA", 1.8, 1.50, 1.8, false}}, 
        {'V',{"VAL", 4.2, 2.00, 2.2, false}}, 
        {'L',{"LEU", 3.8, 2.50, 2.5, false}}, 
        {'I',{"ILE", 4.5, 2.50, 2.5, false}}, 
        {'F',{"PHE", 2.8, 3.00, 2.8, false}}, 
        {'Y',{"TYR",-1.3, 3.20, 2.9, false}}, 
        {'W',{"TRP",-0.9, 3.50, 3.2, false}}, 
        {'S',{"SER",-0.8, 1.80, 1.9, false}}, 
        {'T',{"THR",-0.7, 2.00, 2.1, false}}, 
        {'C',{"CYS", 2.5, 1.80, 1.9, false}}, 
        {'M',{"MET", 1.9, 3.00, 2.7, false}}, 
        {'H',{"HIS",-3.2, 2.80, 2.4, false}}, 
        {'D',{"ASP",-3.5, 2.20, 2.1, false}}, 
        {'E',{"GLU",-3.5, 2.80, 2.4, false}}, 
        {'N',{"ASN",-3.5, 2.20, 2.1, false}}, 
        {'Q',{"GLN",-3.5, 2.80, 2.4, false}}, 
        {'K',{"LYS",-3.9, 3.50, 2.6, false}}, 
        {'R',{"ARG",-4.5, 4.00, 2.7, false}}
    };
    if(m.count(code)) return m[code];
    return {"UNK",0, 1.5, 2.0, false};
}

struct Particle {
    int id; char code; std::string name;
    double hydro, sc_len, sc_radius;
    bool is_breaker;
    Vec3 pos, vel, force;
    Vec3 sc_pos; 
    bool active;
    int spawn_age;
};

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
    if (d < 0.1) d = 0.1;

    double force_mag = 0.0;
    double contact_dist = radius1 + radius2;

    if (d < contact_dist) {
        force_mag = -300.0 * (contact_dist - d); 
    }
    else if (useHydro && d < 12.0) {
        if (p1.hydro > 0 && p2.hydro > 0) {
            force_mag = 80.0 * (p1.hydro + p2.hydro) / d; 
        } else if (p1.hydro < -1.0 && p2.hydro < -1.0) {
            force_mag = -20.0 / (d * d); 
        } else {
            force_mag = -5.0 / (d * d); 
        }
    }

    if (force_mag > 500.0) force_mag = 500.0;
    if (force_mag < -500.0) force_mag = -500.0;

    Vec3 f = normalize(dir) * force_mag;
    p1.force += f;
    p2.force -= f;
}

int main() {
    std::string rawInput;
    std::cout << "--- Project-137 v50.1 [FINAL] ---" << std::endl;
    std::cout << "Sequence: ";
    std::cin >> rawInput;
    std::transform(rawInput.begin(), rawInput.end(), rawInput.begin(), ::toupper);

    std::string aaSeq = rawInput; 
    std::vector<Particle> p(aaSeq.length());
    for(int i=0; i<(int)aaSeq.length(); ++i) {
        AAProps props = getAAProps(aaSeq[i]);
        p[i] = {i+1, aaSeq[i], props.name, props.hydro, props.sc_len, props.sc_radius, props.is_helix_breaker,
                {0,0,0}, {0,0,0}, {0,0,0}, {0,0,0}, false, 0};
    }

    std::mt19937 gen(137);
    std::uniform_real_distribution<> dis(-1.0, 1.0);

    int N = p.size();
    int TOTAL_STEPS = (N * STEPS_PER_AA) + 12000; // 十分な時間を確保
    int spawned_count = 0;

    for (int step = 0; step < TOTAL_STEPS; step++) {
        
        // --- 1. Kinetic Spawning ---
        if (step % STEPS_PER_AA == 0 && spawned_count < N) {
            p[spawned_count].active = true;
            Vec3 spawn_pos = {0,0,0};
            if (spawned_count == 0) spawn_pos = {0,0,0};
            else if (spawned_count == 1) spawn_pos = p[0].pos + Vec3{3.8, 0, 0};
            else {
                Vec3 prev = p[spawned_count-1].pos;
                Vec3 prev2 = p[spawned_count-2].pos;
                Vec3 direction = normalize(prev - prev2);
                direction.x += dis(gen) * 0.2; direction.y += dis(gen) * 0.2; direction.z += dis(gen) * 0.2;
                spawn_pos = prev + normalize(direction) * 3.8;
            }
            p[spawned_count].pos = spawn_pos;
            if (spawned_count > 0) p[spawned_count].vel = normalize(p[spawned_count].pos - p[spawned_count-1].pos) * 0.5;
            spawned_count++;
        }

        // 年齢更新
        for(int i=0; i<spawned_count; i++) p[i].spawn_age++;

        // 側鎖位置更新
        for (int i = 1; i < spawned_count - 1; i++) {
            Vec3 dir = computeSideChainDir(p[i-1], p[i], p[i+1]);
            p[i].sc_pos = p[i].pos + dir * p[i].sc_len;
        }
        if (spawned_count > 1) {
            p[0].sc_pos = p[0].pos + normalize(p[0].pos - p[1].pos) * p[0].sc_len;
            int last = spawned_count - 1;
            Vec3 dir = normalize(p[last].pos - p[last-1].pos);
            Vec3 up = {0,1,0}; if(std::abs(dir.y)>0.9) up={1,0,0};
            p[last].sc_pos = p[last].pos + normalize(cross(dir, up)) * p[last].sc_len;
        }

        for (int i = 0; i < spawned_count; i++) p[i].force = {0,0,0};

        // 2. 熱揺らぎ
        for (int i = 0; i < spawned_count; i++) {
            p[i].force += Vec3{dis(gen), dis(gen), dis(gen)} * BASE_TEMP;
        }

        // 3. 結合バネ (Backbone Spring)
        for (int i = 0; i < spawned_count - 1; i++) {
            Vec3 dir = p[i+1].pos - p[i].pos;
            double d = length(dir);
            double force_mag = 250.0 * (d - CA_DIST); 
            Vec3 f = normalize(dir) * force_mag;
            p[i].force += f;
            p[i+1].force -= f;
        }

        // 4. 剛性制御 (Exoskeleton / Angle Spring)
        for (int i = 0; i < spawned_count - 2; i++) {
            Vec3 dir = p[i+2].pos - p[i].pos;
            double d = length(dir);
            bool is_rigid = (p[i+1].code == 'P' || p[i+1].code == 'G');
            double target = is_rigid ? 6.5 : 5.4; 
            double stiffness = is_rigid ? 500.0 : 40.0;
            double force_mag = stiffness * (d - target); 
            Vec3 f = normalize(dir) * force_mag;
            p[i].force += f;
            p[i+2].force -= f;
        }

        // 5. Alpha-Helix Maker (Hydrogen Bond Emulation)
        for (int i = 0; i < spawned_count - 4; i++) {
            bool helix_compatible = true;
            for(int k=0; k<=4; k++) if (p[i+k].is_breaker) { helix_compatible = false; break; }

            if (helix_compatible) {
                Vec3 dir = p[i+4].pos - p[i].pos;
                double d = length(dir);
                double target = 6.0; 
                double force_mag = 100.0 * (d - target);
                
                // ★修正: Helix力にも安全弁をつける
                if (force_mag > 200.0) force_mag = 200.0;
                if (force_mag < -200.0) force_mag = -200.0;

                Vec3 f = normalize(dir) * force_mag;
                p[i].force += f;
                p[i+4].force -= f;
            }
        }

        // 6. 衝突・相互作用 (Solid Physics)
        for (int i = 0; i < spawned_count; i++) {
            for (int j = i + 3; j < spawned_count; j++) {
                // CA vs CA
                applyForce(p[i], p[j], p[i].pos, p[j].pos, 1.9, 1.9, false);
                // SC vs SC
                if (p[i].sc_len > 0.1 && p[j].sc_len > 0.1) {
                    applyForce(p[i], p[j], p[i].sc_pos, p[j].sc_pos, p[i].sc_radius, p[j].sc_radius, true);
                }
                // SC vs CA
                if (p[i].sc_len > 0.1) applyForce(p[i], p[j], p[i].sc_pos, p[j].pos, p[i].sc_radius, 1.9, false);
                if (p[j].sc_len > 0.1) applyForce(p[i], p[j], p[i].pos, p[j].sc_pos, 1.9, p[j].sc_radius, false);
            }
        }

        // 7. 積分 (Soft Landing Damping)
        for (int i = 0; i < spawned_count; i++) {
            p[i].vel += p[i].force * DT;
            if (length(p[i].vel) > MAX_VEL) p[i].vel = normalize(p[i].vel) * MAX_VEL;
            
            // 若い原子は強く減衰させ、暴れを防ぐ
            double damp = (p[i].spawn_age < 500) ? 0.6 : DAMPING;
            p[i].vel = p[i].vel * damp; 
            p[i].pos += p[i].vel * DT;
        }
    }

    std::ofstream out("Universe_v50_1_FINAL.pdb");
    int atomId = 1;
    for (int i=0; i<N; i++) {
        char buf[100];
        // B-factor欄(最後)にspawn_ageを入れると、後で「どの順でできたか」が見えて面白い
        sprintf(buf, "ATOM  %5d  CA  %3s A%4d    %8.3f%8.3f%8.3f  1.00%6.2f           C\n", 
                atomId++, p[i].name.c_str(), p[i].id, p[i].pos.x, p[i].pos.y, p[i].pos.z, (double)p[i].id); 
        out << buf;
        if (p[i].sc_len > 0.1) {
            sprintf(buf, "ATOM  %5d  CB  %3s A%4d    %8.3f%8.3f%8.3f  1.00%6.2f           C\n", 
                    atomId++, p[i].name.c_str(), p[i].id, p[i].sc_pos.x, p[i].sc_pos.y, p[i].sc_pos.z, (double)p[i].id); 
            out << buf;
        }
    }
    out << "END\n";
    std::cout << "[Complete] Universe_v50_1_FINAL.pdb generated." << std::endl;
    return 0;
}

