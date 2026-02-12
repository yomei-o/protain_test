#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <cmath>
#include <map>
#include <algorithm>
#include <random>
#include <iomanip>

// --- Project-137 v53.0 [Zenith] ---
// 最終実験のための「頂点」バージョン。
// 1. S-S結合時の「立体障害緩和」: システイン間の結合をスムーズにするため、結合時のみ半径を縮小。
// 2. Exponential Cooling: 指数関数的な冷却により、より安定した最良構造へ導く。
// 3. UX向上: 進捗バーの表示で、計算の進行を可視化。

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
    bool isNaN() const { return std::isnan(x) || std::isnan(y) || std::isnan(z); }
};

Vec3 cross(Vec3 a, Vec3 b) { return {a.y*b.z - a.z*b.y, a.z*b.x - a.x*b.z, a.x*b.y - a.y*b.x}; }
double dot(Vec3 a, Vec3 b) { return a.x*b.x + a.y*b.y + a.z*b.z; }
double length(Vec3 v) { return std::sqrt(v.x*v.x + v.y*v.y + v.z*v.z); } 
Vec3 normalize(Vec3 v) { double l = length(v); return (l>1e-6)? v * (1.0/l) : Vec3{0,0,0}; }

struct AAProps { std::string name; double hydro; double sc_len; double sc_radius; bool is_helix_breaker; };

AAProps getAAProps(char code) {
    static std::map<char, AAProps> m = {
        {'G',{"GLY",-0.4, 0.00, 1.7, true}}, 
        {'P',{"PRO",-1.6, 1.50, 1.6, true}}, 
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
    
    // ★Zenith修正: S-S結合するペアなら、衝突半径を小さくして「抱擁」を許可する
    bool is_SS_pair = (p1.code == 'C' && p2.code == 'C' && d < 12.0);
    double effective_contact = is_SS_pair ? (radius1 + radius2) * 0.6 : (radius1 + radius2);

    // 1. 衝突 (Hard Core)
    if (d < effective_contact) {
        force_mag = -300.0 * (effective_contact - d); 
    }
    // 2. ジスルフィド結合
    else if (is_SS_pair) {
        // ターゲット距離 3.0A (CA-CA間距離として妥当な値)へ強く引く
        force_mag = 120.0 * (d - 3.0) / d; 
    }
    // 3. 疎水性・親水性・中性
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
    std::cout << "--- Project-137 v53.0 [Zenith] ---" << std::endl;
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
    int TOTAL_STEPS = (N * STEPS_PER_AA) + 30000; 
    int spawned_count = 0;

    std::cout << "Simulating " << N << " residues..." << std::endl;

    for (int step = 0; step < TOTAL_STEPS; step++) {
        
        // 進捗表示
        if (step % (TOTAL_STEPS / 20) == 0) {
            std::cout << "\rProgress: [" << std::string(step * 20 / TOTAL_STEPS, '#') 
                      << std::string(20 - (step * 20 / TOTAL_STEPS), ' ') << "] " 
                      << int((double)step/TOTAL_STEPS*100) << "%" << std::flush;
        }

        // Kinetic Spawning
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

        for(int i=0; i<spawned_count; i++) p[i].spawn_age++;

        // Side Chain Update
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

        // 1. Cooling (Exponential Decay)
        // 終盤にかけて急激に冷やすことで、構造をより強く固定する
        double progress = (double)step / TOTAL_STEPS;
        double current_temp = BASE_TEMP * std::exp(-5.0 * std::pow(std::max(0.0, progress - 0.5), 2)); 
        if (step > TOTAL_STEPS - 2000) current_temp = 0.0;

        for (int i = 0; i < spawned_count; i++) {
            p[i].force += Vec3{dis(gen), dis(gen), dis(gen)} * current_temp;
        }

        // 2. Backbone Springs
        for (int i = 0; i < spawned_count - 1; i++) {
            Vec3 dir = p[i+1].pos - p[i].pos;
            double d = length(dir);
            double force_mag = 250.0 * (d - CA_DIST); 
            Vec3 f = normalize(dir) * force_mag;
            p[i].force += f;
            p[i+1].force -= f;
        }

        // 3. Exoskeleton (Collagen Rigidity)
        for (int i = 0; i < spawned_count - 2; i++) {
            Vec3 dir = p[i+2].pos - p[i].pos;
            double d = length(dir);
            bool is_rigid = (p[i+1].code == 'P' || p[i+1].code == 'G');
            double target = is_rigid ? 6.5 : 5.4; 
            double stiffness = is_rigid ? 1000.0 : 40.0;
            double force_mag = stiffness * (d - target); 
            Vec3 f = normalize(dir) * force_mag;
            p[i].force += f;
            p[i+2].force -= f;
        }

        // 4. Helix Maker
        for (int i = 0; i < spawned_count - 4; i++) {
            bool helix_compatible = true;
            for(int k=0; k<=4; k++) if (p[i+k].is_breaker) { helix_compatible = false; break; }

            if (helix_compatible) {
                Vec3 dir = p[i+4].pos - p[i].pos;
                double d = length(dir);
                double force_mag = 200.0 * (d - 6.0); 
                if (force_mag > 200.0) force_mag = 200.0;
                if (force_mag < -200.0) force_mag = -200.0;
                Vec3 f = normalize(dir) * force_mag;
                p[i].force += f;
                p[i+4].force -= f;
            }
        }

        // 5. Interactions (Collision & S-S)
        for (int i = 0; i < spawned_count; i++) {
            for (int j = i + 3; j < spawned_count; j++) {
                applyForce(p[i], p[j], p[i].pos, p[j].pos, 1.9, 1.9, false);
                if (p[i].sc_len > 0.1 && p[j].sc_len > 0.1) 
                    applyForce(p[i], p[j], p[i].sc_pos, p[j].sc_pos, p[i].sc_radius, p[j].sc_radius, true);
                if (p[i].sc_len > 0.1) 
                    applyForce(p[i], p[j], p[i].sc_pos, p[j].pos, p[i].sc_radius, 1.9, false);
                if (p[j].sc_len > 0.1) 
                    applyForce(p[i], p[j], p[i].pos, p[j].sc_pos, 1.9, p[j].sc_radius, false);
            }
        }

        // 6. Integration & NaN Guard
        for (int i = 0; i < spawned_count; i++) {
            p[i].vel += p[i].force * DT;
            if (length(p[i].vel) > MAX_VEL) p[i].vel = normalize(p[i].vel) * MAX_VEL;
            double damp = (p[i].spawn_age < 500) ? 0.6 : DAMPING;
            p[i].vel = p[i].vel * damp; 
            
            Vec3 next_pos = p[i].pos + p[i].vel * DT;
            if (!next_pos.isNaN()) {
                p[i].pos = next_pos;
            } else {
                p[i].vel = {0,0,0}; 
            }
        }
    }
    std::cout << "\rProgress: [####################] 100%" << std::endl;

    std::ofstream out("Universe_v53_Zenith.pdb");
    int atomId = 1;
    for (int i=0; i<N; i++) {
        char buf[100];
        // B-factorに半径を出力
        sprintf(buf, "ATOM  %5d  CA  %3s A%4d    %8.3f%8.3f%8.3f  1.00%6.2f           C\n", 
                atomId++, p[i].name.c_str(), p[i].id, p[i].pos.x, p[i].pos.y, p[i].pos.z, p[i].sc_radius); 
        out << buf;
        if (p[i].sc_len > 0.1) {
            sprintf(buf, "ATOM  %5d  CB  %3s A%4d    %8.3f%8.3f%8.3f  1.00%6.2f           C\n", 
                    atomId++, p[i].name.c_str(), p[i].id, p[i].sc_pos.x, p[i].sc_pos.y, p[i].sc_pos.z, p[i].sc_radius); 
            out << buf;
        }
    }
    out << "END\n";
    std::cout << "[Complete] Universe_v53_Zenith.pdb generated." << std::endl;
    return 0;
}

