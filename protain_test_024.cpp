#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <cmath>
#include <iomanip>
#include <map>

// --- Project-137: Universe OS v13.0 "Full-Atom & Ribbon" ---
// コンセプト: CA座標からN, C, Oの原子座標を逆算して「ペプチド平面」を生成。
// さらに、PDBヘッダに「HELIXレコード」を自動記述し、美しいリボンを描画させる。

const double PI = 3.1415926535;
const double CA_DIST = 3.8;
const double MIN_DIST = 4.0;

struct Vec3 {
    double x, y, z;
    Vec3 operator-(const Vec3& o) const { return {x - o.x, y - o.y, z - o.z}; }
    Vec3 operator+(const Vec3& o) const { return {x + o.x, y + o.y, z + o.z}; }
    Vec3 operator*(double s) const { return {x * s, y * s, z * s}; }
};

Vec3 cross(Vec3 a, Vec3 b) {
    return {a.y*b.z - a.z*b.y, a.z*b.x - a.x*b.z, a.x*b.y - a.y*b.x};
}

double length(Vec3 v) { return std::sqrt(v.x*v.x + v.y*v.y + v.z*v.z); }
Vec3 normalize(Vec3 v) { double l = length(v); return (l>1e-6)? Vec3{v.x/l, v.y/l, v.z/l} : Vec3{0,0,0}; }

struct Residue {
    int resSeq;
    char code;
    std::string resName;
    double hydrophobicity;
    double charge;
    bool isHelix; // 螺旋かどうか
    Vec3 ca, n, c, o; // フルアトム座標
};

struct AAProps { std::string name; double hydro; double charge; };
AAProps getAAProps(char code) {
    static std::map<char, AAProps> m = {
        {'I',{"ILE",4.5,0}}, {'V',{"VAL",4.2,0}}, {'L',{"LEU",3.8,0}}, {'F',{"PHE",2.8,0}},
        {'C',{"CYS",2.5,0}}, {'M',{"MET",1.9,0}}, {'A',{"ALA",1.8,0}}, {'G',{"GLY",-0.4,0}},
        {'T',{"THR",-0.7,0}},{'S',{"SER",-0.8,0}},{'W',{"TRP",-0.9,0}},{'Y',{"TYR",-1.3,0}},
        {'P',{"PRO",-1.6,0}},{'H',{"HIS",-3.2,1}}, {'E',{"GLU",-3.5,-1}},{'Q',{"GLN",-3.5,0}},
        {'D',{"ASP",-3.5,-1}},{'N',{"ASN",-3.5,0}},{'K',{"LYS",-3.9,1}}, {'R',{"ARG",-4.5,1}}
    };
    if(m.count(code)) return m[code];
    return {"UNK",0,0};
}

int main() {
    std::string inputSeq;
    std::cout << "Project-137 v13.0 [Full-Atom & Ribbon Viewer]" << std::endl;
    std::cout << "Input: ";
    std::cin >> inputSeq;

    std::vector<Residue> protein;
    for(int i=0; i<(int)inputSeq.length(); ++i) {
        Residue res;
        res.resSeq = i+1; res.code = inputSeq[i];
        AAProps p = getAAProps(res.code);
        res.resName = p.name; res.hydrophobicity = p.hydro; res.charge = p.charge;
        res.isHelix = false;
        protein.push_back(res);
    }

    // 1. CA（骨格）の生成と物理演算 (v11/v12ベース)
    double currX=0, currY=0, currZ=0;
    double angleX=0.6, angleY=0.0;

    for(int i=0; i<(int)protein.size(); ++i) {
        double avgHydro = 0; int count=0;
        for(int k=i-2; k<=i+2; k++) {
            if(k>=0 && k<(int)protein.size()) { avgHydro+=protein[k].hydrophobicity; count++; }
        }
        avgHydro /= (count>0?count:1);

        double localChargeSum = 0;
        for(int k=i-2; k<=i+2; k++) {
            if(k>=0 && k<(int)protein.size()) localChargeSum += std::abs(protein[k].charge);
        }

        if(avgHydro > 1.0) { 
            angleX = angleX*0.8 + 0.6*0.2;
            angleY += 1.0;
            protein[i].isHelix = true; // 螺旋フラグON！
        } else {
            double repulsion = 1.0 + (localChargeSum * 0.2);
            angleX += (PI*(3.0-std::sqrt(5.0))) * 0.5 * repulsion;
            angleY += 0.5 * repulsion;
        }

        Vec3 dir = {std::sin(angleX)*std::cos(angleY), std::sin(angleX)*std::sin(angleY), std::cos(angleX)};
        dir = normalize(dir);
        currX += dir.x*CA_DIST; currY += dir.y*CA_DIST; currZ += dir.z*CA_DIST;
        protein[i].ca = {currX, currY, currZ};
    }

    // 衝突解消 (軽く200回)
    for(int iter=0; iter<200; ++iter) {
        for(int i=0; i<(int)protein.size(); ++i) {
            for(int j=i+1; j<(int)protein.size(); ++j) {
                Vec3 d = protein[j].ca - protein[i].ca;
                double dist = length(d);
                if(dist < MIN_DIST && dist > 0.1) {
                    double push = (MIN_DIST - dist) * 0.5 * 0.1;
                    Vec3 move = normalize(d) * push;
                    protein[i].ca = protein[i].ca - move;
                    protein[j].ca = protein[j].ca + move;
                }
            }
        }
    }

    // 2. フルアトム生成 (N, C, O を捏造してペプチド平面を作る)
    for(int i=0; i<(int)protein.size(); ++i) {
        Vec3 ca_prev = (i > 0) ? protein[i-1].ca : protein[i].ca - Vec3{1,0,0};
        Vec3 ca_next = (i < (int)protein.size()-1) ? protein[i+1].ca : protein[i].ca + Vec3{1,0,0};

        Vec3 v_prev = normalize(protein[i].ca - ca_prev);
        Vec3 v_next = normalize(ca_next - protein[i].ca);

        // NとCをCAの少し手前と奥に配置 (実際の結合距離に近似)
        protein[i].n = protein[i].ca - v_prev * 1.45;
        protein[i].c = protein[i].ca + v_next * 1.52;

        // O(酸素)の向きを決定 (ペプチド平面に垂直な方向へ張り出す)
        Vec3 ortho = cross(v_prev, v_next);
        if(length(ortho) < 0.1) ortho = {0,1,0}; // 直線すぎた場合のフォールバック
        ortho = normalize(ortho);
        protein[i].o = protein[i].c + ortho * 1.23;
    }

    // 3. PDB出力
    std::ofstream out("Universe_v13_Ribbon.pdb");
    out << "HEADER    PROJECT-137 V13.0 RIBBON MASTER\n";

    // 【ハック】HELIXレコードの自動生成
    // 連続してisHelix=trueになっている区間を探し、ビューワーに「ここは螺旋だ」と教える
    int helixId = 1;
    for(int i=0; i<(int)protein.size(); i++) {
        if(protein[i].isHelix) {
            int start = i;
            while(i < (int)protein.size() && protein[i].isHelix) i++;
            int end = i - 1;
            if(end - start >= 3) { // 3アミノ酸以上なら螺旋として登録
                char buf[100];
                sprintf(buf, "HELIX  %3d %3d %3s A %4d  %3s A %4d  1\n", 
                        helixId, helixId, protein[start].resName.c_str(), protein[start].resSeq,
                        protein[end].resName.c_str(), protein[end].resSeq);
                out << buf;
                helixId++;
            }
        }
    }

    // フルアトム出力 (順番が重要: N -> CA -> C -> O)
    int atomId = 1;
    for (auto& r : protein) {
        char buf[100];
        sprintf(buf, "ATOM  %5d  N   %3s A%4d    %8.3f%8.3f%8.3f  1.00 20.00           N\n", atomId++, r.resName.c_str(), r.resSeq, r.n.x, r.n.y, r.n.z); out << buf;
        sprintf(buf, "ATOM  %5d  CA  %3s A%4d    %8.3f%8.3f%8.3f  1.00 20.00           C\n", atomId++, r.resName.c_str(), r.resSeq, r.ca.x, r.ca.y, r.ca.z); out << buf;
        sprintf(buf, "ATOM  %5d  C   %3s A%4d    %8.3f%8.3f%8.3f  1.00 20.00           C\n", atomId++, r.resName.c_str(), r.resSeq, r.c.x, r.c.y, r.c.z); out << buf;
        sprintf(buf, "ATOM  %5d  O   %3s A%4d    %8.3f%8.3f%8.3f  1.00 20.00           O\n", atomId++, r.resName.c_str(), r.resSeq, r.o.x, r.o.y, r.o.z); out << buf;
    }

    // オキシトシン環用 CONECT
    for(int i=0; i<(int)protein.size(); i++) {
        if (protein[i].resName == "CYS") {
            int bestPartner = -1;
            double bestDist = 8.0; 
            for(int j=i+1; j<(int)protein.size(); j++) {
                if (protein[j].resName == "CYS") {
                    double d = length(protein[i].ca - protein[j].ca);
                    if(d < bestDist) { bestDist = d; bestPartner = j; }
                }
            }
            if (bestPartner != -1) {
                // 原子IDは 4*(resSeq-1) + 2 がCA
                int atomI = 4*i + 2; int atomJ = 4*bestPartner + 2;
                out << "CONECT " << std::setw(5) << atomI << std::setw(5) << atomJ << "\n";
            }
        }
    }

    out << "END\n";
    std::cout << "[Complete] 'Universe_v13_Ribbon.pdb' generated." << std::endl;
    std::cout << "Mol* Viewer should automatically show beautiful ribbons now!" << std::endl;

    return 0;
}
