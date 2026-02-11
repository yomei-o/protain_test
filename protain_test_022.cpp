#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <cmath>
#include <iomanip>
#include <map>
#include <algorithm>

// --- Project-137: Universe OS v11.0 "Collision & Relaxation" ---
// コンセプト: 「パーソナルスペース」の実装。
// 生成後に衝突判定を行い、めり込んだ原子同士を物理的に引き離す（斥力）。

const double PI = 3.1415926535;
const double CA_DIST = 3.8;      // 理想的な結合距離
const double MIN_DIST = 4.0;     // 原子同士の最小距離（これより近いと反発）

struct Residue {
    int resSeq;
    char code;
    std::string resName;
    double hydrophobicity;
    double charge;
    double x, y, z;
    bool hasSSBond; // SS結合済みフラグ
};

// アミノ酸データ
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

void normalize(double& x, double& y, double& z) {
    double len = std::sqrt(x*x + y*y + z*z);
    if(len > 1e-6) { x/=len; y/=len; z/=len; }
}

int main() {
    std::string inputSeq;
    std::cout << "Project-137 v11.0 [Collision Solver]" << std::endl;
    std::cout << "Input: ";
    std::cin >> inputSeq;

    std::vector<Residue> protein;
    for(int i=0; i<(int)inputSeq.length(); ++i) {
        Residue res;
        res.resSeq = i+1; res.code = inputSeq[i];
        AAProps p = getAAProps(res.code);
        res.resName = p.name; res.hydrophobicity = p.hydro; res.charge = p.charge;
        res.hasSSBond = false;
        protein.push_back(res);
    }

    // 1. 初期生成（v10の物理ロジック）
    double currX=0, currY=0, currZ=0;
    double angleX=0.6, angleY=0.0; // Helix初期値

    for(int i=0; i<(int)protein.size(); ++i) {
        // 物理パラメータ取得
        double avgHydro = 0; int count=0;
        for(int k=i-2; k<=i+2; k++) {
            if(k>=0 && k<(int)protein.size()) { avgHydro+=protein[k].hydrophobicity; count++; }
        }
        avgHydro /= (count>0?count:1);

        double localChargeSum = 0;
        for(int k=i-2; k<=i+2; k++) {
            if(k>=0 && k<(int)protein.size()) localChargeSum += std::abs(protein[k].charge);
        }

        // 挙動決定
        if(avgHydro > 1.0) { // 疎水性＝硬い
            double targetAngle = 0.6;
            angleX = angleX*0.8 + targetAngle*0.2;
            angleY += 1.0;
        } else { // 親水性＝曲がる
            double repulsion = 1.0 + (localChargeSum * 0.2);
            angleX += (PI*(3.0-std::sqrt(5.0))) * 0.5 * repulsion;
            angleY += 0.5 * repulsion;
        }

        double dirX = std::sin(angleX)*std::cos(angleY);
        double dirY = std::sin(angleX)*std::sin(angleY);
        double dirZ = std::cos(angleX);
        normalize(dirX,dirY,dirZ);
        currX += dirX*CA_DIST; currY += dirY*CA_DIST; currZ += dirZ*CA_DIST;
        
        protein[i].x = currX; protein[i].y = currY; protein[i].z = currZ;
    }

    // 2. 衝突解消ループ (Relaxation) - ここがv11の肝！
    // 「ぐちゃぐちゃ」を直すための整列運動
    std::cout << "Solving collisions..." << std::endl;
    int iterations = 500; // 500回揺する
    double pushStrength = 0.1; // 押し返す強さ

    for(int iter=0; iter<iterations; ++iter) {
        for(int i=0; i<(int)protein.size(); ++i) {
            for(int j=i+1; j<(int)protein.size(); ++j) {
                // 距離計算
                double dx = protein[j].x - protein[i].x;
                double dy = protein[j].y - protein[i].y;
                double dz = protein[j].z - protein[i].z;
                double dist = std::sqrt(dx*dx + dy*dy + dz*dz);

                // もし近すぎたら（衝突！）
                if(dist < MIN_DIST && dist > 0.1) {
                    // 反発ベクトル作成
                    double push = (MIN_DIST - dist) * 0.5 * pushStrength;
                    double px = (dx/dist) * push;
                    double py = (dy/dist) * push;
                    double pz = (dz/dist) * push;

                    // お互いに逆方向へ移動
                    protein[i].x -= px; protein[i].y -= py; protein[i].z -= pz;
                    protein[j].x += px; protein[j].y += py; protein[j].z += pz;
                }
            }
        }
        
        // 鎖が千切れないように結合距離(3.8)を強制補正
        for(int i=0; i<(int)protein.size()-1; ++i) {
             double dx = protein[i+1].x - protein[i].x;
             double dy = protein[i+1].y - protein[i].y;
             double dz = protein[i+1].z - protein[i].z;
             double dist = std::sqrt(dx*dx + dy*dy + dz*dz);
             if(dist > 0.1) {
                 double correction = (dist - CA_DIST) * 0.5; // 伸びすぎ・縮みすぎを直す
                 double cx = (dx/dist) * correction;
                 double cy = (dy/dist) * correction;
                 double cz = (dz/dist) * correction;
                 
                 protein[i].x += cx; protein[i].y += cy; protein[i].z += cz;
                 protein[i+1].x -= cx; protein[i+1].y -= cy; protein[i+1].z -= cz;
             }
        }
    }

    // 3. PDB出力 (SS結合の整理付き)
    std::ofstream out("Universe_v11_Collision.pdb");
    out << "HEADER    PROJECT-137 V11.0 COLLISION SOLVER\n";
    for (auto& r : protein) {
        char buf[100];
        sprintf(buf, "ATOM  %5d  CA  %3s A%4d    %8.3f%8.3f%8.3f  1.00  4.00           C\n",
                r.resSeq, r.resName.c_str(), r.resSeq, r.x, r.y, r.z);
        out << buf;
    }

    // スマートなSS結合判定（1対1対応）
    for(int i=0; i<(int)protein.size(); i++) {
        if (protein[i].resName == "CYS" && !protein[i].hasSSBond) {
            int bestPartner = -1;
            double bestDist = 10.0; // 探索範囲

            // 最も近い未結合のCYSを探す
            for(int j=i+1; j<(int)protein.size(); j++) {
                if (protein[j].resName == "CYS" && !protein[j].hasSSBond) {
                    double dx = protein[i].x - protein[j].x;
                    double dy = protein[i].y - protein[j].y;
                    double dz = protein[i].z - protein[j].z;
                    double dist = std::sqrt(dx*dx + dy*dy + dz*dz);
                    
                    if (dist < bestDist) {
                        bestDist = dist;
                        bestPartner = j;
                    }
                }
            }
            
            // パートナーが見つかったら結合
            if (bestPartner != -1) {
                out << "CONECT " << std::setw(5) << protein[i].resSeq << std::setw(5) << protein[bestPartner].resSeq << "\n";
                protein[i].hasSSBond = true;
                protein[bestPartner].hasSSBond = true;
            }
        }
    }
    out << "END\n";
    
    std::cout << "[Complete] Generated with collision detection." << std::endl;
    return 0;
}
