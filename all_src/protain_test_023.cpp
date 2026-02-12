#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <cmath>
#include <iomanip>
#include <map>

// --- Project-137: Universe OS v12.0 "Visualizer Edition" ---
// 物理計算はv11と同じだが、出力フォーマットを「Viewerフレンドリー」に完全改良。
// 勝手な結合線を防ぎ、一本の鎖として美しく表示させる。

const double PI = 3.1415926535;
const double CA_DIST = 3.8;
const double MIN_DIST = 4.0;

struct Residue {
    int resSeq;
    char code;
    std::string resName;
    double hydrophobicity;
    double charge;
    double x, y, z;
    double stress; // 物理的ストレス（表示色用）
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

void normalize(double& x, double& y, double& z) {
    double len = std::sqrt(x*x + y*y + z*z);
    if(len > 1e-6) { x/=len; y/=len; z/=len; }
}

int main() {
    std::string inputSeq;
    std::cout << "Project-137 v12.0 [Visualizer Edition]" << std::endl;
    std::cout << "Input: ";
    std::cin >> inputSeq;

    std::vector<Residue> protein;
    for(int i=0; i<(int)inputSeq.length(); ++i) {
        Residue res;
        res.resSeq = i+1; res.code = inputSeq[i];
        AAProps p = getAAProps(res.code);
        res.resName = p.name; res.hydrophobicity = p.hydro; res.charge = p.charge;
        res.stress = 0.0;
        protein.push_back(res);
    }

    // 1. 初期生成 (v10 Logic)
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
        } else {
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

    // 2. 衝突解消 (v11 Logic)
    int iterations = 500;
    for(int iter=0; iter<iterations; ++iter) {
        for(int i=0; i<(int)protein.size(); ++i) {
            for(int j=i+1; j<(int)protein.size(); ++j) {
                double dx = protein[j].x - protein[i].x;
                double dy = protein[j].y - protein[i].y;
                double dz = protein[j].z - protein[i].z;
                double dist = std::sqrt(dx*dx + dy*dy + dz*dz);

                if(dist < MIN_DIST && dist > 0.1) {
                    double push = (MIN_DIST - dist) * 0.5 * 0.1;
                    protein[i].x -= (dx/dist)*push; protein[i].y -= (dy/dist)*push; protein[i].z -= (dz/dist)*push;
                    protein[j].x += (dx/dist)*push; protein[j].y += (dy/dist)*push; protein[j].z += (dz/dist)*push;
                    
                    // ストレス値を記録（赤く表示するため）
                    protein[i].stress += 0.1;
                    protein[j].stress += 0.1;
                }
            }
        }
        for(int i=0; i<(int)protein.size()-1; ++i) {
             double dx = protein[i+1].x - protein[i].x;
             double dy = protein[i+1].y - protein[i].y;
             double dz = protein[i+1].z - protein[i].z;
             double dist = std::sqrt(dx*dx + dy*dy + dz*dz);
             if(dist > 0.1) {
                 double correction = (dist - CA_DIST) * 0.5;
                 protein[i].x += (dx/dist)*correction; protein[i].y += (dy/dist)*correction; protein[i].z += (dz/dist)*correction;
                 protein[i+1].x -= (dx/dist)*correction; protein[i+1].y -= (dy/dist)*correction; protein[i+1].z -= (dz/dist)*correction;
             }
        }
    }

    // 3. PDB出力 (Visualization Fix)
    std::ofstream out("Universe_v12_Visual.pdb");
    out << "HEADER    PROJECT-137 V12.0 VISUALIZER\n";
    
    // ATOMレコード
    for (auto& r : protein) {
        // B-factor (tempFactor) にストレス値を入れる -> ビューワーで色分け可能に
        double tempFactor = (r.stress > 10.0) ? 99.0 : r.stress * 5.0; 
        
        char buf[100];
        // フォーマット厳守: %8.3f%8.3f%8.3f%6.2f%6.2f
        sprintf(buf, "ATOM  %5d  CA  %3s A%4d    %8.3f%8.3f%8.3f  1.00%6.2f           C\n",
                r.resSeq, r.resName.c_str(), r.resSeq, r.x, r.y, r.z, tempFactor);
        out << buf;
    }

    // CONECTレコード (Backbone)
    // これを入れるとViewerは「ああ、これは1本の鎖なんだ」と理解して綺麗に繋ぐ
    for(int i=0; i<(int)protein.size()-1; i++) {
        out << "CONECT " << std::setw(5) << protein[i].resSeq << std::setw(5) << protein[i+1].resSeq << "\n";
    }

    // CONECTレコード (SS-Bond)
    // Cys同士の結合だけを追加で記述
    for(int i=0; i<(int)protein.size(); i++) {
        if (protein[i].resName == "CYS") {
            int bestPartner = -1;
            double bestDist = 8.0; 
            for(int j=i+1; j<(int)protein.size(); j++) {
                if (protein[j].resName == "CYS") {
                    double d = std::sqrt(std::pow(protein[i].x-protein[j].x,2) + std::pow(protein[i].y-protein[j].y,2) + std::pow(protein[i].z-protein[j].z,2));
                    if(d < bestDist) { bestDist = d; bestPartner = j; }
                }
            }
            if (bestPartner != -1) {
                out << "CONECT " << std::setw(5) << protein[i].resSeq << std::setw(5) << protein[bestPartner].resSeq << "\n";
            }
        }
    }

    out << "END\n";
    
    std::cout << "[Complete] 'Universe_v12_Visual.pdb' generated." << std::endl;
    std::cout << "Open in Mol* -> It should look like a clean chain now." << std::endl;

    return 0;
}

