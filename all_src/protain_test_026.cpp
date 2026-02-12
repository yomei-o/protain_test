#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <cmath>
#include <iomanip>
#include <map>
#include <algorithm>

// --- Project-137: Universe OS v15.1 "Universal Parser" ---
// コンセプト: 入力配列が塩基配列(DNA/RNA)かアミノ酸配列かを自動判定する。

const double PI = 3.1415926535;
const double CA_DIST = 3.8;
const double MIN_DIST = 4.0;

struct Vec3 {
    double x, y, z;
    Vec3 operator-(const Vec3& o) const { return {x - o.x, y - o.y, z - o.z}; }
    Vec3 operator+(const Vec3& o) const { return {x + o.x, y + o.y, z + o.z}; }
    Vec3 operator*(double s) const { return {x * s, y * s, z * s}; }
};

Vec3 cross(Vec3 a, Vec3 b) { return {a.y*b.z - a.z*b.y, a.z*b.x - a.x*b.z, a.x*b.y - a.y*b.x}; }
double length(Vec3 v) { return std::sqrt(v.x*v.x + v.y*v.y + v.z*v.z); }
Vec3 normalize(Vec3 v) { double l = length(v); return (l>1e-6)? Vec3{v.x/l, v.y/l, v.z/l} : Vec3{0,0,0}; }

struct Residue {
    int resSeq; char code; std::string resName;
    double hydrophobicity; double charge;
    bool isHelix; 
    Vec3 ca, n, c, o;
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
    return {"UNK",0,0}; // 未知のアミノ酸
}

char codonToAA(const std::string& codon) {
     static std::map<std::string, char> m = {
        {"UUU",'F'}, {"UUC",'F'}, {"UUA",'L'}, {"UUG",'L'},
        {"CUU",'L'}, {"CUC",'L'}, {"CUA",'L'}, {"CUG",'L'},
        {"AUU",'I'}, {"AUC",'I'}, {"AUA",'I'}, {"AUG",'M'}, 
        {"GUU",'V'}, {"GUC",'V'}, {"GUA",'V'}, {"GUG",'V'},
        {"UCU",'S'}, {"UCC",'S'}, {"UCA",'S'}, {"UCG",'S'},
        {"CCU",'P'}, {"CCC",'P'}, {"CCA",'P'}, {"CCG",'P'},
        {"ACU",'T'}, {"ACC",'T'}, {"ACA",'T'}, {"ACG",'T'},
        {"GCU",'A'}, {"GCC",'A'}, {"GCA",'A'}, {"GCG",'A'},
        {"UAU",'Y'}, {"UAC",'Y'}, {"UAA",'_'}, {"UAG",'_'}, 
        {"CAU",'H'}, {"CAC",'H'}, {"CAA",'Q'}, {"CAG",'Q'},
        {"AAU",'N'}, {"AAC",'N'}, {"AAA",'K'}, {"AAG",'K'},
        {"GAU",'D'}, {"GAC",'D'}, {"GAA",'E'}, {"GAG",'E'},
        {"UGU",'C'}, {"UGC",'C'}, {"UGA",'_'}, {"UGG",'W'},
        {"CGU",'R'}, {"CGC",'R'}, {"CGA",'R'}, {"CGG",'R'},
        {"AGU",'S'}, {"AGC",'S'}, {"AGA",'R'}, {"AGG",'R'},
        {"GGU",'G'}, {"GGC",'G'}, {"GGA",'G'}, {"GGG",'G'}
    };
    if (m.count(codon)) return m[codon];
    return '?';
}

int main() {
    std::string rawInput;
    std::cout << "--- Project-137 v15.1 [Universal Parser] ---" << std::endl;
    std::cout << "Input Sequence (Amino Acid OR DNA/RNA): ";
    std::cin >> rawInput;
    
    // 入力文字を大文字に統一
    std::transform(rawInput.begin(), rawInput.end(), rawInput.begin(), ::toupper);

    // 【新機能】オートディテクト (DNA/RNA か アミノ酸 か)
    bool isNucleotide = true;
    for (char c : rawInput) {
        if (c != 'A' && c != 'T' && c != 'G' && c != 'C' && c != 'U') {
            isNucleotide = false; // ATGCU以外の文字（MやLなど）があればアミノ酸！
            break;
        }
    }

    std::string aaSeq = "";

    if (isNucleotide && rawInput.length() >= 3) {
        std::cout << "[System] Genetic code detected. Booting Ribosome..." << std::endl;
        // DNAをRNAに転写
        for(char& c : rawInput) { if(c == 'T') c = 'U'; }
        
        // RNAをアミノ酸に翻訳
        for(size_t i = 0; i + 2 < rawInput.length(); i += 3) {
            std::string codon = rawInput.substr(i, 3);
            char aa = codonToAA(codon);
            if(aa == '_' || aa == '?') {
                std::cout << "[Ribosome] Stop codon reached at index " << i << std::endl;
                break;
            }
            aaSeq += aa;
        }
    } else {
        std::cout << "[System] Amino Acid sequence detected. Bypassing Ribosome..." << std::endl;
        aaSeq = rawInput;
    }

    if (aaSeq.empty()) {
        std::cout << "Error: No valid sequence to process." << std::endl;
        return 1;
    }

    std::cout << "[Sequence] " << aaSeq << " (" << aaSeq.length() << " residues)" << std::endl;

    // --- 以降、v14と同等の完璧な物理エンジンとフルアトム生成 ---
    
    std::vector<Residue> protein;
    for(int i=0; i<(int)aaSeq.length(); ++i) {
        Residue res; res.resSeq = i+1; res.code = aaSeq[i];
        AAProps p = getAAProps(res.code);
        res.resName = p.name; res.hydrophobicity = p.hydro; res.charge = p.charge;
        res.isHelix = false;
        protein.push_back(res);
    }

    double currX=0, currY=0, currZ=0;
    double angleX=0.8, angleY=0.0;

    // 骨格計算
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
            protein[i].isHelix = true;
            angleX = 0.9;
            angleY += 1.74533; // 完璧な螺旋角度(100度)
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

    // 剛体ロック付きリラクゼーション
    for(int iter=0; iter<200; ++iter) {
        for(int i=0; i<(int)protein.size(); ++i) {
            for(int j=i+1; j<(int)protein.size(); ++j) {
                if(protein[i].isHelix && protein[j].isHelix && std::abs(i - j) < 6) continue;

                Vec3 d = protein[j].ca - protein[i].ca;
                double dist = length(d);
                if(dist < MIN_DIST && dist > 0.1) {
                    double push = (MIN_DIST - dist) * 0.5 * 0.1;
                    Vec3 move = normalize(d) * push;
                    if(!protein[i].isHelix) protein[i].ca = protein[i].ca - move;
                    if(!protein[j].isHelix) protein[j].ca = protein[j].ca + move;
                }
            }
        }
    }

    // フルアトム計算(ペプチド平面)
    for(int i=0; i<(int)protein.size(); ++i) {
        Vec3 ca_prev = (i > 0) ? protein[i-1].ca : protein[i].ca - Vec3{1,0,0};
        Vec3 ca_next = (i < (int)protein.size()-1) ? protein[i+1].ca : protein[i].ca + Vec3{1,0,0};

        Vec3 v_prev = normalize(protein[i].ca - ca_prev);
        Vec3 v_next = normalize(ca_next - protein[i].ca);

        protein[i].n = protein[i].ca - v_prev * 1.45;
        protein[i].c = protein[i].ca + v_next * 1.52;

        Vec3 ortho = cross(v_prev, v_next);
        if(length(ortho) < 0.1) ortho = {0,1,0};
        ortho = normalize(ortho);
        protein[i].o = protein[i].c + ortho * 1.23;
    }

    // PDB出力
    std::ofstream out("Universe_v15_1_Universal.pdb");
    out << "HEADER    PROJECT-137 V15.1 UNIVERSAL PARSER\n";

    int helixId = 1;
    for(int i=0; i<(int)protein.size(); i++) {
        if(protein[i].isHelix) {
            int start = i;
            while(i < (int)protein.size() && protein[i].isHelix) i++;
            int end = i - 1;
            if(end - start >= 4) {
                char buf[100];
                sprintf(buf, "HELIX  %3d %3d %3s A %4d  %3s A %4d  1\n", 
                        helixId, helixId, protein[start].resName.c_str(), protein[start].resSeq,
                        protein[end].resName.c_str(), protein[end].resSeq);
                out << buf;
                helixId++;
            }
        }
    }

    int atomId = 1;
    for (auto& r : protein) {
        char buf[100];
        sprintf(buf, "ATOM  %5d  N   %3s A%4d    %8.3f%8.3f%8.3f  1.00 20.00           N\n", atomId++, r.resName.c_str(), r.resSeq, r.n.x, r.n.y, r.n.z); out << buf;
        sprintf(buf, "ATOM  %5d  CA  %3s A%4d    %8.3f%8.3f%8.3f  1.00 20.00           C\n", atomId++, r.resName.c_str(), r.resSeq, r.ca.x, r.ca.y, r.ca.z); out << buf;
        sprintf(buf, "ATOM  %5d  C   %3s A%4d    %8.3f%8.3f%8.3f  1.00 20.00           C\n", atomId++, r.resName.c_str(), r.resSeq, r.c.x, r.c.y, r.c.z); out << buf;
        sprintf(buf, "ATOM  %5d  O   %3s A%4d    %8.3f%8.3f%8.3f  1.00 20.00           O\n", atomId++, r.resName.c_str(), r.resSeq, r.o.x, r.o.y, r.o.z); out << buf;
    }

    // Cys-Cysリング結合
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
                int atomI = 4*i + 2; int atomJ = 4*bestPartner + 2;
                out << "CONECT " << std::setw(5) << atomI << std::setw(5) << atomJ << "\n";
            }
        }
    }

    out << "END\n";
    std::cout << "[Complete] Output saved to 'Universe_v15_1_Universal.pdb'" << std::endl;

    return 0;
}

