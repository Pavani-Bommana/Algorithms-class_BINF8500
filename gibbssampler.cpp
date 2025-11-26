#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <cmath>
#include <cstdlib>
#include <ctime>

#define MAXMOTLEN 100        
#define MINMOTLEN 6          
#define NSEEDS 10            
#define ADJCYCLES 25          
#define MAXADJ 3             
#define PLATEAUCYCLES 200    
#define MAXCYCLES 2000       
#define PSEUDO 0.25         

using namespace std;

struct Sequence {
    string header;
    string seq;
    int len;
    int motPos; 
};

struct Result {
    double score;
    int motLen;          
    vector<int> motPos;  
};


vector<Sequence> read_fasta(const string& filename) {
    vector<Sequence> seqs;
    ifstream file(filename.c_str());
    if (!file) {
        cerr << "Error: cannot open " << filename << "\n";
        return seqs;
    }
    string line, header, seq;
    while (getline(file, line)) {
        if (line.empty()) continue;
        if (line[0] == '>') {
            if (!seq.empty()) {
                seqs.push_back({header, seq, (int)seq.size(), 0});
                seq.clear();
            }
            header = line.substr(1);
        } else {
            for (size_t i=0;i<line.size();++i) {
                char c = line[i];
                if (!isspace((unsigned char)c)) {
                    seq.push_back((char)toupper(c));
                }
            }
        }
    }
    if (!seq.empty()) seqs.push_back({header, seq, (int)seq.size(), 0});
    return seqs;
}

void estimate_background(const vector<Sequence>& seqs, double bg[4]) {
    for (int k=0;k<4;++k) bg[k]=PSEUDO;
    for (size_t i=0;i<seqs.size();++i) {
        for (int j=0;j<seqs[i].len;++j) {
            char b = seqs[i].seq[j];
            int idx = (b=='A'?0:b=='C'?1:b=='G'?2:b=='T'?3:-1);
            if (idx>=0) bg[idx]+=1.0;
        }
    }
    double sum = 0.0;
    for (int k=0;k<4;++k) sum += bg[k];
    for (int k=0;k<4;++k) bg[k] /= sum;
}


void build_pwm_leave_one(const vector<Sequence>& seqs, const vector<int>& pos, int motLen, int skip, double pwm[MAXMOTLEN][4]) {
    for (int j=0;j<motLen;++j) for (int k=0;k<4;++k) pwm[j][k]=PSEUDO;
    for (size_t i=0;i<seqs.size();++i) {
        if ((int)i==skip) continue;
        int p = pos[i];
        // ensure valid window
        if (p<0 || p+motLen>seqs[i].len) continue;
        for (int j=0;j<motLen;++j) {
            char b = seqs[i].seq[p+j];
            int idx = (b=='A'?0:b=='C'?1:b=='G'?2:b=='T'?3:-1);
            if (idx>=0) pwm[j][idx] += 1.0;
        }
    }
    
    for (int j=0;j<motLen;++j) {
        double s=0.0;
        for (int k=0;k<4;++k) s+=pwm[j][k];
        if (s==0.0) s=1.0;
        for (int k=0;k<4;++k) pwm[j][k] /= s;
    }
}


double score_kmer_log_odds(const string& s, int start, int motLen, double pwm[MAXMOTLEN][4], const double bg[4]) {
    double sc = 0.0;
    for (int j=0;j<motLen;++j) {
        char b = s[start+j];
        int idx = (b=='A'?0:b=='C'?1:b=='G'?2:b=='T'?3:-1);
        if (idx<0) return -INFINITY; // if non-ACGT, skip
        double p = pwm[j][idx];
        double q = bg[idx];
        if (p<=0.0 || q<=0.0) return -INFINITY;
        sc += log(p/q)/log(2.0); // log2
    }
    return sc;
}


int sample_new_position(const Sequence& seq, int motLen, double pwm[MAXMOTLEN][4], const double bg[4]) {
    int L = seq.len;
    int W = motLen;
    int nwin = L - W + 1;
    if (nwin <= 0) return 0;
    vector<double> cum(nwin, 0.0);
    
    for (int i=0;i<nwin;++i) {
        double s = score_kmer_log_odds(seq.seq, i, W, pwm, bg);
        if (!isfinite(s)) s = -1e9;
        double r = pow(2.0, s); 
        cum[i] = r;
    }
    
    for (int i=1;i<nwin;++i) cum[i] += cum[i-1];
    double total = cum[nwin-1];
    if (total<=0.0) {
        
        double u = (double)rand()/RAND_MAX;
        int pos = (int)(u * nwin);
        if (pos<0) pos=0;
        if (pos>=nwin) pos=nwin-1;
        return pos;
    }
    double R = ((double)rand()/RAND_MAX) * total;
    int pos = 0;
    while (pos<nwin-1 && R >= cum[pos]) ++pos;
    return pos;
}


double total_log_odds_score(const vector<Sequence>& seqs, const vector<int>& pos, int motLen, const double bg[4]) {
    double total = 0.0;
    double pwm[MAXMOTLEN][4];
    for (size_t i=0;i<seqs.size();++i) {
        build_pwm_leave_one(seqs, pos, motLen, (int)i, pwm);
        int p = pos[i];
        if (p<0 || p+motLen>seqs[i].len) continue;
        double sc = score_kmer_log_odds(seqs[i].seq, p, motLen, pwm, bg);
        if (isfinite(sc)) total += sc;
    }
    return total;
}


void gibbs_cycle(vector<Sequence>& seqs, vector<int>& pos, int motLen, const double bg[4]) {
    double pwm[MAXMOTLEN][4];
    for (size_t i=0;i<seqs.size();++i) {
        build_pwm_leave_one(seqs, pos, motLen, (int)i, pwm);
        int newp = sample_new_position(seqs[i], motLen, pwm, bg);
        pos[i] = newp;
    }
}


int adjust_motif_length(vector<Sequence>& seqs, vector<int>& pos, int motLen, const double bg[4]) {
    int bestLen = motLen;
    double bestScore = total_log_odds_score(seqs, pos, motLen, bg);


    for (int delta=-MAXADJ; delta<=MAXADJ; ++delta) {
        if (delta==0) continue;
        int newLen = motLen + delta;
        if (newLen < MINMOTLEN || newLen > MAXMOTLEN) continue;


        vector<int> newPos = pos;
        for (size_t i=0;i<seqs.size();++i) {
            int maxStart = seqs[i].len - newLen;
            if (maxStart < 0) { newPos[i] = 0; continue; }
            if (newPos[i] > maxStart) newPos[i] = maxStart;
            if (newPos[i] < 0) newPos[i] = 0;
        }


        gibbs_cycle(seqs, newPos, newLen, bg);

        double sc = total_log_odds_score(seqs, newPos, newLen, bg);
        if (sc > bestScore) {
            bestScore = sc;
            bestLen = newLen;
            pos = newPos;
        }
    }
    return bestLen;
}


Result run_single_seed(vector<Sequence>& seqs, int motLen0, const double bg[4]) {
    int N = (int)seqs.size();
    vector<int> pos(N, 0);

    for (int i=0;i<N;++i) {
        int maxStart = seqs[i].len - motLen0;
        if (maxStart < 0) maxStart = 0;
        pos[i] = (maxStart>0) ? (rand() % (maxStart+1)) : 0;
    }

    int motLen = motLen0;
    double bestScore = -1e9;
    vector<int> bestPos = pos;
    int bestLen = motLen;

    int NCyc = 0;
    int NPlat = 0;
    double seedBest = -1e9;

    while (NCyc < MAXCYCLES && NPlat < PLATEAUCYCLES) {
        
        gibbs_cycle(seqs, pos, motLen, bg);

        
        double sc = total_log_odds_score(seqs, pos, motLen, bg);
        if (sc > seedBest) { seedBest = sc; NPlat = 0; } else { NPlat++; }

        if (sc > bestScore) {
            bestScore = sc;
            bestPos = pos;
            bestLen = motLen;
        }


        if (NCyc > 0 && (NCyc % ADJCYCLES == 0)) {
            int newLen = adjust_motif_length(seqs, pos, motLen, bg);
            if (newLen != motLen) {
                motLen = newLen; 
                double sc2 = total_log_odds_score(seqs, pos, motLen, bg);
                if (sc2 > bestScore) {
                    bestScore = sc2;
                    bestPos = pos;
                    bestLen = motLen;
                }
                NPlat = 0; 
            }
        }

        NCyc++;
    }

    Result r;
    r.score = bestScore;
    r.motLen = bestLen;
    r.motPos = bestPos;
    return r;
}

Result run_multi_seed(vector<Sequence>& seqs, int motLen0, const double bg[4]) {
    Result best;
    best.score = -1e9;
    best.motLen = motLen0;
    best.motPos.assign(seqs.size(), 0);

    for (int s=0; s<NSEEDS; ++s) {
        Result r = run_single_seed(seqs, motLen0, bg);
        if (r.score > best.score) best = r;
    }
    return best;
}

int main(int argc, char* argv[]) {
    if (argc < 3) {
        cerr << "Usage: ./gibbs_sampler <fasta_file> <estimated_motif_length>\n";
        return 1;
    }
    srand((unsigned)time(NULL));

    string filename = argv[1];
    int motLen0 = atoi(argv[2]);
    if (motLen0 < MINMOTLEN || motLen0 > MAXMOTLEN) {
        cerr << "Estimated motif length must be between " << MINMOTLEN << " and " << MAXMOTLEN << "\n";
        return 1;
    }

    vector<Sequence> seqs = read_fasta(filename);
    if (seqs.empty()) {
        cerr << "No sequences read. Check FASTA input.\n";
        return 1;
    }

    double bg[4];
    estimate_background(seqs, bg);

    Result best = run_multi_seed(seqs, motLen0, bg);

    cout << "\nSimple Gibbs motif sampler output:\n\n";
    cout << "     Input file: " << filename << "\n";
    cout << "     Initial motif length: " << motLen0 << "\n";
    cout << "     Final motif length:   " << best.motLen << "\n";
    cout << "     Final score:          " << best.score << "\n\n";
    cout << "Motif sequences and locations:\n\n";
    for (size_t i=0;i<seqs.size();++i) {
        int start = best.motPos[i] + 1; // 1-based
        int end = best.motPos[i] + best.motLen;
        string motif = seqs[i].seq.substr(best.motPos[i], best.motLen);
        // pad motif to align (optional)
        cout << motif;
        int pad = 22 - (int)motif.size();
        if (pad < 1) pad = 1;
        for (int k=0;k<pad;++k) cout << ' ';
        cout << start << "-" << end << " " << seqs[i].header << "\n";
    }

    return 0;
}