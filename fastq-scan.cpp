#include <unistd.h>
#include <iostream>
#include <iomanip>
#include <string>
#include <fstream>
#include <algorithm>
#include <vector>
#include <math.h>
#include <sstream>
using namespace std;
const string VERSION = "1.0.1";
const int MAX_READ_LENGTH = 50000000;

class Stats {
    public:
        unsigned long int total_bp;

        // Read stats
        vector<unsigned int> read_length;
        vector<unsigned int> read_length_count;
        unsigned long int read_total;
        unsigned int read_min;
        unsigned int read_max;
        double read_mean;
        double read_std;
        double read_median;
        double read_25th;
        double read_75th;

        // Qual stats
        vector<double> per_read_qual;
        vector<unsigned long int> per_base_qual;
        vector<unsigned long int> per_base_count;
        int phred;
        unsigned int qual_min;
        unsigned int qual_max;
        double qual_sum;
        double qual_mean;
        double qual_median;
        double qual_25th;
        double qual_75th;
        double qual_std;

        void init(void) {
            read_total = 0;
            total_bp = 0;
            qual_sum = 0;
            read_length_count.resize(MAX_READ_LENGTH,0);
            per_base_qual.resize(MAX_READ_LENGTH,0);
            per_base_count.resize(MAX_READ_LENGTH,0);
            phred = 33;
        }

        double get_percentile(vector<double> array, size_t size, float percentile) {
            if (size % 2 == 0) {
                int p = int(size * percentile - 1);
                return (array[p] + array[p]) / 2.0;
            } else {
                int p = int(size * percentile);
                return array[p];
            }
        }

        double get_std(vector<double> array, double mean) {
            double temp = 0;
            for (unsigned int i = 0; i < read_total; i++) {
                temp += pow(((array[i]-phred) - mean), 2);
            }
            return sqrt(temp / read_total);
        }

        double get_percentile(vector<unsigned int> array, size_t size, float percentile) {
            if (size % 2 == 0) {
                int p = int(size * percentile - 1);
                return (array[p] + array[p]) / 2.0;
            } else {
                int p = int(size * percentile);
                return array[p];
            }
        }

        double get_std(vector<unsigned int> array, double mean) {
            double temp = 0;
            for (unsigned int i = 0; i < read_total; i++) {
                temp += pow((array[i] - mean), 2);
            }
            return sqrt(temp / read_total);
        }

        int transform_quality(string qual) {
            unsigned int total = 0;
            total_bp += qual.length();
            read_length_count[qual.length()]++;
            for (unsigned int i = 0; i < qual.length(); i++) {
                unsigned int qual_val = (unsigned int)qual[i];
                per_base_qual[i] += qual_val;
                per_base_count[i]++;
                total += qual_val;
            }

            if (qual.length() > 0) {
                double avg_qual = total / qual.length();
                qual_sum += avg_qual;
                per_read_qual.push_back(avg_qual);
                return 0;
            } else {
                return 1;
            }
        }

        void read_stats(void) {
            sort(read_length.begin(), read_length.end());
            read_min = read_length.front();
            read_mean = total_bp / float(read_total);
            read_std = get_std(read_length, read_mean);
            read_max = read_length.back();
            read_25th = get_percentile(read_length, read_length.size(), 0.25);
            read_median = get_percentile(read_length, read_length.size(), 0.50);
            read_75th = get_percentile(read_length, read_length.size(), 0.75);
        }

        void qual_stats(void) {
            sort(per_read_qual.begin(), per_read_qual.end());
            qual_min = per_read_qual.front() - phred;
            qual_mean = (qual_sum / read_total) - phred;
            qual_std = get_std(per_read_qual, qual_mean);
            qual_max = per_read_qual.back() - phred;
            qual_25th = get_percentile(per_read_qual, per_read_qual.size(), 0.25) - phred;
            qual_median = get_percentile(per_read_qual, per_read_qual.size(), 0.50) - phred;
            qual_75th = get_percentile(per_read_qual, per_read_qual.size(), 0.75) - phred;
        }

        void jsonify_stats(float GENOME_SIZE, bool QC_ONLY) {
            string t1 = "    ";
            string t2 = "        ";
            cout << "{" << endl;
            cout << t1 << "\"qc_stats\": {" << endl;
            cout << t2 << "\"total_bp\": " << total_bp << "," << endl;
            if (GENOME_SIZE <= 1)  {
                cout << t2 << "\"coverage\": 0.00," << endl;
            } else {
                cout << t2 << "\"coverage\": " << total_bp / GENOME_SIZE << "," << endl;
            }
            cout << t2 << "\"read_total\": " << read_total << "," << endl;
            cout << t2 << "\"read_min\": " << read_min << "," << endl;
            cout << t2 << "\"read_mean\": " << read_mean << "," << endl;
            cout << t2 << "\"read_std\": " << read_std << "," << endl;
            cout << t2 << "\"read_median\": " << read_median << "," << endl;
            cout << t2 << "\"read_max\": " << read_max << "," << endl;
            cout << t2 << "\"read_25th\": " << read_25th << "," << endl;
            cout << t2 << "\"read_75th\": " << read_75th << "," << endl;
            cout << t2 << "\"qual_min\": " << qual_min << "," << endl;
            cout << t2 << "\"qual_mean\": " << qual_mean << "," << endl;
            cout << t2 << "\"qual_std\": " << qual_std << "," << endl;
            cout << t2 << "\"qual_max\": " << qual_max << "," << endl;
            cout << t2 << "\"qual_median\": " << qual_median << "," << endl;
            cout << t2 << "\"qual_25th\": " << qual_25th << "," << endl;
            cout << t2 << "\"qual_75th\": " << qual_75th << endl;
            
            
            if (QC_ONLY) {
                cout << t1 << "}" << endl;
            } else {
                cout << t1 << "}," << endl;
                cout << t1 << "\"read_lengths\": {" << endl;
                for (unsigned int i = read_min; i <= read_max; i++) {
                    if (i % 5 == 0) {
                        cout << endl;
                    }
                    cout << t2 << "\"" << i << "\": " << read_length_count[i];
                    if (i < read_max) {
                        cout << ",";
                    }
                }
                cout << endl << t1 << "}," << endl;
                cout << t1 << "\"per_base_quality\": {" << endl;
                for (unsigned int i = 0; i < read_max; i++) {
                    if (i % 5 == 0 && i != 0) {
                        cout << endl;
                    }
                    cout << t2 << "\"" << i + 1 << "\": " << (per_base_qual[i] / float(per_base_count[i])) - phred;
                    if (i < read_max - 1) {
                        cout << ",";
                    }
                }
                cout << endl << t1 << "}" << endl;
            }
            cout << "}" << endl;
        }
};

static int usage()
{
    cout << "Usage: cat FASTQ | fastq-scan [options]" << endl;
    cout << "Version: " << VERSION << endl;
    cout << endl;
    cout << "Optional arguments:" << endl;
    cout << "    -g INT   Genome size for calculating estimated sequencing coverage. (Default 1)" << endl;
    cout << "    -p INT   ASCII offset for input quality scores, can be 33 or 64. (Default 33)" << endl;
    cout << "    -q       Print only the QC stats, do not print read lengths or per-base quality scores" << endl;
    cout << "    -v       Print version information and exit" << endl;
    cout << "    -h       Show this message and exit" << endl;
    cout << endl;
    return 0;
}

static int version()
{
    cout << "fastq-scan " << VERSION << endl;
    return 0;
}

int main(int argc, char **argv) {
    // Read command line
    float GENOME_SIZE =  1.0;
    int PHRED_OFFSET = 33;
    bool QC_ONLY = false;
    int opt;
    while ((opt = getopt(argc, argv, "g:p:qvh")) >= 0) {
        switch (opt) {
            case 'g': GENOME_SIZE = atof(optarg); break;
            case 'p': PHRED_OFFSET = atoi(optarg); break;
            case 'q': QC_ONLY = true; break;
            case 'v': return version();
            case 'h': return usage();
        }
    }
    if (!(PHRED_OFFSET == 33 || PHRED_OFFSET == 64)) {
        cerr << "Invalid value for -p (" << PHRED_OFFSET << "), only 33 or 64 are valid" << endl;
        return 1;
    } else if (GENOME_SIZE < 0) {
        cerr << "Invalid value for -g (" << GENOME_SIZE << "), value muse be >= 0" << endl;
        return 1;
    }
    
    if (isatty(0)) return usage();

    // Parse FASTQ
    Stats stats;
    stats.init();
    stats.phred = PHRED_OFFSET;
    string name, seq, plus, qual;
    ifstream in("/dev/stdin", ios::in);
    while(true) {
        if(!getline(in, name, '\n')) break;
        if(!getline(in, seq, '\n')) break;
        if(!getline(in, plus, '\n')) break;
        if(!getline(in, qual, '\n')) break;
        stats.read_length.push_back(seq.length());
        stats.read_total++;
        int missing_qual = stats.transform_quality(qual);
        if (missing_qual == 1){
            cerr << "WARNING: Missing quality for a read.\n";
            cerr << name << "\n";
            cerr << seq << "\n";
            cerr << plus << "\n";
            cerr << qual << "\n";
            cerr << "Please fix it to continue.\n";
            exit(50);
        }
    }
    in.close();

    if (stats.read_total > 0) {
        // Determine Stats
        stats.read_stats();
        stats.qual_stats();
        stats.jsonify_stats(GENOME_SIZE, QC_ONLY);
    } else {
        // Empty file, or nothing to parse
        cerr << "Nothing to parse in STDIN, please verify your file is not empty." << endl;
    }
    return 0;
}
