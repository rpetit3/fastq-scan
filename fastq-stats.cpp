#include <iostream>
#include <iomanip>
#include <string>
#include <fstream>
#include <algorithm>
#include <vector>
#include <math.h>
#include <sstream>
using namespace std;
const float VERSION = 0.1;
const int MAX_READ_LENGTH = 50000;

class Stats {
    public:
        unsigned long int total_bp;

        // Read stats
        vector<unsigned int> read_length;
        vector<unsigned int> read_length_count;
        unsigned long int read_total;
        unsigned int read_min;
        double read_mean;
        double read_std;
        double read_median;
        unsigned int read_max;
        double read_25th;
        double read_75th;

        // Qual stats
        vector<double> per_read_qual;
        vector<unsigned int> per_base_qual;
        vector<unsigned int> per_base_count;
        int phred;
        unsigned int min_phred;
        unsigned int max_phred;
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
            min_phred = 1000;
            max_phred = -1000;
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

        void transform_quality(string qual) {
            unsigned int total = 0;
            total_bp += qual.length();
            read_length_count[qual.length()]++;
            for (unsigned int i = 0; i < qual.length(); i++) {
                unsigned int qual_val = (unsigned int)qual[i];
                if (qual_val < min_phred) {
                    min_phred = qual_val;
                } else if (qual_val > max_phred) {
                    max_phred = qual_val;
                }
                per_base_qual[i] += qual_val;
                per_base_count[i]++;
                total += qual_val;
            }
            double avg_qual = total / qual.length();
            qual_sum += avg_qual;
            per_read_qual.push_back(avg_qual);
        }

        void guess_phred() {
            if (max_phred > 74 && min_phred > 58) {
                phred = 64;
            } else if (max_phred <= 74 && min_phred >= 33) {
                phred = 33;
            } else {
                phred = 33;
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
            qual_mean = (qual_sum / read_total) - phred;
            qual_std = get_std(per_read_qual, qual_mean);
            qual_25th = get_percentile(per_read_qual, per_read_qual.size(), 0.25) - phred;
            qual_median = get_percentile(per_read_qual, per_read_qual.size(), 0.50) - phred;
            qual_75th = get_percentile(per_read_qual, per_read_qual.size(), 0.75) - phred;
        }

        void jsonify_stats(float GENOME_SIZE) {
            string t1 = "    ";
            string t2 = "        ";
            cout << "{" << endl;
            cout << t1 << "\"qc_stats\": {" << endl;
            cout << t2 << "\"total_bp\":" << total_bp << "," << endl;
            cout << t2 << "\"coverage\":" << total_bp / GENOME_SIZE << "," << endl;
            cout << t2 << "\"read_total\":" << read_total << "," << endl;
            cout << t2 << "\"read_min\":" << read_min << "," << endl;
            cout << t2 << "\"read_mean\":" << read_mean << "," << endl;
            cout << t2 << "\"read_std\":" << read_std << "," << endl;
            cout << t2 << "\"read_median\":" << read_median << "," << endl;
            cout << t2 << "\"read_max\":" << read_max << "," << endl;
            cout << t2 << "\"read_25th\":" << read_25th << "," << endl;
            cout << t2 << "\"read_75th\":" << read_75th << "," << endl;
            cout << t2 << "\"qual_mean\":" << qual_mean << "," << endl;
            cout << t2 << "\"qual_std\":" << qual_std << "," << endl;
            cout << t2 << "\"qual_median\":" << qual_median << "," << endl;
            cout << t2 << "\"qual_25th\":" << qual_25th << "," << endl;
            cout << t2 << "\"qual_75th\":" << qual_75th << endl;
            cout << t1 << "}," << endl;
            cout << t1 << "\"read_lengths\": {" << endl;
            for (unsigned int i = read_min; i <= read_max; i++) {
                if (i % 5 == 0) {
                    cout << endl;
                }
                cout << t2 << "\"" << i << "\":" << read_length_count[i];
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
                cout << t2 << "\"" << i + 1 << "\":" << (per_base_qual[i] / float(per_base_count[i])) - phred;
                if (i < read_max - 1) {
                    cout << ",";
                }
            }
            cout << endl << t1 << "}" << endl;
            cout << "}" << endl;
        }
};

int main(int argc,char **argv) {
    Stats stats;
    stats.init();
    string name, seq, plus, qual;
    ifstream in("/dev/stdin", ios::in);
    while(true) {
        if(!getline(in,name,'\n')) break;
        if(!getline(in,seq,'\n')) break;
        if(!getline(in,plus,'\n')) break;
        if(!getline(in,qual,'\n')) break;
        stats.read_length.push_back(seq.length());
        stats.read_total++;
        stats.transform_quality(qual);
    }
    in.close();

    // Determine Stats
    float GENOME_SIZE = argc == 1 ? 1.0 : atof(argv[1]);
    stats.guess_phred();
    stats.read_stats();
    stats.qual_stats();
    stats.jsonify_stats(GENOME_SIZE);
    return 0;
}
