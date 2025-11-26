#include <iostream>
#include <fstream>
#include <vector>
#include <sstream>
#include <cmath>
#include <cstdlib>
#include <ctime>
#include <limits>
#include <iomanip>
#include <string>

using namespace std;

double euclidean(const vector<double>& a, const vector<double>& b) {
    double sum = 0.0;
    for (size_t i = 0; i < a.size(); ++i)
        sum += pow(a[i] - b[i], 2);
    return sqrt(sum);
}

vector<vector<double>> initialize_centroids(const vector<vector<double>>& data, int k) {
    vector<vector<double>> centroids;
    for (int i = 0; i < k; ++i)
        centroids.push_back(data[rand() % data.size()]);
    return centroids;
}

vector<int> assign_clusters(const vector<vector<double>>& data, const vector<vector<double>>& centroids) {
    vector<int> assignments(data.size());
    for (size_t i = 0; i < data.size(); ++i) {
        double min_dist = numeric_limits<double>::max();
        int cluster = 0;
        for (size_t j = 0; j < centroids.size(); ++j) {
            double dist = euclidean(data[i], centroids[j]);
            if (dist < min_dist) {
                min_dist = dist;
                cluster = j;
            }
        }
        assignments[i] = cluster;
    }
    return assignments;
}

vector<vector<double>> update_centroids(const vector<vector<double>>& data, const vector<int>& assignments, int k) {
    vector<vector<double>> centroids(k, vector<double>(data[0].size(), 0.0));
    vector<int> counts(k, 0);

    for (size_t i = 0; i < data.size(); ++i) {
        int cluster = assignments[i];
        for (size_t j = 0; j < data[i].size(); ++j)
            centroids[cluster][j] += data[i][j];
        counts[cluster]++;
    }

    for (int i = 0; i < k; ++i)
        for (size_t j = 0; j < centroids[i].size(); ++j)
            centroids[i][j] /= max(1, counts[i]);

    return centroids;
}

double compute_wcss(const vector<vector<double>>& data, const vector<vector<double>>& centroids, const vector<int>& assignments) {
    double wcss = 0.0;
    for (size_t i = 0; i < data.size(); ++i)
        wcss += pow(euclidean(data[i], centroids[assignments[i]]), 2);
    return wcss;
}

double compute_mean_distance(const vector<vector<double>>& data, const vector<vector<double>>& centroids, const vector<int>& assignments) {
    double total = 0.0;
    for (size_t i = 0; i < data.size(); ++i)
        total += euclidean(data[i], centroids[assignments[i]]);
    return total / data.size();
}

int main(int argc, char* argv[]) {
    if (argc != 4) {
        cout << "Usage: " << argv[0] << " <input_file> <k_min> <k_max>\n";
        return 1;
    }

    string filename = argv[1];
    int k_min = stoi(argv[2]);
    int k_max = stoi(argv[3]);

    ifstream file(filename);
    if (!file.is_open()) {
        cout << "Error opening file: " << filename << "\n";
        return 1;
    }

    string line, token;
    getline(file, line); 

    vector<vector<double>> data;
    vector<string> genome_labels;

    while (getline(file, line)) {
        istringstream ss(line);
        getline(ss, token, '\t');
        genome_labels.push_back(token);
        vector<double> row;
        while (getline(ss, token, '\t'))
            row.push_back(stod(token));
        data.push_back(row);
    }

    int n = data.size();
    int m = data[0].size();


    vector<double> means(m, 0.0), stddevs(m, 0.0);
    for (int j = 0; j < m; ++j)
        for (int i = 0; i < n; ++i)
            means[j] += data[i][j];
    for (int j = 0; j < m; ++j)
        means[j] /= n;
    for (int j = 0; j < m; ++j)
        for (int i = 0; i < n; ++i)
            stddevs[j] += pow(data[i][j] - means[j], 2);
    for (int j = 0; j < m; ++j)
        stddevs[j] = sqrt(stddevs[j] / n);
    for (int i = 0; i < n; ++i)
        for (int j = 0; j < m; ++j)
            data[i][j] = (data[i][j] - means[j]) / stddevs[j];

    cout << "K\tMean Distance\tWCSS\t\tAIC\t\tBIC\n";
    srand(time(0)); // seed once

    for (int k = k_min; k <= k_max; ++k) {
        int num_seeds = static_cast<int>(300 * k * sqrt(k));

        double best_wcss = numeric_limits<double>::max();
        vector<int> best_assignments;
        vector<vector<double>> best_centroids;

        for (int seed = 0; seed < num_seeds; ++seed) {
            vector<vector<double>> centroids = initialize_centroids(data, k);
            vector<int> assignments;

            for (int iter = 0; iter < 500; ++iter) {
                assignments = assign_clusters(data, centroids);
                centroids = update_centroids(data, assignments, k);
            }

            double wcss = compute_wcss(data, centroids, assignments);
            if (wcss < best_wcss) {
                best_wcss = wcss;
                best_assignments = assignments;
                best_centroids = centroids;
            }
        }

        double mean_dist = compute_mean_distance(data, best_centroids, best_assignments);
        double aic = best_wcss + 2 * m * k;
        double bic = best_wcss + log(n) * m * k;

        cout << fixed << setprecision(5);
        cout << k << "\t" << mean_dist << "\t" << best_wcss << "\t" << aic << "\t" << bic << "\n";

        vector<vector<pair<double, string>>> cluster_contents(k);
        for (size_t i = 0; i < data.size(); ++i) {
            int cluster_id = best_assignments[i];
            double dist = euclidean(data[i], best_centroids[cluster_id]);
            cluster_contents[cluster_id].emplace_back(dist, genome_labels[i]);
        }

        string base_name = filename.substr(0, filename.find_last_of('.'));
        ofstream out(base_name + "_k" + to_string(k) + ".txt");
        if (out.is_open()) {
            out << "Detailed output for k=" << k << "\n\n";
            out << "List of data points in each cluster preceded by distance to centroid:\n\n";

            for (int cluster_id = 0; cluster_id < k; ++cluster_id) {
                out << "Cluster " << cluster_id + 1 << ":\n\n";
                for (const auto& entry : cluster_contents[cluster_id]) {
                    out << "    " << fixed << setprecision(4) << entry.first << "\t" << entry.second << "\n";
                }
                out << "\n";
            }

            out << "Centroid coordinates:\n\n";
            for (int cluster_id = 0; cluster_id < k; ++cluster_id) {
                out << "Cluster " << cluster_id + 1 << ":\t";
                for (double val : best_centroids[cluster_id])
                    out << fixed << setprecision(5) << val << "\t";
                out << "\n";
            }

            out << "\nMutual pairwise distances among centroids:\n\n";
            for (int i = 0; i < k; ++i) {
                for (int j = 0; j < k; ++j) {
                    double dist = euclidean(best_centroids[i], best_centroids[j]);
                    out << fixed << setprecision(4) << dist << "\t";
                }
                out << "Cluster " << i + 1 << "\n";
            }

            out.close();
        }
    }

    return 0;
}