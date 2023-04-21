#include <iostream>
#include <algorithm> 
#include <cstdlib>
#include <vector> 
#include <cmath>
#include <random>
#include <fstream>
#include <string>
#include <sstream>
#include <seqan/arg_parse.h>
#include <seqan/index.h>
#include <seqan/seq_io.h>
#include <seqan/sequence.h>
#include <seqan/basic.h>
#include <seqan/stream.h>
#include <seqan/file.h>
#include <seqan/modifier.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <sstream>
#include <map>

using namespace std;
using namespace seqan;

vector<vector<double>> normalizeMatrix(const vector<vector<double>>& matrix) {
    // Find the maximum value in the matrix
    double maxVal = matrix[0][0];
    for (int i = 0; i < matrix.size(); i++) {
        for (int j = 0; j < matrix[i].size(); j++) {
            if (matrix[i][j] > maxVal) {
                maxVal = matrix[i][j];
            }
        }
    }

    // Normalize the matrix by dividing each element by the maximum value
    vector<vector<double>> normalizedMatrix;
    for (int i = 0; i < matrix.size(); i++) {
        vector<double> row;
        for (int j = 0; j < matrix[i].size(); j++) {
            row.push_back(matrix[i][j] / maxVal);
        }
        normalizedMatrix.push_back(row);
    }

    return normalizedMatrix;
}

void printMatrix(const vector<vector<double>>& matrix) {
    for (int i = 0; i < matrix.size(); i++) {
        for (int j = 0; j < matrix[i].size(); j++) {
            cout << matrix[i][j] << " ";
        }
        cout << endl;
    }
}

//Function to generate a randomized string based on the input matrix
std::string generateRandomMotifFromMatrix(std::vector<std::vector<double>>& matrix) {
    std::string dnaString;
    std::random_device rd;
    std::mt19937 gen(rd());

    for (int i = 0; i < matrix[0].size(); i++) {
        std::vector<double> probabilities = {matrix[0][i], matrix[1][i], matrix[2][i], matrix[3][i]};
        std::discrete_distribution<int> distribution(probabilities.begin(), probabilities.end());
        char nucleotide;
        switch (distribution(gen)) {
            case 0:
                nucleotide = 'A';
                break;
            case 1:
                nucleotide = 'C';
                break;
            case 2:
                nucleotide = 'G';
                break;
            case 3:
                nucleotide = 'T';
                break;
        }
        dnaString += nucleotide;
    }
    return dnaString;
}

std::string generateRandomSequence(int n) {
    std::string dnaString;
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_int_distribution<int> distribution(0, 3); // A, C, G, T

    for (int i = 0; i < n; i++) {
        int nucleotideIndex = distribution(gen);
        char nucleotide;
        switch (nucleotideIndex) {
            case 0:
                nucleotide = 'A';
                break;
            case 1:
                nucleotide = 'C';
                break;
            case 2:
                nucleotide = 'G';
                break;
            case 3:
                nucleotide = 'T';
                break;
        }
        dnaString += nucleotide;
    }
    return dnaString;
}
std::string getConsensusSequence(const std::vector<std::string>& motifs) {
    // Determine length of motifs
    const int motifLength = motifs[0].length();

    // Create a matrix to hold the counts of each nucleotide at each position
    std::vector<std::vector<int>> counts(motifLength, std::vector<int>(4, 0)); // 4 = number of nucleotides

    // Fill the counts matrix
    for (const auto& motif : motifs) {
        for (int i = 0; i < motifLength; i++) {
            switch (motif[i]) {
                case 'A': counts[i][0]++; break;
                case 'C': counts[i][1]++; break;
                case 'G': counts[i][2]++; break;
                case 'T': counts[i][3]++; break;
                default: break;
            }
        }
    }

    // Determine the most frequent nucleotide at each position
    std::string consensus(motifLength, 'N'); // 'N' = undefined nucleotide
    for (int i = 0; i < motifLength; i++) {
        const int maxCount = *std::max_element(counts[i].begin(), counts[i].end());
        if (maxCount > 0) {
            const int index = std::distance(counts[i].begin(), std::find(counts[i].begin(), counts[i].end(), maxCount));
            switch (index) {
                case 0: consensus[i] = 'A'; break;
                case 1: consensus[i] = 'C'; break;
                case 2: consensus[i] = 'G'; break;
                case 3: consensus[i] = 'T'; break;
                default: break;
            }
        }
    }

    return consensus;
}

int main(int argc, char const ** argv)
{

    //t sequences each of length n
    int t = stoi(argv[1]);
    int n = stoi(argv[2]);
    //int t = 20;
    //int n = 600;

    vector<vector<double>> gcn4_jaspar_matrix = {
            {30, 45, 0, 1, 87, 0, 1, 7, 82, 6, 16},
            {15, 15, 0, 0, 0, 62, 0, 78, 3, 20, 31},
            {26, 19, 0, 84, 2, 28, 2, 5, 2, 13, 13},
            {19, 11, 90, 5, 1, 0, 87, 0, 3, 51, 30}
    };
    int motif_length = gcn4_jaspar_matrix[0].size();

    vector<vector<double>> normalized_gcn4_jaspar_matrix = normalizeMatrix(gcn4_jaspar_matrix);

    printMatrix(normalized_gcn4_jaspar_matrix);


    // Generate motifs
    vector<string> motifs;
    for (int i = 0; i < t; i++) {
        motifs.push_back(generateRandomMotifFromMatrix(normalized_gcn4_jaspar_matrix));
    }

    // Generate random sequences
    vector<string> sequences;
    for (int i = 0; i < t; i++) {
        sequences.push_back(generateRandomSequence(n));
    }

    // Implant motifs in each sequence and record positions and motifs
    vector<tuple<int, int, string>> motif_pos;
    for (int i = 0; i < t; i++) {
        int pos = rand() % (n - 20); // choose a random position to implant the motif

        for (int i = 0; i < motifs.size(); i++) {
            string motif = motifs[i];
            sequences[i].replace(pos, 20, motif); // implant motif
            motif_pos.push_back(make_tuple(i, pos, motif)); // record position and motif
        }

    }

    // Write sequences to FASTA file
    ofstream fastaFile("gcn4_sequences.fasta");
    for (int i = 0; i < t; i++) {
        fastaFile << ">seq_" << i << endl;
        fastaFile << sequences[i] << endl;
    }
    fastaFile.close();

    string csvFilePath = "motifs.csv";
    string consensus = getConsensusSequence(motifs);

    // Write motifs to CSV file
    ofstream csvFile("implanted_motifs.csv");
    csvFile << "seq\tpos\tmotif" << endl;
    csvFile << "X\tX\t" << consensus << endl; // write consensus sequence
    for (auto m : motif_pos) {
        csvFile << get<0>(m) << "\t" << get<1>(m) << "\t" << get<2>(m) << endl;
    }
    csvFile.close();





    return 0;
}
