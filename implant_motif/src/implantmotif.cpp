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

bool freePos(int& i, vector<int>& p) {
    for(int a : p)
        if(i == a) 
            return false;
    return true;
}

// from https://www.geeksforgeeks.org/hamming-distance-two-strings/
int hammingDist(DnaString str1, DnaString str2) {
    int i = 0, count = 0;
    for(int x = 0; x < length(str1); x++) {
        if (str1[i] != str2[i])
            count++;
        i++;
    }
    return count;
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

    vector<vector<double>> normalized_gcn4_jaspar_matrix = normalizeMatrix(gcn4_jaspar_matrix);

    printMatrix(normalized_gcn4_jaspar_matrix);

    // Generate 20 randomized DNA strings based on the input matrix
    for (int i = 0; i < t; i++) {
        std::string randomString = generateRandomMotifFromMatrix(normalized_gcn4_jaspar_matrix);
        std::cout << randomString << std::endl;
    }


    string filename = "implanted_gcn4.fasta"; // Output filename

    ofstream outfile(filename);
    if (!outfile) {
        cerr << "Error: could not open file \"" << filename << "\" for writing.\n";
        return 1;
    }

    for (int i = 0; i < t; i++) {
        string seq = generateRandomSequence(n);
        outfile << ">seq_" << i << "\n";
        outfile << seq << "\n";
    }

    outfile.close();
    cout << "Generated " << t << " random sequences of length " << n << " and saved them to file \"" << filename << "\".\n";



    return 0;
}
