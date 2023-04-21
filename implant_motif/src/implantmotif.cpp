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
std::string generateRandomString(std::vector<std::vector<double>>& matrix) {
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


    vector<vector<double>> gcn4_jaspar_matrix = {
            {30, 45, 0, 1, 87, 0, 1, 7, 82, 6, 16},
            {15, 15, 0, 0, 0, 62, 0, 78, 3, 20, 31},
            {26, 19, 0, 84, 2, 28, 2, 5, 2, 13, 13},
            {19, 11, 90, 5, 1, 0, 87, 0, 3, 51, 30}
    };

    vector<vector<double>> normalized_gcn4_jaspar_matrix = normalizeMatrix(gcn4_jaspar_matrix);

    printMatrix(normalized_gcn4_jaspar_matrix);

    // Generate 20 randomized DNA strings based on the input matrix
    for (int i = 0; i < 20; i++) {
        std::string randomString = generateRandomString(normalized_gcn4_jaspar_matrix);
        std::cout << randomString << std::endl;
    }

    //t sequences each of length n
    int t = stoi(argv[1]);
    int n = stoi(argv[2]);
    //int t = 20;
    //int n = 600;

    //motif
    int motif_length = stoi(argv[3]);
    int d = stoi(argv[4]);
    //int motif_length = 11; //l
    //int d = 2;




    /*
    //for the sake of reproducibility
    std::mt19937 rng;

    //from https://stackoverflow.com/questions/7560114/random-number-c-in-some-range
    std::random_device rd; // obtain a random number from hardware
    std::mt19937 eng(rd()); // seed the generator

    //randomly chosing at how many positions the planted motif will have mutations
    std::uniform_int_distribution<> mut_num_distr(0, d); //how many mutations
    std::uniform_int_distribution<> mut_pos_distr(0, motif_length-1); //at which position
    std::uniform_int_distribution<> motif_pos_distr(0, n - motif_length); //at which position
    std::uniform_int_distribution<> dna(0, 3); //Dna nucleotide

    DnaString planted_motif;
    for (int i = 0; i < motif_length; ++i)
            appendValue(planted_motif, Dna(dna(eng)));
    //randomly chose a planted motif
    cout << "planted_motif " << planted_motif << endl;

    //save the planted motif in a file
    ofstream planted_motif_File;
    ostringstream ossP;
    ossP << "syn_planted_motif_" << motif_length << "_" << d << ".csv";
    string planted_motif_FileName = ossP.str();
    planted_motif_File.open (planted_motif_FileName);
    planted_motif_File << "seq\tpos\tmotif\tmutations\tmut_pos" << endl;
    planted_motif_File << "X\tX\t" << planted_motif << "\t0\t";


    //save the sequences in a file
    ofstream fastaFile;
    ostringstream oss;
    oss << "syn_planted_motif_" << motif_length << "_" << d << ".fasta";
    string fastaFileName = oss.str();
    fastaFile.open (fastaFileName);
    
    //mutate the planted motif for each of the t sequences
    for (int seqNum = 0; seqNum < t; seqNum++){
        vector<int> mutation_positions;
        DnaString planted_motif_in_seq;
        planted_motif_in_seq = planted_motif;
        //randomly chose the number of positions at which the planted motif will have mutations
        //int number_of_mut = mut_num_distr(eng);
        int number_of_mut = d;
        if(number_of_mut > 0){
            //for each mutation
            for(int i=0; i<number_of_mut; ++i){
                //randomly chose the positions at which the mutation will occur
                int mutation_pos = mut_pos_distr(eng);
                while(!freePos(mutation_pos, mutation_positions)) {
                    mutation_pos = mut_pos_distr(eng);
                };
                mutation_positions.push_back(mutation_pos);
                DnaString r = Dna(dna(eng));
                while(planted_motif_in_seq[mutation_pos] == r) 
                    r = Dna(dna(eng));
                planted_motif_in_seq[mutation_pos] = r[0]; //the mutation
            }
        }
        //randomly choose the positions at which the motif will be in the sequence
        int pos_of_motif = motif_pos_distr(eng);
        DnaString sequence;
        for (int j = 0; j <= n; j++){
            if(j == pos_of_motif){
                sequence += planted_motif_in_seq;
            } else if(j > pos_of_motif &&(j < (pos_of_motif+length(planted_motif)))){
                continue;
            } else if (j < pos_of_motif || (j > (pos_of_motif+length(planted_motif)))){
                appendValue(sequence, Dna(dna(eng)));
            }
        }
        
        
        
        //saving the sequences in a fasta file
        ostringstream oss2;
        oss2 << ">seq" << seqNum << "\n";
        string seq_identifier = oss2.str();
        fastaFile << seq_identifier;

        fastaFile << sequence;

        ostringstream oss3;
        oss3 << "\n";
        string new_line = oss3.str();
        fastaFile << new_line;

        planted_motif_File << new_line;
        planted_motif_File << seqNum << "\t" << pos_of_motif << "\t";
        planted_motif_File << planted_motif_in_seq << "\t" <<  hammingDist(planted_motif, planted_motif_in_seq) << "\t";
        for(int i : mutation_positions)
            planted_motif_File << "|" << i;

        
    }
    
    //fastaFile << "Writing this to a file.\n";
    fastaFile.close();
    planted_motif_File.close();*/
    return 0;
}
