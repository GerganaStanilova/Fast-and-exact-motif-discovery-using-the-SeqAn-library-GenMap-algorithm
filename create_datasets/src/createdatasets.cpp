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

using namespace std;
using namespace seqan;

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

int main()
{
    

    //t sequences each of length n
    int t = 20;
    int n = 600;
    
    //(15,4)
    //int motif_length = 15; //l
    //int d = 4;

    //(14,4)
    //int motif_length = 14; //l
    //int d = 4;
    
    //(16,5) 
    //int motif_length = 16; //l
    //int d = 5;

    //(18,6)
    //int motif_length = 18; //l
    //int d = 6;
    
    //smaller dataset
    int motif_length = 11; //l
    int d = 2;

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
    ossP << "n_planted_motif_" << motif_length << "_" << d << ".csv";
    string planted_motif_FileName = ossP.str();
    planted_motif_File.open (planted_motif_FileName);
    planted_motif_File << "seq\tpos\tmotif\tmutations\tmut_pos" << endl;
    planted_motif_File << "X\tX\t" << planted_motif << "\t0\t";


    //save the sequences in a file
    ofstream fastaFile;
    ostringstream oss;
    oss << "n_exampleDataset_" << motif_length << "_" << d << ".fasta";
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
    planted_motif_File.close();
    return 0;
}
