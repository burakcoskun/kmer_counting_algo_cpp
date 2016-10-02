#include <cstdio>
#include <cstdlib>
#include <vector>
#include <string>
#include <fstream>
#include <iostream>
#include <unordered_map>
#include "kmer_algo.h"

#define mod_prime_num 	10000019

using namespace std;

char conv[128];


int k_mer_length;
int topCount;
long long tot_length;
string all_sequences;

vector<vector<pair<long long ,long long> > > m_hash_map;
vector<pair<long long, long long> > res;

void set_conv(){
	conv['A'] = 0;
	conv['C'] = 1;
	conv['T'] = 2;
	conv['G'] = 3;
}

void insert_hash(long long int pos,long long int hash_val,string &sequence){
	// printf("%lld\n", hash_val);
	for(int i = 0 ; i < m_hash_map[hash_val].size(); ++i){
		int j , k = m_hash_map[hash_val][i].second;
		for(j = 0 ; j < k_mer_length ; ++j)
			if(sequence[pos+j] != all_sequences[k + j])
				break;

		if(j == k_mer_length){
			m_hash_map[hash_val][i].first ++; 
			return ;
		}
	}

	m_hash_map[hash_val].push_back(make_pair(1,tot_length + pos));


}

void calc_hash(string &sequence ){
	if(sequence.length() < k_mer_length)
		return ; 
	long long int pow4 = 1;
	long long int hash_val = 0;
	for(int i = k_mer_length -1 ;  i >= 0 ; --i){
		hash_val += ((pow4 * conv[sequence[i]]) % (mod_prime_num));
		hash_val %= mod_prime_num;

		if( i > 0){
			pow4 *= 4;
			pow4 %= mod_prime_num;
		}
	}

	insert_hash((long long)0,hash_val,sequence);

	for(int i = k_mer_length ; i< sequence.length() ; ++i){
		hash_val -= ((conv[sequence[i-k_mer_length]] * pow4)%(mod_prime_num));
		hash_val += mod_prime_num;
		hash_val %= mod_prime_num; 
		hash_val *= 4;
		hash_val %= mod_prime_num;
		hash_val += conv[sequence[i]];
		hash_val %= (mod_prime_num);
		insert_hash(i-k_mer_length+1,hash_val,sequence);
	}

}

void read_input(char *file_name){
	string sequence;
	ifstream input(file_name);

	while( !input.eof() ){

		getline(input,sequence,'\n');
		getline(input,sequence,'\n');
		all_sequences += sequence;
		calc_hash(sequence);
		// cout << all_sequences[all_sequences.size()-1] << '\n';

		// cout << sequence + '\n';
		tot_length += sequence.size();

		getline(input,sequence,'\n');
		getline(input,sequence,'\n');

		// if(c % 10000 == 0)
			// printf("%d\n", c);

	}
	// all_sequences.pop_back();
}

void get_results(){
	for(int i = 0 ; i < mod_prime_num ; ++i)
		for(int j = 0 ; j < m_hash_map[i].size(); ++j)
			res.push_back(m_hash_map[i][j]);
	sort(res.begin(),res.end(),greater<pair<int,int> >());
	for(int i = 0 ; i < topCount ; ++i){
		for(int j = 0 ; j < k_mer_length; ++j)
			printf("%c",all_sequences[res[i].second+j]);
		printf("\n");
	}
}

int main(int argc, char **argv){
	m_hash_map.resize(mod_prime_num);
	set_conv();
	k_mer_length = atoi(argv[2]);
	topCount = atoi(argv[3]);
	read_input(argv[1]);
	get_results();
	return 0;
}