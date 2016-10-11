#include <vector>
#include <fstream>
#include <string>
#include <queue>

using namespace std;

#include "kmer_algo.h"

#define mod_prime_num 	10000019


char conv[128];


int k_mer_length;
int topCount;
int tot_length;
string all_sequences;

typedef pair<int,int>  PAIR;

vector<vector<PAIR > > m_hash_map;
priority_queue<PAIR , vector<PAIR >, greater<PAIR > > res;

void set_conv(){
	conv['A'] = 0;
	conv['C'] = 1;
	conv['T'] = 2;
	conv['G'] = 3;
}

void insert_hash(int  pos,int  hash_val,string &sequence){
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
	int  pow4 = 1;
	int  hash_val = 0;
	for(int i = k_mer_length -1 ;  i >= 0 ; --i){
		hash_val += ((pow4 * conv[sequence[i]]) % (mod_prime_num));
		hash_val %= mod_prime_num;

		if( i > 0){
			pow4 *= 4;
			pow4 %= mod_prime_num;
		}
	}

	insert_hash((int)0,hash_val,sequence);

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
		tot_length += sequence.size();

		getline(input,sequence,'\n');
		getline(input,sequence,'\n');

	}
}

void get_results(){
	for(int i = 0 ; i < mod_prime_num ; ++i)
		for(int j = 0 ; j < m_hash_map[i].size(); ++j){
			if(res.size() < topCount) 
				res.push(m_hash_map[i][j]);
			else if(res.top().first < m_hash_map[i][j].first){
				res.pop();
				res.push(m_hash_map[i][j]);
			}
		}
	while( !res.empty() ){
		for(int j = 0 ; j < k_mer_length; ++j)
			printf("%c",all_sequences[res.top().second+j]);
		printf("\n");
		res.pop();
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