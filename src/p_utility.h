#include <iostream>
#include <fstream>
#include <sstream>
#include <stdlib.h>
#include <algorithm>
#include <vector>

const int MAX_STR = 1024;
using namespace std;
typedef unsigned long long kmer_int_type_t;
int g_window_length = 10;
int kmer_length = 15;

vector< pair<int,int> > dataidx_length_vec;
vector<string> data;
vector<string> dataid;
map<string,string> dataID_detailedID_map;
map<string,int> dataID_idx_map;
unsigned long long Length = 0;

string Nx = "";
double N = 0;
bool sorter(pair<int,int>p1, pair<int,int> p2)
{
    return (p1.second>p2.second || 
		    (p1.second == p2.second) && (p1.first<p2.second));
}
void load_data_fasta(char* file)
{

    ifstream in;
    in.open(file);
    time_t beg = time(NULL);
    std::cerr << "Begin loading reads ..." << std::endl;
    istringstream istr;
    string s,temp,id;
    int i=0;
    while(getline(in,id))
    {
	getline(in,s);
	//getline(in,temp);
	//getline(in,temp);
	istr.str(id);
	string id_;
	istr>>id_;
	id_=id_.substr(1);
	istr.clear();
	dataID_detailedID_map[id_] = id;
	data.push_back(s);
	dataid.push_back(id_);
	pair<int,int> p = make_pair(i,s.length());
	dataidx_length_vec.push_back(p);
	dataID_idx_map[id_] = i;
	Length += s.length();
	i++;

    }
    time_t end = time(NULL);
    cerr << data.size()<<"("<<dataidx_length_vec.size()<<" "<<Length<<")" 
	    <<" reads have been loaded, "
         << "elapsed time: " << (end-beg) << " s)" << std::endl;
    in.close();
}
void load_data_fastq(char* file) 
{
    ifstream in;
    in.open(file);
    time_t beg = time(NULL);
    std::cerr << "Begin loading reads ..." << std::endl;
    istringstream istr;
    string s,temp,id;
    int i=0;
    while(getline(in,id))
    {
	getline(in,s);
	getline(in,temp);
	getline(in,temp);
	istr.str(id);
	string id_;
	istr>>id_;
	id_=id_.substr(1);
	dataID_detailedID_map[id_] = id;
	istr.clear();
	data.push_back(s);
	dataid.push_back(id_);
	pair<int,int> p = make_pair(i,s.length());
	dataidx_length_vec.push_back(p);
	dataID_idx_map[id_] = i; //new
	Length += s.length();
	i++;

    }
    time_t end = time(NULL);
    cerr << data.size()<<"("<<dataidx_length_vec.size()<<" "<<Length<<")" 
	    <<" reads have been loaded, "
         << "elapsed time: " << (end-beg) << " s)" << std::endl;
    in.close();
    return;
}
