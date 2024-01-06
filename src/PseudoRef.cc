#include <iostream>
#include <fstream>
#include <sstream>
#include <stdlib.h>
#include <algorithm>
#include <vector>
#include <unordered_map>
#include "utility.h"
#include "p_utility.h"

using namespace std;
map<string,bool> reads_in_ref;
unordered_map<kmer_int_type_t, int> minimizers_in_ref;
map<kmer_int_type_t, vector<int> > minimizers_in_ref2;
string output_dir="";
int current_ref_idx = 0;

bool check_sequence_drna(int idx, string& sequence)
{
    //return true;
    //for simulation perfect reference
    /*
    string id_ = dataid[idx];
    string::size_type pos = id_.find("_");
    if(pos == id_.npos) return false;
    string tid = id_.substr(0,pos);
    if(reads_in_ref.find(tid) != reads_in_ref.end()) return false;
    reads_in_ref[tid] = true;
    
    return true;
    */
    //if(sequence.length() < 2500) return false;
    if(sequence.length() < 1500) return false;
    //if(sequence.length() < 3000) return false;
    //cout<<"Check sequence..."<<'\n';
    vector<kmer_int_type_t> all_kmers;
    vector<kmer_int_type_t> all_minimizers;
    for (size_t j = 0; j <= sequence.length()-kmer_length; j++) 
    {
        const std::string& kmer = sequence.substr(j, kmer_length);

        kmer_int_type_t kmer_val = kmer_to_intval(kmer, kmer_length);
    /*
    	if (DS_MODE)
          kmer_val = get_DS_kmer_val(kmer_val, kmer_length);
*/
        all_kmers.push_back(kmer_val);

        if(j < g_window_length - 1) continue;
        vector<kmer_int_type_t>::iterator min_it
                = min_element(all_kmers.end()-g_window_length, all_kmers.end());

        kmer_val = *min_it;

        //kmer_to_reads[kmer_val].push_back(i);
        //reads_to_kmer[i].push_back(kmer_val);
	all_minimizers.push_back(kmer_val);
     }
    
    sort(all_minimizers.begin(),all_minimizers.end());
    all_minimizers.erase(unique(all_minimizers.begin(),all_minimizers.end()),all_minimizers.end());
    int in = 0, out = 0;
    for(size_t i=0;i<all_minimizers.size();i++)
    {
        if(minimizers_in_ref.find(all_minimizers[i]) == minimizers_in_ref.end()) 
		out++;
	else 
		in++;
    }
    double r = (1.0*all_minimizers.size())/(1.0*sequence.length());
    if(in<40 && in<0.3*out && r>0.15)
    {

     for(size_t i=0;i<all_minimizers.size();i++)
     {
        //minimizers_in_ref[all_minimizers[i]] = 1;
	if(minimizers_in_ref.find(all_minimizers[i]) == minimizers_in_ref.end())
		minimizers_in_ref[all_minimizers[i]] = 1;
     }
     /*
     cout<<dataid[idx]<<" minimizers-vs-length: "
	<<all_minimizers.size()<<" / "<<sequence.length()<<" "
	<<(1.0*all_minimizers.size())/(1.0*sequence.length())
	<<" in-vs-out: "<<in<<" "<<out
	<<'\n';
	*/
     return true;
    }
    return false;
}
bool check_sequence3(int idx, string& sequence)
{
    //return true;
    //for simulation perfect reference
    /*
    string id_ = dataid[idx];
    string::size_type pos = id_.find("_");
    if(pos == id_.npos) return false;
    string tid = id_.substr(0,pos);
    if(reads_in_ref.find(tid) != reads_in_ref.end()) return false;
    reads_in_ref[tid] = true;
    
    return true;
    */
    if(sequence.length() > 5000) return false;
    if(sequence.length() < 500) return false;
    //cout<<"Check sequence..."<<'\n';
    vector<kmer_int_type_t> all_kmers;
    vector<kmer_int_type_t> all_minimizers;
    for (size_t j = 0; j <= sequence.length()-kmer_length; ++j) 
    {

        const std::string& kmer = sequence.substr(j, kmer_length);

        kmer_int_type_t kmer_val = kmer_to_intval(kmer, kmer_length);
    
    	if (DS_MODE)
          kmer_val = get_DS_kmer_val(kmer_val, kmer_length);

        all_kmers.push_back(kmer_val);

        if(j < g_window_length - 1) continue;
        vector<kmer_int_type_t>::iterator min_it
                = min_element(all_kmers.end()-g_window_length, all_kmers.end());

        kmer_val = *min_it;

        //kmer_to_reads[kmer_val].push_back(i);
        //reads_to_kmer[i].push_back(kmer_val);
	all_minimizers.push_back(kmer_val);
     }
    sort(all_minimizers.begin(),all_minimizers.end());
    all_minimizers.erase(unique(all_minimizers.begin(),all_minimizers.end()),all_minimizers.end());

    int in = 0, out = 0;
    for(size_t i=0;i<all_minimizers.size();i++)
    {
        if(minimizers_in_ref.find(all_minimizers[i]) == minimizers_in_ref.end()) 
		out++;
	else 
		in++;
    }
    double r = (1.0*all_minimizers.size())/(1.0*sequence.length());
    bool c0 = in<40 && r>0.15;
    bool c1 = in>=40 && in<100 && out>3*in && r>0.15;
    //if(in<100 && in<0.3*out && r>0.15)
    if(c0 || c1)
    {

     for(size_t i=0;i<all_minimizers.size();i++)
     {
        minimizers_in_ref[all_minimizers[i]] = 1;
     }
     /*
     cout<<dataid[idx]<<" minimizers-vs-length: "
	<<all_minimizers.size()<<" / "<<sequence.length()<<" "
	<<(1.0*all_minimizers.size())/(1.0*sequence.length())
	<<" in-vs-out: "<<in<<" "<<out
	<<'\n';
     */
     return true;
    }
    return false;
}
bool check_sequence2(int idx, string& sequence)
{
    //return true;
    //for simulation perfect reference
    /*
    string id_ = dataid[idx];
    string::size_type pos = id_.find("_");
    if(pos == id_.npos) return false;
    string tid = id_.substr(0,pos);
    if(reads_in_ref.find(tid) != reads_in_ref.end()) return false;
    reads_in_ref[tid] = true;
    
    return true;
    */
    if(sequence.length() > 5000) return false;
    if(sequence.length() < 800) return false;
    //cout<<"Check sequence..."<<'\n';
    vector<kmer_int_type_t> all_kmers;
    vector<kmer_int_type_t> all_minimizers;
    for (size_t j = 0; j <= sequence.length()-kmer_length; ++j) 
    {

        const std::string& kmer = sequence.substr(j, kmer_length);

        kmer_int_type_t kmer_val = kmer_to_intval(kmer, kmer_length);
    
    	if (DS_MODE)
          kmer_val = get_DS_kmer_val(kmer_val, kmer_length);

        all_kmers.push_back(kmer_val);

        if(j < g_window_length - 1) continue;
        vector<kmer_int_type_t>::iterator min_it
                = min_element(all_kmers.end()-g_window_length, all_kmers.end());

        kmer_val = *min_it;

        //kmer_to_reads[kmer_val].push_back(i);
        //reads_to_kmer[i].push_back(kmer_val);
	all_minimizers.push_back(kmer_val);
     }
    sort(all_minimizers.begin(),all_minimizers.end());
    all_minimizers.erase(unique(all_minimizers.begin(),all_minimizers.end()),all_minimizers.end());

    int in = 0, out = 0;
    for(size_t i=0;i<all_minimizers.size();i++)
    {
        if(minimizers_in_ref.find(all_minimizers[i]) == minimizers_in_ref.end()) 
		out++;
	else 
		in++;
    }
    double r = (1.0*all_minimizers.size())/(1.0*sequence.length());
    bool c0 = in<30 && r>0.15;
    bool c1 = in>=30 && in<100 && out>3*in && r>0.15;
    //if(in<100 && in<0.3*out && r>0.15)
    if(c0 || c1)
    {

     for(size_t i=0;i<all_minimizers.size();i++)
     {
        minimizers_in_ref[all_minimizers[i]] = 1;
     }
     /*
     cout<<dataid[idx]<<" minimizers-vs-length: "
	<<all_minimizers.size()<<" / "<<sequence.length()<<" "
	<<(1.0*all_minimizers.size())/(1.0*sequence.length())
	<<" in-vs-out: "<<in<<" "<<out
	<<'\n';
	*/
     return true;
    }
    return false;
}


int num_of_most_element(vector<int>& ref_shared_m)
{
    //cout<<"---"<<endl;
    if(ref_shared_m.empty()) return 0;
    vector<int> temp = ref_shared_m;
    sort(ref_shared_m.begin(),ref_shared_m.end());

    int numb = 1;
    int max = 0;
    //cout<<ref_shared_m[0]<<" ";
    for(size_t i=1;i<ref_shared_m.size();i++)
    {
	//cout<<ref_shared_m[i]<<" ";
        if(ref_shared_m[i] == ref_shared_m[i-1])
	{
	    numb++;
	}
	else{

	    if(numb > max ) max = numb;
	    numb = 1;
	}
    }
    if(numb > max ) max = numb;
    return max;
    //cout<<endl;
    ref_shared_m.erase(unique(ref_shared_m.begin(),ref_shared_m.end()),ref_shared_m.end());
    vector<int> num(ref_shared_m.size(),0);
    for(size_t i=0;i<ref_shared_m.size();i++)
    {
	//cout<<ref_shared_m[i]<<" ";
        for(size_t j=0;j<temp.size();j++)
	{
	    if(temp[j] == ref_shared_m[i]) num[i]++;
	}
    }
    //cout<<" &"<<endl;
    int a = *max_element(num.begin(),num.end());
    //cout<<max<<" "<<a<<endl;
    return a;

}
bool check_sequence(int idx, string& sequence)
{
    //return true;
    //for simulation perfect reference
    /*
    string id_ = dataid[idx];
    string::size_type pos = id_.find("_");
    if(pos == id_.npos) return false;
    string tid = id_.substr(0,pos);
    if(reads_in_ref.find(tid) != reads_in_ref.end()) return false;
    reads_in_ref[tid] = true;
    
    return true;
    */
    //if(sequence.length() < 2500) return false;
    if(sequence.length() < 1500) return false;
    //if(sequence.length() < 3000) return false;
    //cout<<"Check sequence..."<<'\n';
    vector<kmer_int_type_t> all_kmers;
    vector<kmer_int_type_t> all_minimizers;
    for (size_t j = 0; j <= sequence.length()-kmer_length; ++j) 
    {

        const std::string& kmer = sequence.substr(j, kmer_length);

        kmer_int_type_t kmer_val = kmer_to_intval(kmer, kmer_length);
    
    	if (DS_MODE)
          kmer_val = get_DS_kmer_val(kmer_val, kmer_length);

        all_kmers.push_back(kmer_val);

        if(j < g_window_length - 1) continue;
        vector<kmer_int_type_t>::iterator min_it
                = min_element(all_kmers.end()-g_window_length, all_kmers.end());

        kmer_val = *min_it;

        //kmer_to_reads[kmer_val].push_back(i);
        //reads_to_kmer[i].push_back(kmer_val);
	all_minimizers.push_back(kmer_val);
     }
    sort(all_minimizers.begin(),all_minimizers.end());
    all_minimizers.erase(unique(all_minimizers.begin(),all_minimizers.end()),all_minimizers.end());

    int in = 0, out = 0;
    vector<int> ref_shared_m;
    for(size_t i=0;i<all_minimizers.size();i++)
    {

	map<kmer_int_type_t, vector<int> >::iterator it = minimizers_in_ref2.find(all_minimizers[i]);
        if( it == minimizers_in_ref2.end()) 
		out++;
	else{
       		for(size_t j=0;j<it->second.size();j++)
		{
		    ref_shared_m.push_back(it->second[j]);
		}
		in++;
		
	}
    }
    double r = (1.0*all_minimizers.size())/(1.0*sequence.length());
    int max = num_of_most_element(ref_shared_m);//最多能与某一个ref share 多少个minimizer
    if(in<40 && in<0.3*out && r>0.15 && max <= 10)
    {

     for(size_t i=0;i<all_minimizers.size();i++)
     {
        if(minimizers_in_ref2.find(all_minimizers[i])== minimizers_in_ref2.end()){
	    vector<int> v(1,current_ref_idx);
	     minimizers_in_ref2[all_minimizers[i]] = v;
	}
	else minimizers_in_ref2[all_minimizers[i]].push_back(current_ref_idx);
        //minimizers_in_ref[all_minimizers[i]] = 1;
     }
     //sort(ref_shared_m.begin(),ref_shared_m.end());
     //ref_shared_m.erase(unique(ref_shared_m.begin(),ref_shared_m.end()),ref_shared_m.end());
     //int max = num_of_most_element(ref_shared_m);//最多能与某一个ref share 多少个minimizer
     /*
     cout<<dataid[idx]<<" minimizers-vs-length: "
	<<all_minimizers.size()<<" / "<<sequence.length()<<" "
	<<(1.0*all_minimizers.size())/(1.0*sequence.length())
	<<" in-vs-out: "<<in<<" "<<out<<" "
	<<(1.0*in)/(1.0*out)<<" " //6.26
	<<ref_shared_m.size()<<" "
	<<max<<" "
	<<'\n';
	*/
     for(size_t i=0;i<ref_shared_m.size();i++)
     {
	 break;
         cout<<ref_shared_m[i]<<" ";
     }
     //cout<<'\n';
     return true;
    }
    return false;
}
void get_seudo_reference()
{
    time_t beg = time(NULL);
    //std::cerr << "Begin getting seudo-reference ..." << std::endl;
    sort(dataidx_length_vec.begin(),dataidx_length_vec.end(),sorter);
    unsigned long long Length_ = 0;
    int Idx = 0;
    for(size_t i=0;i<dataidx_length_vec.size();i++)
    {
        Length_+= dataidx_length_vec[i].second;
	if(double((1.0*Length_)/(1.0*Length))>N)
	{
	   Idx = i;
	   break;
	}
    }
    string info="";
    string seq="";
    int i1=0,i2=0;
    int seudo_idx=0;
    ofstream out((output_dir+"/"+iFlag+".pseudo.fasta").c_str());
    ofstream out2((output_dir+"/"+iFlag+".pseudoref-reads.info").c_str());
    for(size_t i=0;i<=Idx;i++){
	int idx = dataidx_length_vec[i].first; //index 
	//if(i%4 !=0 ) continue;//trick
   	//if(i>0 && i%10 == 0)
   	if(i>0 && i%1 == 0 && seq != "")
	{
	  out<<">pseudo_"<<seudo_idx<<" "<<info<<'\n';
	  out<<seq<<'\n';
	  info="";
	  seq="";
	  seudo_idx++;
	  i1=0;
	  i2=0;
	}
	bool sequence_ok = false;
	if(DS_MODE){
         if(iFlag == "1st")sequence_ok = check_sequence(idx,data[idx]);
	 else if(iFlag == "2nd") sequence_ok = check_sequence2(idx,data[idx]);
	 else if(iFlag == "3rd") sequence_ok = check_sequence3(idx,data[idx]);
	}
	else
	{
	  if(iFlag == "1st")sequence_ok = check_sequence_drna(idx,data[idx]);
	  else if(iFlag == "2nd") sequence_ok = check_sequence2(idx,data[idx]);
	  else if(iFlag == "3rd") sequence_ok = check_sequence3(idx,data[idx]);
	}
	if(sequence_ok)
	{
	        seq.append(data[idx]);
	        i2=seq.length();
	        info.append(dataid[idx]);
	    	info.append("_"+to_string(i1)+"_"+to_string(i2)+"-");
		out2<<"pseudo_"<<seudo_idx<<" "<<dataID_detailedID_map[ dataid[idx] ]<<" : "<<i1<<" "<<i2<<'\n';
		i1=seq.length();
		current_ref_idx++;
	}
	
    }

    out.close();
    out2.close();

}
int main(int argc,char*argv[])
{
    if(argc == 1) {
        cout<<"-"<<endl;
        cout<<"This is a simple program to get the pseudo reference from a fasta(or fastq) file."<<endl;
        cout<<"./exe reads.fasta fa(or fq for .fastq file) iFlag(e.g. 1st,2nd,3rd,4th..) flag(dran or cdna) outputdir"<<endl;
        cout<<"Two files will be generated: pseudo-iFlag.fasta and reads_in_pseudo.info"<<endl;
        cout<<"-"<<endl;
        return 0;
    }
    string filetype=argv[2];
    iFlag = argv[3];
    set();

    seqFlag = argv[4];
    output_dir=argv[5];

    if(seqFlag == "dRNA") DS_MODE = false;
    else if(seqFlag == "cDNA") DS_MODE = true;
    else if(seqFlag == "PacBio"){
	    PB_MODE = true;
	    DS_MODE = false;
    }

    if(filetype == "fq")
	load_data_fastq(argv[1]);
    else if(filetype == "fa")
	load_data_fasta(argv[1]);

    get_seudo_reference();
}
