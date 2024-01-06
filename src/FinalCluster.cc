#include <iostream>
#include <map>
#include <string>
#include <vector>
#include <vector>
#include <algorithm>
#include <fstream>
#include <sstream>


using namespace std;
ofstream out;
ofstream out2;//expression matrix
map<string, vector<int> > gene_cluster_read_map;
map<string, vector<string> > readid_readseq_map;
map<string, bool> NS_reads_map;
int gidx = 0,tidx=0;

char* reads;
bool Mflag = 0;
string outdir;
int sampleNum;
string readsFormat="fq";
void load_clusterFile(char*file)
{
    ifstream in(file);
    string s;
    int maxGid, maxTid;
    while(getline(in,s))
    {
	string::size_type pos=s.find("gene_cluster_"); //s: 06T19:20:09Z,gene_cluster_0

	string read = s.substr(0,pos);

        string s0=s.substr(pos),id;//s0:gene_cluster_0
       // string::size_type pos2=s0.find(",transcript_cluster_");

    	//string s1 = s0.replace(pos2,20," ");
        string s2 = s0.replace(0,13,"");
        istringstream istr_,istr;
        istr_.str(s2);
        int gid, tid;
        istr_>>gid;
	//cout<<read<<" "<<gid + gidx<<" "<<gid + gidx<<endl;
	out<<read<<"gene_cluster_"<<gid + gidx<<endl;//",transcript_cluster_"<<gid+gidx<<endl;
	//out<<gid + gidx<<"	"<<read.substr(1,read.length()-3)<<endl;//",transcript_cluster_"<<gid+gidx<<endl;
	string rid = read.substr(1),rid_;
	istringstream istr1;
	istr1.str(rid);
	istr1>>rid_;
	ofstream out;
	if(readsFormat == "fq")
	  	out.open((outdir+"/fastq_files/"+to_string(gid+gidx)+".fastq"),ios::app);
	else if (readsFormat == "fa")
		out.open((outdir+"/fastq_files/"+to_string(gid+gidx)+".fasta"),ios::app);

	NS_reads_map[rid_] = true;

	map<string, vector<string> >::iterator it=readid_readseq_map.find(rid_);
	if(it != readid_readseq_map.end()){
	    for(size_t i = 0;i<it->second.size();i++)
	    {
	        out<<it->second[i]<<'\n';
	    }
	}

	if(Mflag)
	{
	  string Gcluster = "gene_cluster_" + to_string(gid+gidx);
	  if(gene_cluster_read_map.find(Gcluster) == gene_cluster_read_map.end())
	  {
	    vector<int>v(sampleNum,0);
	    size_t lastPos = read.find_last_of('_');
	    int number = stoi(read.substr(lastPos+7));//***_sample0
	    v[number] = 1;
	    gene_cluster_read_map[Gcluster] = v;
	  }
	  else{
	    size_t lastPos = read.find_last_of('_');
	    int number = stoi(read.substr(lastPos+7));//***_sample0
	    gene_cluster_read_map[Gcluster][number] += 1;
	  }
	
	 }



	maxGid=gid; 
    }
    gidx += maxGid;
}
int get_sample_number(char*reads)
{
    ifstream in(reads);
    istringstream istr;
    string s1,s2,s3,s4;
    while(getline(in,s1))
    {
	getline(in,s2);
	getline(in,s3);
	getline(in,s4);
    }
    size_t lastPos = s3.find_last_of('_');
    int number = stoi(s3.substr(lastPos+7));//***_sample0
    return number;
}
void FQorFA(char*reads)
{
    bool flag = true;
    ifstream in(reads);
    string s;
    getline(in,s);
    if(s[0] == '>') readsFormat = "fa";
    else if(s[0] == '@') readsFormat = "fq";
    in.close();
    return;
}
void get_singleton_reads()
{
    ofstream out(outdir+"/GeLuster_singleton.tsv");
    map<string, vector<string> >::iterator it = readid_readseq_map.begin();
    int i = 0;
    for(;it != readid_readseq_map.end();it++)
    {
        if(NS_reads_map.find(it->first) == NS_reads_map.end())
	{
	    out<<it->second.front()<<" ,singleton_cluster"<<endl;
	    i++;
	}
    }
    return;
}
void load_read(char*reads)
{
    ifstream in(reads);
    istringstream istr;
    string s1,s2,s3,s4;
    if(readsFormat == "fq")
    {
      while(getline(in,s1))
      {
	getline(in,s2);
	getline(in,s3);
	getline(in,s4);
	string id = s1.substr(1),id_;
	istringstream istr;
	istr.str(id);
	istr>>id_;
	vector<string> v;
	v.push_back(s1);
	v.push_back(s2);
	v.push_back(s3);
	v.push_back(s4);
	readid_readseq_map[id_] = v;
      }
    }
    else if(readsFormat == "fa")
    {
      while(getline(in,s1))
      {
          getline(in,s2);
	  string id = s1.substr(1),id_;
	  istringstream istr;
	  istr.str(id);
	  istr>>id_;
	  vector<string> v;
	  v.push_back(s1);
	  v.push_back(s2);
	  readid_readseq_map[id_] = v;
      }
    }
    return;
}
int main(int argc,char*argv[])
{
    char* reads = argv[argc-1];//reads 
    string F = argv[argc-2];//one files or more
    if(F == "0") Mflag = false; 
    else Mflag = true;
    outdir = argv[argc-3];

    if(!Mflag)//one input reads
    {
      FQorFA(reads);//ensure the format of the reads;
      load_read(reads);
      out.open((outdir+"/GeLuster.tsv").c_str());//GeLuster.tsv
      out<<"";
      for(int i=1;i<argc-3;i++)
	    load_clusterFile(argv[i]);
      get_singleton_reads();
    }
    else//multiple input reads
    {
	sampleNum = get_sample_number(reads)+1;
	cerr<<sampleNum<<endl;
	out.open((outdir+"/GeLuster.tsv").c_str());//GeLuster.tsv
	out<<"";
	for(int i=1;i<argc-3;i++)
		load_clusterFile(argv[i]);
	string Matrix = argv[argc-2];
	Matrix = outdir + "/GeLuster.gene-expression-proxy.matrix";

	out2.open(Matrix.c_str());
	out2<<"clusters";
	for(int i=0;i<sampleNum;i++) out2<<" sample_"<<i;
	out2<<endl;
	map<string, vector<int> >::iterator it = gene_cluster_read_map.begin();
	for(;it != gene_cluster_read_map.end();it++)
	{
	    out2<<it->first;
	    for(int i=0;i<it->second.size();i++) 
		    out2<<" "<<it->second[i];
	    out2<<endl;
	}
    }
    return 0;
}
