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
int gidx = 0,tidx=0;
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

	maxGid=gid; 
    }
    gidx += maxGid;
}
int main(int argc,char*argv[])
{
    out.open(argv[argc-1]);
    out<<"";
    for(int i=1;i<argc-1;i++)
	    load_clusterFile(argv[i]);
    return 0;
}
