#include <iostream>
#include <fstream>
#include <sstream>
#include <stdlib.h>
#include <algorithm>
#include <vector>

using namespace std;
//./exe a.fastq,b.fastq,c.fastq output_dir/
//To generate a merged file for multiple input files
int main(int argc,char*argv[])
{
    if(argc == 1) {
        cerr<<"Please provide input files!(e.g.,file1.fastq,file2.fastq)"<<endl;
	return 0;
    }
    string files="";	
    files = argv[1];
    for(size_t i=0;i<files.length();i++) {
       if(files[i] == ',')
	       files[i] = ' ';
    }
    ofstream out(argv[2]);
    istringstream istr;
    istr.str(files);
    string file;
    int Index = 0;
    while(istr>>file)
    {
	cout<<file<<endl;
        ifstream in(file.c_str());
	string s1,s2;
	while(getline(in,s1))
	{
	       getline(in,s2);
	 	out<<s1<<"_sample"<<Index<<'\n';
		out<<s2<<'\n';
	}
	Index++;
	in.close();
    }
    out.close();

    return 0;
}
