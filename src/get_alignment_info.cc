
#include "GArgs.h"
#include "GStr.h"
#include "gff.h"
#include "GBam.h"
#include "GBitVec.h"
#include "time.h"
//#include "tablemaker.h"
#include "GIntHash.hh"
#include "tmerge.h"

#include <iostream>
#include <string>
#include <fstream>
using namespace std;
//#define GMEMTRACE 1  //debugging mem allocation
ofstream outReadInfo;
TInputFiles bamreader;
GStr tmp_path;
bool debugMode=false;
bool verbose=false;
bool mergeMode = false; //--merge option
bool keepTempFiles = false; //--keeptmp
bool forceBAM = false; //useful for stdin
int main(int argc, char* argv[]) {
	//cout<<"Begin getting read-to-seudoRef alignment..."<<endl;
	bamreader.Add(argv[1]);
	outReadInfo.open(argv[2]);
	outReadInfo<<"read_id chromosome strand start_pos end_pos alignment_flag"<<endl;
/* yuting test
 *
 GVec<int> g_vec;
 g_vec.cPush(10);
 g_vec.cPush(15);
 for(int i=0;i<g_vec.Count();i++)
 {
    cout<<g_vec[i]<<endl;
 }
return 0;
*/
 // == Process arguments.
 /* yuting
*/

 int currentstart=0, currentend=0;
 int ng_start=0;
 int ng_end=-1;
 int ng=0;
 GStr lastref;
 bool no_ref_used=true;
 int lastref_id=-1; //last seen gseq_id

int bamcount=bamreader.start(); ////setup and open input files
 GBamRecord* brec=NULL;
 bool more_alns=true;
 int prev_pos=0;
 bool skipGseq=false;
 while (more_alns) {   //AAA
	 bool chr_changed=false;
	 int pos=0;
	 const char* refseqName=NULL;
	 char xstrand=0;
	 int nh=1;
	 int hi=0;
	 int gseq_id=lastref_id;  //current chr id
	 bool new_bundle=false;
	 //delete brec;
	 if ((brec=bamreader.next())!=NULL) {
		 if (brec->isUnmapped()) continue;

		 //if (!brec->isPrimary()) continue; //yuting
		 //if(brec->refId() != brec->mate_refId()) continue; //yuting
		 if(brec->isSupplementary()) continue;
		 string align_flag="primary";
		 if(!brec->isPrimary()) align_flag = "secondary";

		 if (brec->start<1 || brec->mapped_len<10) {
			 if (verbose) GMessage("Warning: invalid mapping found for read %s (position=%d, mapped length=%d)\n",
					 brec->name(), brec->start, brec->mapped_len);
			 continue;
		 }

		 refseqName=brec->refName();
		 int s1=brec->tag_int("s1");
		 int s2=brec->tag_int("s2");
		 int ec =0;
		 ec = brec->exons.Count();
		 int mc=brec->tag_int("cm");//minimizer count
        	 //outReadInfo<<brec->name()<<" "<<refseqName<<" "<<xstrand<<" "<<brec->start<<" "<<brec->end<<" "<<align_flag<<endl;
        	 outReadInfo<<brec->name()<<" "<<refseqName<<" "<<xstrand<<" "<<brec->start<<" "<<brec->end<<" "<<align_flag
			 <<" "<<s1<<" "<<s2<<" "<<ec<<" "<<mc<<" - "<<brec->getCIGAR()<<endl;
		 continue;
	 }
	 else { //no more alignments
                 more_alns=false;
	 
	 }

 }
}
