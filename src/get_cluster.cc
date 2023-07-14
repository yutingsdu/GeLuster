#include <iostream>
#include <vector>
#include <stdlib.h>
#include <algorithm>
#include <math.h>
#include <fstream>
#include <sstream>
#include <map>
#include "p_utility.h"
using namespace std;
string output_dir="";
ofstream out_my_cluster;//single
ofstream out_my_cluster2;//multi
ofstream out_singleton;
ofstream out_post;

int ClusterIdx = 0;
map<string,vector<string> > seudoref_alignedreads_map;//N50 //reads_in_seudoRef(id) -> alined_reads(id)
map<string,vector<string> > seudoref_alignedreads_map2; //N30 //reads_in_seudoRef(id) -> alined_reads(id)
map<string,vector<string> > reads_seudoref_map;
map<string, bool> reads_in_seudoref;
map<string,vector<string> > seudoref_seudoreads_map;//reads_in_seudoRef(id1) -> vec( reads_in_seudoRef)(id2) id2 mapped to id1
map<string,vector<string> > seudoreads_seudoref_map;//reads_in_seudoRef(id2) -> vec( reads_in_seudoRef)(id1) id2 mapped to id1
map<string,pair<string,double> > seudoreads_seudoref_map2;//reads_in_seudoRef(id2) -> id2 mapped to id1(where id1 has max mapping score)
struct seudoref_info
{
    string rid;
    int spos;
    int epos;
};
map<string, vector<seudoref_info> > cluster_seudorefinfo_map; //seudo_x -> readid
//map<string, vector<seudoref_info> > cluster_seudorefinfo_map;
void add_reads(map<string,vector<string> >& seudoref_alignedreads_map,string seudorefid,string readsid)
{
    if(seudoref_alignedreads_map.find(seudorefid) == seudoref_alignedreads_map.end())
    {
        vector<string> v(1,readsid);
	seudoref_alignedreads_map[seudorefid] = v;
    }
    else seudoref_alignedreads_map[seudorefid].push_back(readsid);
/*
    if(reads_seudoref_map.find(readsid) == reads_seudoref_map.end())
    {
        vector<string> v(1,seudorefid);
	reads_seudoref_map[readsid] = v;
    }
    else reads_seudoref_map[readsid].push_back(seudorefid);
*/
    return;
}

void get_cluster_pro(map<string,vector<string> >& seudoref_alignedreads_map,//N50
		     map<string,vector<string> >& seudoref_alignedreads_map2) //N30
{
    cerr<<"Begin getting cluster..."<<endl;
    map<string,bool> used_reads;
    map<string,bool> current_reads;
    vector<string> reads_in_current_cluster;

    int idx = 0;
    size_t i=0;
    while(1)
    {
	current_reads.clear();
	for(;i<dataidx_length_vec.size();i++)
	{
	    int readidx = dataidx_length_vec[i].first;
	    string readid = dataid[readidx];

	    if(!used_reads.empty() && used_reads.find(readid) != used_reads.end()) continue; 
	    reads_in_current_cluster.push_back(readid);
	    used_reads[readid] = true;
	    break;
	}
	if(reads_in_current_cluster.empty()) break;

	/* //temporarily don't consider
	for(size_t i=0;i<reads_in_current_cluster.size();i++)
	{
	    string r_ = reads_in_current_cluster[i];
	    //cout<<"a: "<<r_<<endl;
	    //if(used_reads.find(r_) != used_reads.end()) continue;
	    used_reads[r_] = true;
	    map<string,vector<string> >::iterator it_ = seudoref_alignedreads_map.find(r_);
	    if(it_ != seudoref_alignedreads_map.end())
	    {
	        for(size_t j=0;j<it_->second.size();j++)
		{
		    if(seudoref_alignedreads_map.find(it_->second[j]) != seudoref_alignedreads_map.end())
		    {
			//cout<<"B: "<<it_->second[j]<<endl;
			if(used_reads.find(it_->second[j]) == used_reads.end())
			{
			    reads_in_current_cluster.push_back(it_->second[j]);
			    used_reads[it_->second[j]] = true;
			}

		    }
		}
	    }

	}
	*/
	size_t size_ = reads_in_current_cluster.size();
	for(size_t i=0;i<size_;i++) //find in N30 Ref;remember: N50ref has more refs than that in N30Ref
	{
	    string r_ = reads_in_current_cluster[i];
	    map<string,vector<string> >::iterator it_ = seudoref_alignedreads_map2.find(r_);
	    if(it_ != seudoref_alignedreads_map2.end())
	    {
		for(size_t j=0;j<it_->second.size();j++)
		{
		    if(used_reads.find(it_->second[j]) == used_reads.end())
		    {
		        reads_in_current_cluster.push_back(it_->second[j]);
			used_reads[it_->second[j]] = true;
		    }
		}
	    }
	}
	for(size_t i=size_;i<reads_in_current_cluster.size();i++) //find in N5 Ref
	{
	    string r_ = reads_in_current_cluster[i];
	    map<string,vector<string> >::iterator it_ = seudoref_alignedreads_map.find(r_);
	    if(it_ != seudoref_alignedreads_map.end())
	    {
	        for(size_t j=0;j<it_->second.size();j++)
		{
		    if(used_reads.find(it_->second[j]) == used_reads.end())
		    {
		        reads_in_current_cluster.push_back(it_->second[j]);
			used_reads[it_->second[j]] = true;
		    }
		}
	    }
	}
	//sort(reads_in_current_cluster.begin(),reads_in_current_cluster.end());
	//reads_in_current_cluster.erase(unique(reads_in_current_cluster.begin(),reads_in_current_cluster.end()),reads_in_current_cluster.end());
	cout<<"Cluster_"<<idx<<" "<<reads_in_current_cluster.size()<<": ";
	for(size_t i=0;i<reads_in_current_cluster.size();i++) cout<<reads_in_current_cluster[i]<<" ";
      	cout<<endl;

	for(size_t i=0;i<reads_in_current_cluster.size();i++) {
	    out_my_cluster<<dataID_detailedID_map[reads_in_current_cluster[i]]<<",gene_cluster_"<<idx<<",transcript_cluster_"<<idx<<endl;

	}
	reads_in_current_cluster.clear();
	idx++;
    }
}
void get_cluster3(map<string,vector<string> >& seudoref_alignedreads_map)
{
     map<string,bool> used_ref;
     map<string,bool> used_reads;
     int idx = 0;
     if(Nx == "N80") idx = 100000;
     else if(Nx == "N90") idx = 500000;
     /*seudoref_seudoreads_map
      * key	value
      *seudo_0: seudo_1, seudo_2;
      *seudo_1: seudo_3, seudo_4;
      *I will get seudo_0 seudo_1 seudo_2 seudo_3 seudo_4 in one cluster
      *
      */
    map<string,vector<string> >::iterator it_ = seudoref_seudoreads_map.begin();
    for(;it_ != seudoref_seudoreads_map.end();it_++)
    {
	//break;
	string r1,r2;//r2 mapped to r1 (r1(ref)->r2(reads) r1<-r2)
	r1 = it_->first;
	//cout<<endl<<"Here:"<<r1<<endl;
	if(seudoreads_seudoref_map.find(r1) != seudoreads_seudoref_map.end()) continue;
	vector<string> seudo_ref;
	seudo_ref.push_back(r1);
	for(size_t k=0;k<seudo_ref.size();k++)
	{
	  string r1_a = seudo_ref[k];
	  map<string,vector<string> >::iterator it_a = seudoref_seudoreads_map.find(r1_a);
	  if(it_a == seudoref_seudoreads_map.end()) continue;
	  for(size_t i=0;i<it_a->second.size();i++)
	  {
 	    string x = it_a->second[i];
	    map<string,pair<string,double> >::iterator it__ = seudoreads_seudoref_map2.find(x);
	    if(it__->second.first != r1_a) continue;
	    seudo_ref.push_back(it_a->second[i]);
	  }
	}
	/*
	if(it_->second.size() != 1) continue; //new2
	r2 = it_->second.front();
	
	map<string,vector<string> >::iterator it__=seudoreads_seudoref_map.find(r2);
	if(it__->second.size() != 1) continue;
	string r_temp = it__->second.front();
	if(r_temp != r1) continue;
	if(r1 == r2) continue;
	
	if(seudoreads_seudoref_map.find(r1) != seudoreads_seudoref_map.end()
		|| seudoref_seudoreads_map.find(r2) != seudoref_seudoreads_map.end()
	  )
		continue;

	//cout<<"good: "<<r1<<" "<<r2<<" "<<r_temp<<endl;
	if(used_ref.find(it_->first) != used_ref.end()) continue;
	//vector<string> seudo_ref(1,it_->first);
	vector<string> seudo_ref;
	seudo_ref.push_back(r1);
	seudo_ref.push_back(r2);
	*/
	/*
	for(size_t i=0;i<it_->second.size();i++)
		seudo_ref.push_back(it_->second[i]);
	*/
	/*
	for(size_t i=0;i<seudo_ref.size();i++)
	{
	    if(used_ref.find(seudo_ref[i]) != used_ref.end()) continue;
	    used_ref[seudo_ref[i]] = true;
	    map<string,vector<string> >::iterator it1 = seudoref_seudoreads_map.find(seudo_ref[i]);
	    if(it1 != seudoref_seudoreads_map.end())
	    {
		if(used_ref.find(it1->first) == used_ref.end())
		    seudo_ref.push_back(it1->first);
	        for(size_t j=0;j<it1->second.size();j++)
		{
		    if(used_ref.find(it1->second[j]) == used_ref.end())
			    seudo_ref.push_back(it1->second[j]);
		    used_ref[it1->second[j]] = true;
		}
	    }
	    break;
	}


	sort(seudo_ref.begin(),seudo_ref.end());
	seudo_ref.erase(unique(seudo_ref.begin(),seudo_ref.end()),seudo_ref.end());
 	*/
	/*
	 *for each seudo get the reads that mapped to it
	 */
	for(size_t i=0;i<seudo_ref.size();i++)
	{
	    //cout<<seudo_ref[i]<<" * ";
	    map<string,vector<string> >::iterator it1 = seudoref_alignedreads_map.find(seudo_ref[i]);
	    if(it1 != seudoref_alignedreads_map.end())
	    {
		if(it1->second.size()>1)
		{
	          for(size_t j=0;j<it1->second.size();j++)
		  {
			if(used_ref.find(it1->second[j]) != used_ref.end()){
				cout<<"WARNING:"<<it1->second[j]<<endl;
			       	continue;
			}
			out_my_cluster2<<dataID_detailedID_map[it1->second[j]]<<" ,gene_cluster_"<<idx<<",transcript_cluster_"<<idx<<'\n';
			used_reads[it1->second[j]] = true;
			used_ref[it1->second[j]] = true;
		  }
		}
	    }
	}

	used_ref[r1] = true;
	used_ref[r2] = true;
	//cout<<endl;
	idx++;
    }
    
   //return;

    map<string,vector<string> >::iterator it = seudoref_alignedreads_map.begin();
    for(;it != seudoref_alignedreads_map.end();it++)
    {
	if(used_ref.find(it->first) != used_ref.end()) continue;
	 used_ref[it->first] = true;
	vector<string> seudo_ref;
	if(it->second.size()>1)
	{
	  for(size_t i=0;i<it->second.size();i++){
	
	    out_my_cluster2<<dataID_detailedID_map[it->second[i]]<<" ,gene_cluster_"<<idx<<",transcript_cluster_"<<idx<<'\n';
	    used_reads[it->second[i]] = true;
	    //if(seudoref_seudoreads_map.find(it->second[i]) != seudoref_seudoreads_map.end())
	    //{
	        //seudo_ref.push_back(it->second[i]);
		//used_ref[it->second[i]] = true;
	    //}
	  }
	}
	/*
	cout<<it->first<<endl;
	map<string,vector<string> >::iterator it_ = seudoref_seudoreads_map.find(it->first);
	if(it_ != seudoref_seudoreads_map.end())
	{
	  for(size_t i=0;i<it_->second.size();i++) seudo_ref.push_back(it_->second[i]);
	  for(size_t i=0;i<seudo_ref.size();i++)
	  {
	     if(used_ref.find(seudo_ref[i]) != used_ref.end()) continue;
	     used_ref[seudo_ref[i]] = true;
	     map<string,vector<string> >::iterator it1 = seudoref_alignedreads_map.find(seudo_ref[i]);
	     if(it1 != seudoref_alignedreads_map.end())
	     {
		 cout<<"hhh"<<endl;
	         for(size_t i=0;i<it1->second.size();i++){
		     out_my_cluster2<<dataID_detailedID_map[it1->second[i]]<<" ,gene_cluster_"<<idx<<",transcript_cluster_"<<idx<<'\n';
		 }
	     }
	  }
	}
	*/
	idx++;
    }
    if(Nx == "N90")
    {
     for(size_t i=0;i<dataidx_length_vec.size();i++)
     {
         int readidx = dataidx_length_vec[i].first;
	 string readid = dataid[readidx];
	 if(used_reads.find(readid) != used_reads.end()) continue;
	 int idx_ = dataID_idx_map[ readid ];

	 int L = dataidx_length_vec[i].second;
	 if(L > 500)
	 {
	     out_my_cluster<<dataID_detailedID_map[readid]<<" ,gene_cluster_"<<idx<<",transcript_cluster_"<<idx<<'\n';
	     out_singleton<<">"<<dataID_detailedID_map[readid].substr(1)<<'\n';
	     out_singleton<<data[idx_]<<'\n';
	     idx++;
	 }
	 else {
	     //out_post<<">"<<dataID_detailedID_map[readid].substr(1)<<'\n';
	     //out_post<<data[idx_]<<'\n';
	     out_post<<dataID_detailedID_map[readid]<<" ,gene_cluster_"<<idx<<",transcript_cluster_"<<idx<<'\n';
	     idx++;
	 }


     }
    }
    else{
     for(size_t i=0;i<dataidx_length_vec.size();i++)
     {
	int readidx = dataidx_length_vec[i].first;
	string readid = dataid[readidx];
        if(used_reads.find(readid) != used_reads.end()) continue;
        out_my_cluster<<dataID_detailedID_map[readid]<<" ,gene_cluster_"<<idx<<",transcript_cluster_"<<idx<<'\n';

        int idx_ = dataID_idx_map[ readid ];
        out_singleton<<">"<<dataID_detailedID_map[readid].substr(1)<<'\n';
        out_singleton<<data[idx_]<<'\n';
        idx++;

     }
    }

}
void get_cluster2(map<string,vector<string> >& seudoref_alignedreads_map)
{
    cerr<<"Begin getting cluster..."<<endl;
    map<string,bool> used_reads;
    map<string,bool> current_reads;
    vector<string> reads_in_current_cluster;

    int idx = 0;
    size_t i=0;
    while(1)
    {
	current_reads.clear();
	/*
	 *
	 *BUG:BUG: long read may mapped to short ones, suppose a mapped to b( a is a common read, b is seudo ref), 
	 *here a will be encountered first anx will be as a used read first; 
	 *
	 */
	for(;i<dataidx_length_vec.size();i++)
	{
	    int readidx = dataidx_length_vec[i].first;
	    string readid = dataid[readidx];
	    cout<<"  "<<readid<<endl;
	    if(!used_reads.empty() && used_reads.find(readid) != used_reads.end()) continue; 
	    reads_in_current_cluster.push_back(readid);
	    used_reads[readid] = true;
	    break;
	}
	if(reads_in_current_cluster.empty()) break;

	/* //temporarily don't consider
	for(size_t i=0;i<reads_in_current_cluster.size();i++)
	{
	    string r_ = reads_in_current_cluster[i];
	    //cout<<"a: "<<r_<<endl;
	    //if(used_reads.find(r_) != used_reads.end()) continue;
	    used_reads[r_] = true;
	    map<string,vector<string> >::iterator it_ = seudoref_alignedreads_map.find(r_);
	    if(it_ != seudoref_alignedreads_map.end())
	    {
	        for(size_t j=0;j<it_->second.size();j++)
		{
		    if(seudoref_alignedreads_map.find(it_->second[j]) != seudoref_alignedreads_map.end())
		    {
			//cout<<"B: "<<it_->second[j]<<endl;
			if(used_reads.find(it_->second[j]) == used_reads.end())
			{
			    reads_in_current_cluster.push_back(it_->second[j]);
			    used_reads[it_->second[j]] = true;
			}

		    }
		}
	    }

	}
	*/
	size_t size_ = reads_in_current_cluster.size();
	for(size_t i=0;i<size_;i++)
	{
	    string r_ = reads_in_current_cluster[i];
	    cout<<"hh: "<<r_<<endl;
	    //cout<<r_<<" "<<size_<<endl;
	    map<string,vector<string> >::iterator it_ = seudoref_alignedreads_map.find(r_);//readid->vector<readid>
	    if(it_ != seudoref_alignedreads_map.end())
	    {
		for(size_t j=0;j<it_->second.size();j++)
		{
		    if(used_reads.find(it_->second[j]) == used_reads.end())// && reads_in_seudoref.find(it_->second[j]) == reads_in_seudoref.end()) //NEW 
		    {
		        reads_in_current_cluster.push_back(it_->second[j]);
			used_reads[it_->second[j]] = true;
		    }
		}
	    }
	    size_ = reads_in_current_cluster.size();//NEW
	}
	//return;
	//sort(reads_in_current_cluster.begin(),reads_in_current_cluster.end());
	//reads_in_current_cluster.erase(unique(reads_in_current_cluster.begin(),reads_in_current_cluster.end()),reads_in_current_cluster.end());
	/*
	cout<<"Cluster_"<<idx<<" "<<reads_in_current_cluster.size()<<": ";
	for(size_t i=0;i<reads_in_current_cluster.size();i++) cout<<reads_in_current_cluster[i]<<" ";
      	cout<<endl;
	*/
	for(size_t i=0;i<reads_in_current_cluster.size();i++) {
	    out_my_cluster<<dataID_detailedID_map[reads_in_current_cluster[i]]<<" ,gene_cluster_"<<idx<<",transcript_cluster_"<<idx<<'\n';
	}
	if(reads_in_current_cluster.size() > 1)
	{
	  for(size_t i=0;i<reads_in_current_cluster.size();i++) {
	    out_my_cluster2<<dataID_detailedID_map[reads_in_current_cluster[i]]<<" ,gene_cluster_"<<idx<<",transcript_cluster_"<<idx<<'\n';
	  }
	}
	else if(reads_in_current_cluster.size() == 1)
	{
	    int idx = dataID_idx_map[ reads_in_current_cluster[0] ];
	    out_singleton<<">"<<dataID_detailedID_map[reads_in_current_cluster[0]].substr(1)<<'\n';
	    out_singleton<<data[idx]<<'\n';
	}
	reads_in_current_cluster.clear();
	idx++;
        //break;
    }
}
void get_cluster()
{
    int idx = 0;
    map<string,bool> used_reads_as_ref;
    map<string,bool> used_reads_as_query;
    map<string,bool> current_reads;
    vector<string> reads_in_current_cluster;

    map<string,vector<string> >::iterator it = reads_seudoref_map.begin();
    while(1)
    {
      current_reads.clear();
      for(;it != reads_seudoref_map.end();it++)
      {
        if(used_reads_as_query.find(it->first) != used_reads_as_query.end()) continue;
	reads_in_current_cluster.push_back(it->first);
	//used_reads_as_query[it->first] = true;
	current_reads[it->first] = true;
	break;
      }

      if(reads_in_current_cluster.empty()) break;

      for(size_t i=0;i<reads_in_current_cluster.size();i++)
      {
	map<string,vector<string> >::iterator it_ = seudoref_alignedreads_map.find(reads_in_current_cluster[i]);
	map<string,bool >::iterator it1 = used_reads_as_ref.find(reads_in_current_cluster[i]);
	map<string,bool >::iterator it2 = used_reads_as_query.find(reads_in_current_cluster[i]);

	if(it1 == used_reads_as_ref.end() && it_ != seudoref_alignedreads_map.end()) //as ref
	{ 
	    used_reads_as_ref[reads_in_current_cluster[i]] = true;
	    for(size_t j=0;j<it_->second.size();j++)
	    {
	     	if(current_reads.find(it_->second[j]) == current_reads.end())
	    	{
	            reads_in_current_cluster.push_back(it_->second[j]);
		    current_reads[it_->second[j]] = true;
	    	}
	    }
	}
	it_ = reads_seudoref_map.find(reads_in_current_cluster[i]);
	if(it2 == used_reads_as_query.end() && it_ != reads_seudoref_map.end()) //as reads
	{
	    used_reads_as_query[reads_in_current_cluster[i]] = true;
	    for(size_t j=0;j<it_->second.size();j++)
	    {
	        if(current_reads.find(it_->second[j]) == current_reads.end())
		{
		    reads_in_current_cluster.push_back(it_->second[j]);
		    current_reads[it_->second[j]] = true;
		}
	    }
	}

      }
      cout<<"Cluster_"<<idx<<" "<<reads_in_current_cluster.size()<<": ";
      for(size_t i=0;i<reads_in_current_cluster.size();i++) cout<<reads_in_current_cluster[i]<<" "<<endl;
      cout<<endl;
      break;
    }
}
void get_seudo_alignment_info(char* file,map<string,vector<string> >& seudoref_alignedreads_map)
{
    cerr<<"Begin getting sedu-alignment info from "<<file<<"..."<<endl;
    ifstream in(file);
    string s;
    istringstream istr;
    string currChr;
    vector<string> currReads;
    vector<string> currReads_flag;
    vector<int> currReads_left, currReads_right;
    getline(in,s);
    string id, chr,temp,flag;
    int start,end,s1,s2;//s1 s2: chaining score of the primary and best secondary alignment score
    int exon_num, minimizer_num;
    int M, NM, I,D,soft_l,soft_r;

    currChr = chr;
    int II = 0;
    while(getline(in,s))
    {
	istr.str(s);
	istr>>id>>chr>>temp>>start>>end>>flag
	    >>s1>>s2>>exon_num>>minimizer_num
	    >>temp>>M>>NM>>I>>D>>soft_l>>soft_r;
	    ;
	istr.clear();
	if( flag != "primary")
	{
	  istr.clear();
	  continue;
	}
	bool c1 = M > I && M > D && M > soft_l && M > soft_r;
	bool c2 = (M <= soft_l || M <= soft_r) && s2==0;
	if(!c1 && !c2)
	{
	    continue;
	}
	II++;
	map<string, vector<seudoref_info> >::iterator it = cluster_seudorefinfo_map.find(chr);
	if(it == cluster_seudorefinfo_map.end()) continue;
	for(size_t j=0;j<it->second.size();j++)
	{
	    if(start>=it->second[j].spos && end <= it->second[j].epos)
	    {
		
	        add_reads(seudoref_alignedreads_map,it->second[j].rid,id);
		/*
		cout<<it->second[j].rid<<"("<<it->second[j].spos<<" "<<it->second[j].epos<<") -> "
			<<id<<"("<<start<<" "<<end<<")"<<'\n';
		*/
	    }
	    else cout<<"BAD"<<endl;
	}

	istr.clear();
    }
    cout<<"check: "<<II<<endl;

}
map<string,string> seudoref_readid_map;
void get_reads_in_seudo_ref(char*file)
{
    cerr<<"Begin getting ref-reads info from "<<file<<"..."<<endl;
    ifstream in(file);
    string s;
    istringstream istr;
    while(getline(in,s))
    {
        string cidx,id,id_,temp;
	int pos1,pos2;
	istr.str(s);
	istr>>cidx>>id;
	while(istr>>temp){
	    if(temp == ":") break;
	}
	istr>>pos1>>pos2;
	istr.clear();

	istringstream istr1;
	istr1.str(id);
	istr1>>id_;
	id_ = id_.substr(1);
	reads_in_seudoref[id_] = true;
	seudoref_readid_map[cidx] = id_;
	//cout<<"--- "<<id_<<endl;
	seudoref_info sri={id_,pos1,pos2};

	//cout<<cidx<<" "<<id_<<" "<<pos1<<" "<<pos2<<endl;
	if(cluster_seudorefinfo_map.find(cidx) == cluster_seudorefinfo_map.end()){
	     vector<seudoref_info> vinfo(1,sri);
	     cluster_seudorefinfo_map[cidx] = vinfo;
	}
	else cluster_seudorefinfo_map[cidx].push_back(sri);
    }
    cout<<"seudo ref number: "<<cluster_seudorefinfo_map.size()<<" "<<reads_in_seudoref.size()<<endl;
    return;
}
void get_self_alignment_info(char* file)
{
     cerr<<"Begin getting sedu-alignment info from "<<file<<"..."<<endl;
    ifstream in(file);
    string s;
    istringstream istr;
    string currChr;
    vector<string> currReads;
    vector<string> currReads_flag;
    vector<int> currReads_left, currReads_right;
    getline(in,s);
    string id1,id2, chr,temp,flag;
    int start,end;
    double score;
    currChr = chr;
    while(getline(in,s))
    {
        istr.str(s);
        istr>>id1>>id2>>temp>>start>>end>>flag>>score;
	istr.clear();
	if( flag != "primary"){
	    //istr.clear();
	    //continue;
	}
	if(seudoref_readid_map.find(id1) != seudoref_readid_map.end() 
			&& seudoref_readid_map.find(id2) != seudoref_readid_map.end())
	{
	    string r1 = seudoref_readid_map[id1], r2 = seudoref_readid_map[id2]; //r1 mapped to r2
	    if(r1 == r2) continue;
	    int idx1 = dataID_idx_map[r1],idx2 = dataID_idx_map[r2];
	    int max=data[idx1].length();
	    if(max < idx2) max - data[idx2].length();

	    //if((1.0*(end-start + 1))/(1.0*max) <= 0.3) continue; //new1

	    if(seudoref_seudoreads_map.find(r2) == seudoref_seudoreads_map.end())//r2->r1 r1 mapped to r2
	    {
	        vector<string> v(1,r1);
		seudoref_seudoreads_map[r2] = v;
	    }
	    else seudoref_seudoreads_map[r2].push_back(r1);

	    if(seudoreads_seudoref_map.find(r1) == seudoreads_seudoref_map.end())//r1->r2 r1 mapped to r2
	    {
	        vector<string> v(1,r2);
		seudoreads_seudoref_map[r1] = v;
	    }
	    else seudoreads_seudoref_map[r1].push_back(r2);

	    if(seudoreads_seudoref_map2.find(r1) == seudoreads_seudoref_map2.end())//r1->r2 r1 mapped to r2
	    {
	        pair<string,double> p = make_pair(r2,score);
		seudoreads_seudoref_map2[r1] = p;
	    }
	    else
	    {
	        if(seudoreads_seudoref_map2[r1].second<score)
		{
		    pair<string,double> p = make_pair(r2,score);
		    seudoreads_seudoref_map2[r1] = p;
		}
	    }

	}
	

    }
    for(map<string,vector<string> > ::iterator it = seudoref_seudoreads_map.begin();it != seudoref_seudoreads_map.end();it++)//ref -> read unique
    {
        sort(it->second.begin(),it->second.end());
	it->second.erase(unique(it->second.begin(),it->second.end()),it->second.end());
    }
    cout<<"self-to-self size: "<<seudoref_seudoreads_map.size()<<" "<<seudoreads_seudoref_map.size()<<endl;
}
int main(int argc,char* argv[])
{
    if(argc == 1) {
        cout<<"-"<<endl;
        cout<<"This is a simple program to get the cluster information from seudo_read_alignment.info(id chr strand start_pos end_pos flag)."<<endl;
        cout<<"./exe reads.fastq fq (or fa for .fasta file) Nx_ref-read.info Nx_read_alignment.info Nx_self.alignment.info Nx (e.g. N10) "<<endl;
        cout<<"-"<<endl;
	return 0;
    }
     string filetype=argv[2];

    if(filetype == "fq")
        load_data_fastq(argv[1]);
    else if(filetype == "fa")
        load_data_fasta(argv[1]);
    sort(dataidx_length_vec.begin(),dataidx_length_vec.end(),sorter);
    Nx = argv[6];
    output_dir = argv[7];
    double x = atof(Nx.substr(1).c_str());
    N = (1.0*x)/100.0;

    get_reads_in_seudo_ref(argv[3]);
    get_seudo_alignment_info(argv[4],seudoref_alignedreads_map);
    get_self_alignment_info(argv[5]);

    string file1 = output_dir + "/" + Nx+".myCluster_singleton.tsv";
    string file2 = output_dir + "/" + Nx+".myCluster_multi.tsv";
    string file3 = output_dir + "/" + Nx+".singleton.fasta";
    string file4 = output_dir + "/" + "unconsidered.tsv";
    out_my_cluster.open(file1.c_str()); //single
    out_my_cluster2.open(file2.c_str()); //multi
    out_singleton.open(file3.c_str());//singleton_read
    if(Nx == "N90") 
	    out_post.open(file4.c_str());

    get_cluster3(seudoref_alignedreads_map);//get cluster from on rean_alignment.info(N50 or N30 or N10)
    out_my_cluster.close();
    out_my_cluster2.close();
    return 0;
}
