![logo](geluster.png)
Description
================

GeLuster is a highly efficient tool for clustering long-read transcriptomic data.


Prerequisites
================

  g++ with support for C++11 (e.g. 4.7.2)
 
  [minimap2][minimap] (In the GeLuster paper we tested with version 2.24-r1122)
  
  [samtools][samtools] (In the GeLuster paper we tested with version 1.10)

# 1. Installing and test
===========================================================================
    
    (A) Please ensure that minimap2 and samtools are properly installed. 
        
        Set the environment variables or type the commands when using GeLuster:

          $ export PATH=/path_to_minimap2/:$PATH
          $ export PATH=/path_to_samtools/:$PATH

        Please change path_to_minimap2(path_to_samtools) to the directory of minimap2(samtools).
     
        e.g. export PATH=/home/yuting/minimap2-2.24/:$PATH
        
    (B) Change to the GeLuster/src directory and make.
    
          $ cd src
          $ make release
          
    (C) Test GeLuster on the demo data set.
        
        Change to GeLuster/sample_test/, and type the following command:
        
          $ ./run_me.sh
          
        If you get the geluster_outdir/GeLuster.tsv, congrats, you succesfully installed GeLuster.
      
        
===========================================================================

# 2. Usage 
===========================================================================
    
    GeLuster v.1.1 usage:

    ** Required **
    
    --reads/-r <string>		: path to the read file
---------------------------------------------------------------------------

    ** Options **
    
    --help/-h               : Output GeLuster Help Information.

    --version/-v            : Print current version of GeLuster.

    --iteration/-i <int>    : Number of GeLuster iterations ([3,9], default: 3).

    --seqType/-s <string>   : 'cDNA' for ONT cDNA data, 'dRNA' for ONT direct RNA data, or 'PacBio' for pacbio data (default:cDNA).

    --rform/-f <string>     : 'fq' for FASTQ format reads, 'fa' for FASTA format reads (default: fq).
   
    --threads/-t <int>      : Number of threads to be used (default: 10).

    --multi/-M              : To generate a proxy of gene expression matrix for multiple RNA-seq samples. Input files should be separated by commas.

    --output_dir/-o <string>: Output path, default: geluster_outdir.

---------------------------------------------------------------------------

    ** Typical Command **
    
    A typical GeLuster command might be:

    GeLuster -r reads.fastq -f fq -s cDNA -o geluster_outdir

---------------------------------------------------------------------------

===========================================================================

# 3. Input and Output 
===========================================================================
  
  The GeLuster's input is long reads in FASTQ or FASTA format, e.g.,
    
    @SRR14181741.1 length=787
    GGTATTACTTCGTTCAGTTACGTATTGCTAAGGTTAACACAAAGACACCATTCTTTCTTCAGCA...
    +SRR14181741.1 length=787
    %%1+)+7<<3.)$'6(5'2,)/6,,/00203&9$-1-233JG209:3./$'&##$%&)(/,**+...
    @SRR14181741.2 length=597
    TGTTATGCACTTGTTCAGTTACGTATTACTAGGGTTAACACAAAGACACCATTTCATAACTTTCT...
    +SRR14181741.2 length=597
    )$#%$%($%$$$*/1422++/04:334($'#$(($%('(&&98+#('&$&'''')$#'$$%')++...

---------------------------------------------------------------------------  
  
  The output is a text file with one line per read.
  Each line specifies the name of the read and the cluster to which the read belongs, e.g.,
    
    @SRR14181741.319454 length=598 ,gene_cluster_0
    @SRR14181741.811316 length=399 ,gene_cluster_0
    @SRR14181741.929273 length=2216 ,gene_cluster_0
    @SRR14181741.1036973 length=470 ,gene_cluster_0
    @SRR14181741.354389 length=3828 ,gene_cluster_0
    @SRR14181741.831192 length=469 ,gene_cluster_0
    @SRR14181741.252153 length=317 ,gene_cluster_0
    @SRR14181741.100390 length=1934 ,gene_cluster_1
    @SRR14181741.1514937 length=783 ,gene_cluster_1
    @SRR14181741.722422 length=793 ,gene_cluster_1
    @SRR14181741.140922 length=2179 ,gene_cluster_1
    @SRR14181741.543032 length=1002 ,gene_cluster_1
    
    
===========================================================================
    
Authors: Ting Yu and Junchi Ma designed and wrote GeLuster.
 
Contact: Any questions, problems, bugs are welcome and should be dumped to Ting Yu <yutingsdu@163.com>
 
Updated on January 4, 2024.

 [minimap]: https://github.com/lh3/minimap2/releases
 [samtools]: https://github.com/samtools/samtools/releases
 