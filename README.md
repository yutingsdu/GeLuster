![logo](geluster.png)
Description
================

GeLuster is a fast and accurate tool for clustering long-read transcriptomic data.

It is free to use, modify, redistribute without any restrictions, except including the license provided with the distribution.


Installation
================

Prerequisites

 g++ with support for C++11 (e.g. 4.7.2)

# 1. Installing 
===========================================================================

===========================================================================

# 2. Usage 
===========================================================================

    GeLuster v.1.0 usage:

    ** Required **


---------------------------------------------------------------------------

    ** Options **
    --help/-h			  : Output GeLuster Help Information.

    --version/-v			  : Print current version of GeLuster.

    --iteration/-i <int>		  : Number of GeLuster iterations ([3,9], default: 3).

    --seqType/-s <string>		  : dRNA for direct RNA (default:cDNA).

    --rform/-f <string>		  : fa for .fasta(.fa) format, fq for .fastq(.fq) format (default: fq).

    --output_dir/-o <string>	  : Output path, default: geluster_outdir.


---------------------------------------------------------------------------

    ** Typical commands **
    
    A typical GeLuster command might be:

    GeLuster -r reads.fastq -f fq -s cDNA -o geluster_outdir



---------------------------------------------------------------------------

===========================================================================


Authors: Ting Yu and Junchi Ma designed and wrote GeLuster.

Contact:
 Any questions, problems, bugs are welcome and should be dumped to
 Ting Yu <yutingsdu@163.com>
 Created on July 15, 2023.

