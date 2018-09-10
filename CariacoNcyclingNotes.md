##Notes for Revisions on N Cycling Manuscript  


##**July 3rd 2018** 
### Run through Xander Tutorial from EDAMAME workshop so that we can use it on Cariaco dataset

Following Taylor's [tutorial](https://github.com/edamame-course/Xander/blob/master/Xander.md).  
Launched my own Amazon EC2 instance using EDAMAME community AMI so I will be charged.  
Public DNS: ec2-54-174-228-12.compute-1.amazonaws.com  
Using same amazon key as downloaded in workshop.  

**NOTE**- reference sequences for tutorial have been manually examined. When doing this for the paper, take a look at the alignment later (file: ref_aligned.faa).  

**NOTE2**- requires high quality, trimmed reads that have already been merged as input. Maria has already done this and will share tomorrow.

<u>Running through tutorial...</u>

First time got errors. Re-did edits to xander_setenv.sh file and it worked so probably typed something wrong. Didn't change any parameters.  

<u>Output files?</u>  
Tutorial explains some. This is also from the xander [website](https://github.com/rdpstaff/Xander_assembler):  

* *Output 1: contig coverage (coverage.txt, can be used to estimate gene abundance and adjust sequence abundance)*.  

* *Output: taxonomic abundance adjusted by coverage, group by lineage (phylum/class) (taxonabund.txt)*.  

AND: 

* *A script in bin/get_OTUabundance.sh is provided to create coverage-adjusted OTU abundance data matrix from contigs of same gene from multiple samples. The data matrix can then imported to R or PhyloSeq for more extensive analysis and visualization functions (see http://rdp.cme.msu.edu/tutorials/stats/RDPtutorial_statistics.html)*

	```
 * Input 1: aligned protein contig files (final_prot_aligned.fasta)
 * Input 2: contig coverage (coverage.txt)
 * Output: data matrix file with the OTU abundance at the each distance
```

Also [Rossmassler et al. 2016](https://academic.oup.com/femsle/article/363/16/fnw162/2197704) say they use this script for their read counts to the contigs. See their table 2.  

<u>Need to start meandering from tutorial a bit.  </u> 

First, copy the get_OTUabundance.sh script into the working directory

```
cp /home/ubuntu/tools/RDPTools/Xander_assembler/bin/get_OTUabundance.sh /home/ubuntu/tools/RDPTools/Xander_assembler/rplB_demo
```

[Website for script](https://github.com/rdpstaff/Xander_assembler/blob/master/bin/get_OTUabundance.sh) says to modfy line 17 but I don't think I have to. It seems like it directs to the right place.  
Run script and specify input/ output according to documentation. These are required parameters:  

```
coverage_file outdir start_dist end_dist aligned_files 

## Input 1: a contig coverage file (used to adjust the sequence abundance)
## Input 2: output directory
## Input 2: start and end distance cutoff
## Input 4: takes the aligned protein contig files (_final_prot_aligned.fasta), must be from the same gene 
## Output: one data matrix containing the OTU abundance at the each distance cutoff
```

Make a directory for results. Then run with input according to above

```
mkdir OTUabundanceresults

./get_OTUabundance.sh k45/rplB/cluster/test_rplB_45_coverage.txt  OTUabundanceresults 0 0.1 k45/rplB/cluster/test_rplB_45_final_prot_aligned.fasta
```

Ok... this make a bunch of files called

```
-rw-rw-r-- 1 ubuntu ubuntu 2.2K Jul  3 18:38 rformat_dist_0.01.txt
-rw-rw-r-- 1 ubuntu ubuntu 1.7K Jul  3 18:38 rformat_dist_0.02.txt
-rw-rw-r-- 1 ubuntu ubuntu 1.4K Jul  3 18:38 rformat_dist_0.03.txt
-rw-rw-r-- 1 ubuntu ubuntu 1.2K Jul  3 18:38 rformat_dist_0.04.txt
-rw-rw-r-- 1 ubuntu ubuntu  824 Jul  3 18:38 rformat_dist_0.05.txt
-rw-rw-r-- 1 ubuntu ubuntu  735 Jul  3 18:38 rformat_dist_0.06.txt
-rw-rw-r-- 1 ubuntu ubuntu  629 Jul  3 18:38 rformat_dist_0.07.txt
-rw-rw-r-- 1 ubuntu ubuntu  548 Jul  3 18:38 rformat_dist_0.08.txt
-rw-rw-r-- 1 ubuntu ubuntu  513 Jul  3 18:38 rformat_dist_0.09.txt
-rw-rw-r-- 1 ubuntu ubuntu 2.4K Jul  3 18:38 rformat_dist_0.0.txt
-rw-rw-r-- 1 ubuntu ubuntu  450 Jul  3 18:38 rformat_dist_0.1.txt
```

Each of this is an OTU table with counts assigned to the different OTUs. Try to understand this.  



##**July 5th 2018**   

### Starting to work with Cariaco Data



- Maria shared metatranscriptomes from November. F and R reads have already been merged and quality controlled.

- I will run final analyses on FULL dataset (metaT and metaG) incase something is in low abundance in metaT, the full contig should be built using the metaG?

- Meanwhile, just downloading and moving datasets is taking a long time. I am making a local copy of the Nov metaT samples (19 fasta files) that is backed up in Wagner's Google drive (unlimited storage space)

- While that is transferring, tried to understand the getOTUabudance.sh function and its output in more detail
	- From git page: *This script clusters the aligned protein contigs from multiple samples and creates a data matrix file with the OTU abundance at the each distance cutoff*
	- From Rossmassler et al.: *Operational taxonomic units (OTUs) for each gene were generated by the get OTUabundance.sh script in Xander using a distance of 0.1 (Wang et al. 2015).*
	- So from this I think the script is generating OTUs at different distance cutoffs (0.01, 0.02, etc. until the maximum I set at 0.1) and generated an OTU table for each.
	- The OTU table as sample ID in rows (in the case of the demo, only 1) and OTUs as columns
	- OTUs are numbered and are taken from the aligned protein contig files (_final_prot_aligned.fasta)
	- However, I can't figure out which contigs are clustered into each of the defined OTUs. They are just labeled (1, 2, 3, etc)
		- I wrote to the RDP staff asking them this question.
		- For now there is still some very useful info in this. For example, before running the ```get_OTUabundance.sh function``` [in just the output from ```xander_setenv.sh``` ] there are several files:
			- ```final_prot.fasta and final_prot_aligned.fasta``` are quality filtered (and aligned) protein contigs
			- ```framebot.txt``` is the the nearest reference seq and % aa identity for each contig
			- ```coverage.txt``` is the contig coverage (from each sample?) and can be used to estimate gene abundance and adjust sequence abundance (they say... but I think that info actually comes from ```get_OTUabundance.sh```)
			- ```taxonabund.txt``` is the taxonomic abundance (taken from ```framebot.txt``` adjusted by coverage, and grouped by lineage (phylum/class)
				-Check out the column "Fraction Abundance" in this output file. I think this is the most useful data so far- this is the relative abundance of the reads matching each of those contigs. 


- Still waiting for data upload
- Meanwhile make directory for practice on Cariaco metaT
	- In ```~/tools/RDPTools/Xander_assembler ``` make new directory ```mkdir nirK_Cariaco_practice```
	- Copy functions into this folder

		```
cp /home/ubuntu/tools/RDPTools/Xander_assembler/bin/xander_setenv.sh /home/ubuntu/tools/RDPTools/Xander_assembler/nirK_Cariaco_practice 
cp /home/ubuntu/tools/RDPTools/Xander_assembler/bin/run_xander_skel.sh /home/ubuntu/tools/RDPTools/Xander_assembler/nirK_Cariaco_practice 
cp /home/ubuntu/tools/RDPTools/Xander_assembler/bin/get_OTUabundance.sh /home/ubuntu/tools/RDPTools/Xander_assembler/nirK_Cariaco_practice
```
- Also copying fasta files from Maria. Just used Filezilla to drag them into ```nirK_Cariaco_practice``` folder. This takes a while. Starting with one file for now (R2b267A_merged.fa).
	- Try running through pipeline to see if I can do the analysis for nirK for this one sample
	- Use the nirK references already in the tutorial (this may need to be updated in the real analysis- they are from 2015?)
	- First modify the ```xander_setenv.sh``` script. Nano and replace file paths with the following

	```
SEQFILE=home/ubuntu/Cariaco_data/R2b267A_merged.fa    
WORKDIR=/home/ubuntu/tools/RDPTools/Xander_assembler/nirK_Cariaco_practice
REF_DIR=/home/ubuntu/tools/RDPTools/Xander_assembler
JAR_DIR=/home/ubuntu/tools/RDPTools
UCHIME=/home/ubuntu/tools/third_party_tools/uchime4.2.40_i86linux32
HMMALIGN=/usr/local/bin/hmmalign
```

- Change security of the file ```chmod 755 xander_setenv.sh```

- Ready to run but first make a tmux session

```
tmux new -s Cariaco_nirKpractice 
./run_xander_skel.sh xander_setenv.sh "build find search" "nirK nirS rplB"
```
- Started at 5:35pm (21:35 on ubuntu's clock). Let's see how long this takes...
- Failed:```search contigs failed for nirK```
- Tried running nirS and rplB to see if I can get any to work. Start at 6:01 pm
	```./run_xander_skel.sh xander_setenv.sh "build find search" "nirK nirS rplB"
```



  
Didn't work. Tried for nirK first, then nirK, nirS and rplB. This is error output:  

```
File "/home/ubuntu/tools/RDPTools/Xander_assembler/pythonscripts/getUniqueStarts.py", line 37, in <module>

    getUnique(infile)

  File "/home/ubuntu/tools/RDPTools/Xander_assembler/pythonscripts/getUniqueStarts.py", line 27, in getUnique

    kmer_pos = lexems[3] + "_" + lexems[7]

IndexError: list index out of range

get uniq starting kmers failed for nirS

get uniq starting kmers failed for rplB

### Search contigs nirK

java -Xmx2G -jar /home/ubuntu/tools/RDPTools/hmmgs.jar search -p 20 1 100 ../k45.bloom /home/ubuntu/tools/RDPTools/Xander_assembler/gene_resource/nirK/for_enone.hmm /home/ubuntu/tools/RDPTools/Xander_assembler/gene_resource/nirK/rev_enone.hmm gene_starts.txt 1> stdout.txt 2> stdlog.txt

search contigs failed for nirK

```

##**July 9th 2018** 
### Try again using metaG file, which will hopefully have more coverage and allow for building of contigs 

- Removed metaT file "R2b267A_merged.fa"
- Added metaG file in same place, "D2b267A_merged.fa.gz"
- Uploaded as zipped file (.gz)
- Therefore need to unzip with gunzip ```gunzip D2b267A_merged.fa.gz```
- Already copied the functions into my nirK_practice folder. Just need to edit xander_setenv so it cd's to correct directory

```
SEQFILE=/home/ubuntu/Cariaco_data/D2b267A_merged.fa
WORKDIR=/home/ubuntu/tools/RDPTools/Xander_assembler/nirK_Cariaco_practice
REF_DIR=/home/ubuntu/tools/RDPTools/Xander_assembler
JAR_DIR=/home/ubuntu/tools/RDPTools
UCHIME=/home/ubuntu/tools/third_party_tools/uchime4.2.40_i86linux32
HMMALIGN=/usr/local/bin/hmmalign

```

- And change permissions of file ```chmod 755 xander_setenv.sh```
- And run Xander. Only use nirS for the first try ```./run_xander_skel.sh xander_setenv.sh "build find search" "nirS"```
- I think it worked! This is output in terminal

```
### Build bloom filter
java -Xmx2G -jar /home/ubuntu/tools/RDPTools/hmmgs.jar build /home/ubuntu/Cariaco_data/D2b267A_merged.fa k45.bloom 45 32 2 4 30 >& k45_bloom_stat.txt
### Find starting kmers for nirS
java -Xmx2G -jar /home/ubuntu/tools/RDPTools/KmerFilter.jar fast_kmer_filter -a -o temp_starts_1.txt -t 1 45 /home/ubuntu/Cariaco_data/D2b267A_merged.fa nirS=/home/ubuntu/tools/RDPTools/Xander_assembler/gene_resource/nirS/ref_aligned.faa 
Starting kmer mapping at Mon Jul 09 14:10:28 UTC 2018
*  Number of threads:       1
*  References:              [nirS]
*  Reads file:              /home/ubuntu/Cariaco_data/D2b267A_merged.fa
*  Kmer length:             15
*  Kmer Refset Size:        40774
Processed 1000001 sequences in 75911 ms
Processed 1999999 sequences in 150785 ms
Processed 3000000 sequences in 232755 ms
Processed 4000000 sequences in 312972 ms
Processed 5000001 sequences in 387762 ms
Finished Processed 5827663 sequences in 451411 ms
### Search contigs nirS
java -Xmx2G -jar /home/ubuntu/tools/RDPTools/hmmgs.jar search -p 20 1 100 ../k45.bloom /home/ubuntu/tools/RDPTools/Xander_assembler/gene_resource/nirS/for_enone.hmm /home/ubuntu/tools/RDPTools/Xander_assembler/gene_resource/nirS/rev_enone.hmm gene_starts.txt 1> stdout.txt 2> stdlog.txt
### Merge contigs
java -Xmx2G -jar /home/ubuntu/tools/RDPTools/hmmgs.jar merge -a -o merge_stdout.txt -s test -b 50 --min-length 150 /home/ubuntu/tools/RDPTools/Xander_assembler/gene_resource/nirS/for_enone.hmm stdout.txt gene_starts.txt_nucl.fasta
Read in 321 contigs, wrote out 138 merged contigs in 1.208s
Warning, no derep mode specified, using unaligned
Processing prot_merged.fasta
Total sequences: 138
Unique sequences: 12
Dereplication complete: 169
### Cluster
Merging file alignment/aligned.stk
Using #=GC_RF as mask sequence
Processing aligned.fasta
Total sequences: 12
Unique sequences: 12
Dereplication complete: 92
Reading sequences(memratio=0.0010564018821559154)...
Using distance model edu.msu.cme.rdp.alignment.pairwise.rna.IdentityDistanceModel
Read 12 Protein sequences (memratio=0.0014101074017590462)
Reading ID Mapping from file /home/ubuntu/tools/RDPTools/Xander_assembler/nirK_Cariaco_practice/k45/nirS/cluster/ids
Read mapping for 12 sequences (memratio=0.0014101074017590462)
Starting distance computations, predicted max edges=144, at=Mon Jul 09 14:18:21 UTC 2018
Dumping 22 edges to partial_matrix0 FINAL EDGES (memory ratio=0.0010144265920570413)
Matrix edges computed: 26
Maximum distance: 0.625
Splits: 1
Partition files merged: 13
Doing complete linkage clustering with step 0.009999999776482582 (realstep=100)
Clustering complete: 20
Converted 12 sequences from [complete.clust_rep_seqs.fasta] (FASTA) to fasta in 0 s
### Chimera removal
uchime v4.2.40
by Robert C. Edgar
http://drive5.com/uchime
This code is donated to the public domain.

00:00  18Mb  100.0% Reading test_nirS_45_nucl_rep_seqs.fasta
00:00  18Mb 12 sequences                                    
00:00  19Mb  100.0% Reading /home/ubuntu/tools/RDPTools/Xander_assembler/gene_resource/nirS/originaldata/nucl.fa
00:00  19Mb 233 sequences                                                                                       
00:00 8.0Mb  100.0% 0/11 chimeras found (0.0%)
### FrameBot
java -jar /home/ubuntu/tools/RDPTools/FrameBot.jar framebot -N -l 150 -o nirS_45 /home/ubuntu/tools/RDPTools/Xander_assembler/gene_resource/nirS/originaldata/framebot.fa nucl_rep_seqs_rmchimera.fasta
### Kmer abundance
java -Xmx2g -jar /home/ubuntu/tools/RDPTools/KmerFilter.jar kmer_coverage -t 1 -m test_nirS_45_match_reads.fa 45 test_nirS_45_final_nucl.fasta test_nirS_45_coverage.txt test_nirS_45_abundance.txt /home/ubuntu/Cariaco_data/D2b267A_merged.fa
```

- Output files are in "cluster" folder. Downloaded some results files and pasted in excel file
	-  Check out CompleteClust. Shows different distance cutoffs for grouping contigs. Using distance of 0, there are 12 unique contigs. Using distance of 0.1, there are 7. This file shows which contigs are in which cluster.
		- This probably matches the "get_OTUabundance.sh" output, which then counts the reads from each sample which belong in each OTU (or cluster). Using these two files together, I could pull out the identity (based on contigs in each cluster) and read counts within that identity for each sample. 

Try again for nirK  
``` ./run_xander_skel.sh xander_setenv.sh "build find search" "nirK" ```

- NOTE this takes less time because it already built all contigs (data in folder "k45" for kmer 45). Now it is just scanning them for nirK
- Not able to find nirK contigs in this sample. Output is:

```
File k45.bloom exists, SKIPPING build (manually delete if you want to rerun)
### Find starting kmers for nirK
java -Xmx2G -jar /home/ubuntu/tools/RDPTools/KmerFilter.jar fast_kmer_filter -a -o temp_starts_1.txt -t 1 45 /home/ubuntu/Cariaco_data/D2b267A_merged.fa nirK=/home/ubuntu/tools/RDPTools/Xander_assembler/gene_resource/nirK/ref_aligned.faa 
Starting kmer mapping at Mon Jul 09 14:44:35 UTC 2018
*  Number of threads:       1
*  References:              [nirK]
*  Reads file:              /home/ubuntu/Cariaco_data/D2b267A_merged.fa
*  Kmer length:             15
*  Kmer Refset Size:        30824
Processed 1000001 sequences in 78263 ms
Processed 2000000 sequences in 153445 ms
Processed 2999999 sequences in 236171 ms
Processed 4000000 sequences in 319531 ms
Processed 5000001 sequences in 395038 ms
Finished Processed 5827663 sequences in 456326 ms
### Search contigs nirK
java -Xmx2G -jar /home/ubuntu/tools/RDPTools/hmmgs.jar search -p 20 1 100 ../k45.bloom /home/ubuntu/tools/RDPTools/Xander_assembler/gene_resource/nirK/for_enone.hmm /home/ubuntu/tools/RDPTools/Xander_assembler/gene_resource/nirK/rev_enone.hmm gene_starts.txt 1> stdout.txt 2> stdlog.txt
### Merge contigs
java -Xmx2G -jar /home/ubuntu/tools/RDPTools/hmmgs.jar merge -a -o merge_stdout.txt -s test -b 50 --min-length 150 /home/ubuntu/tools/RDPTools/Xander_assembler/gene_resource/nirK/for_enone.hmm stdout.txt gene_starts.txt_nucl.fasta
Read in 2 contigs, wrote out 0 merged contigs in 0.188s
```


  
  
  
Try again for amoA_AOA

- Worked! Only assembled into 1 contig with 31 reads. Save output files in excel sheet
- The contig is a Thaumarchaeota Marine Group I!!!!



###I need to investigate how to scale this up.  

- Space on this instance is limited. I tried to add metaG file without removing metaT file. Got this error: 

	```
~/Cariaco_data$ mkdir metaT
mkdir: cannot create directory ‘metaT’: No space left on device
```

- Do I need my own instance (volume?) to which I can upload all the Cariaco files, download the xander algorithms and RDP gene resources, and run virtually. That is going to cost a bit more money...  
- Wrote to Taylor to ask.
	- Taylor linked an ebooklet she is writing [here](https://github.com/dunivint/RDP_Tutorials).
	- Also this is text of her email

```
Hi Liz!
It is good to hear from you, and I am excited to hear you are interested in using
 Xander!! 

It is absolutely possible to re-launch the EDAMAME-xander instance with a larger
 computer, so it is not necessary to build your own AMI. Since you have lots of
  files and will soon have lots of output files, I recommend looking into
   different storage options on amazon, which are explained here and here. With
    Amazon storage and amazon instances, you can "curl and unzip files" in an
     instance OR "attach volumes" when launching instances depending on the
      storage you select.  

The above advice is assuming that your home institution does not have a high
 performance computer (HPC). If it does, I recommend using this instead since it
  will likely be cheaper in the long run. Additionally, Xander comes with 
  RDPTools, which most university HPCs have. 

I'll also point you to a more detailed tutorial that I have been working on with 
the RDP. Eventually, this will be an e-boolket of sorts, so if you run into any 
issues with it or have comments, I would love to hear it! https://github.com/
dunivint/RDP_Tutorials

Let me know if you have any further questions with merging results from different
 samples, selecting parameters, etc!

Best, 
Taylor
```

- I read through the explanation of different storage options. I think Amazon EBS is best for my purposes.
	- Cost is $0.10 per GB-month of provisioned storage. Read this example:
		- For example, let's say that you provision a 2000 GB volume for 12 hours (43,200 seconds) in a 30 day month. In a region that charges $0.10 per GB-month, you would be charged $3.33 for the volume ($0.10 per GB-month * 2000 GB * 43,200 seconds / (86,400 seconds/day * 30 day-month)).
	- Need to make sure to do 2 things
		- Pick an EBS-enabled instance and configure it when launching (follow [these instructions](https://docs.aws.amazon.com/AWSEC2/latest/UserGuide/EBSOptimized.html)).
		- Also this advice: If you are creating a volume for a high-performance storage scenario, you should make sure to use a Provisioned IOPS SSD (io1) volume and attach it to an instance with enough bandwidth to support your application, such as an EBS-optimized instance or an instance with 10-Gigabit network connectivity. The same advice holds for Throughput Optimized HDD (st1) and Cold HDD (sc1) volumes. For more information, see Amazon EC2 Instance Configuration.

### Creating Amazon EBS volume

Steps from [here](https://docs.aws.amazon.com/AWSEC2/latest/UserGuide/ebs-creating-volume.html).	 

- Open the Amazon EC2 console at https://console.aws.amazon.com/ec2/.
- From the navigation bar, select the region in which you would like to create your volume. This choice is important because some Amazon EC2 resources can be shared between regions, while others can't. For more information, see Resource Locations.
- In the navigation pane, choose ELASTIC BLOCK STORE, Volumes.
- Choose Create Volume.
- For Volume Type, choose a volume type. For more information, see [Amazon EBS Volume Types](https://docs.aws.amazon.com/AWSEC2/latest/UserGuide/EBSVolumeTypes.html).
	- Went with General Purpose SSD (gp2) because it seems most standard
- For Size (GiB), type the size of the volume.
	- Chose 1000GiB. This is about 1TB (I think).
- With a Provisioned IOPS SSD volume, for IOPS, type the maximum number of input/output operations per second (IOPS) that the volume should support.
- For Availability Zone, choose the Availability Zone in which to create the volume. EBS volumes can only be attached to EC2 instances within the same Availability Zone.
	- Edamame instance is in us-east-1d


OK. I think it worked. Says Volume created. Volume ID is vol-0c81dbc7cf75e7b97


### Now attach EBS volume to Instance
Instructions [here](https://docs.aws.amazon.com/AWSEC2/latest/UserGuide/ebs-attaching-volume.html).

Searched for EDAMAME instance, which has ID i-0add16999f132ea0e.
Found it and attached. Device is ```/dev/sdf```.

- Also is telling me Linux Devices: /dev/sdf through /dev/sdp.  Note: Newer Linux kernels may rename your devices to /dev/xvdf through /dev/xvdp internally, even when the device name entered here (and shown in the details) is /dev/sdf through /dev/sdp.
- Next sign into instance in Terminal and Filezilla
- Use the ```lsblk``` command to view your available disk devices and their mount points (if applicable) to help you determine the correct device name to use.

```
NAME    MAJ:MIN RM   SIZE RO TYPE MOUNTPOINT
xvda    202:0    0     8G  0 disk 
└─xvda1 202:1    0     8G  0 part /
xvdf    202:80   0  1000G  0 disk 
```
It's there!

- The output of lsblk removes the /dev/ prefix from full device paths. 
- New volumes are raw block devices, and you must create a file system on them before you can mount and use them. Volumes that have been restored from snapshots likely have a file system on them already; if you create a new file system on top of an existing file system, the operation overwrites your data. Use the sudo file -s device command to list special information, such as file system type.

```
sudo file -s /dev/xvdf
/dev/xvdf: data
```

- The ouput "data" means there is no file system.
- Use the following command to create an ext4 file system on the volume. Substitute the device name (such as /dev/xvdf) for device_name. Depending on the requirements of your application or the limitations of your operating system, you can choose a different file system type, such as ext3 or XFS.

```
sudo mkfs -t ext4 /dev/xvdf
```

- Use the following command to create a mount point directory for the volume. The mount point is where the volume is located in the file system tree and where you read and write files to after you mount the volume. Substitute a location for mount_point, such as /data.

```
sudo mkdir Cariaco
```

- Use the following command to mount the volume at the location you just created.
	- Need to do this every time turn instance on/ off

```
sudo mount /dev/xvdf Cariaco
```

- Review the file permissions of your new volume mount to make sure that your users and applications can write to the volume. For more information about file permissions, see File security at The Linux Documentation Project.

```
mv Cariaco_data/D2b267A_merged.fa Cariaco
mv: cannot create regular file ‘Cariaco/D2b267A_merged.fa’: Permission denied
```

- I needed to use sudo to do this

```
sudo mkdir metaG
sudo mv Cariaco_data/D2b267A_merged.fa Cariaco/merged_seqs/metaG
```


- I think I have to re-mount every time I log on. Move a file into here and then log off/ on to see if it is still there.
	- NOTE you don't stop the EBS volume like you stop an instance. As long as you are storing stuff on there, you are paying for it. Check out this [link](https://aws.amazon.com/premiumsupport/knowledge-center/ebs-charge-stopped-instance/) about making the EBS into a snapshot once the info is there so the charges will be minimized.
	- Started up instance. Used ```lsblk``` to see my volumes and xvdf is there with 1000G. 
	- Try mounting again ```sudo mount /dev/xvdf Cariaco```
	- Data file is there.
	- Run ```df -t``` to see disc space
		- Mounted disc is ~1% in use
	- Try adding more through Filezilla
		- Can't because don't have permission to modify these folders in mounted disc. I really need to figure that out
		- Found this on the internet. Seems to work to give me permission from terminal without using sudo ```sudo chmod 777 Cariaco```
		- But still can't make directories or files or upload data in Filezilla
		-  Also tried this from internet ```sudo chown ubuntu /dev/xvdf```
		-  OK finally got it- I had to give permission at every folder level
```
sudo chmod 777 merged_seqs
sudo chmod 777 metaG
sudo chmod 777 metaT
```

- OK so I uploaded another metaG file into metaG folder. Used gunzip to unzip.
- Worked! I am still only using 1% of space. This volume is too big and I will end up paying too much for what I need.
- ###Next steps:###
	- Delete this volume (now that I kind of know how to set it up) and make one for 100GiB
		- Upload files and make it a snapshot to save $$$
	- Temrinate the EDAMAME instance and start new one with fresh files..
	- See if I can run Xander from files in volume
	- Add more gene resources (dsr, etc)



##**July 10th 2018** 
### Build new Instance from EDAMAME RDP instance and attach storage space with Amazon EBS

1. Chose community isntance, RDP-Edamame-2015, ID: ami-e973b782.  
	a. Launched with m4.xlarge (4 CPUS, 16GB RAM). Instance storage: EBS only.  
	b. using region: us-east-1c

2. 	Follow same instructions as yesterday for creating and attaching volume but make 100GiB instead of 1000.
	a. Volume ID: vol-065535bdd665fdc90

3. ssh into instance. Also Filezilla

4. Attached volume to instance in EC2 console. Device name is /dev/sdf.  
	a. However when I ```lsblk``` in terminal, it shows me a device named xvdf that has 100G

5. Create file system on mounted device.  
	a. ```sudo file -s /dev/xvdf```  
		i. output is ```/dev/xvdf: data```  
	b. create an ext4 file system on the volume
		i. ```sudo mkfs -t ext4 /dev/xvdf```  
	c. Create a mount point directory on the volume.  
		i. ```sudo mkdir Cariaco```  
	d. mount the volume at the mount point location.  
		i. ```sudo mount /dev/xvdf Cariaco```  
	e. Give permission to modify files in Cariaco for all users (so I can work from Filezilla)  
		i. ```sudo chmod 777 Cariaco```
		ii. Check permisions with ```ls -lah```  
			```drwxrwxrwx 3 root   root   4.0K Jul 10 13:05 Cariaco```  
	f. Run ```df -t``` to see disc space and location of Cariaco  
	
	```
	Filesystem     1K-blocks    Used Available Use% Mounted on
/dev/xvda1       8115168 2918068   4761824  38% /
none                   4       0         4   0% /sys/fs/cgroup
udev             8211616      12   8211604   1% /dev
tmpfs            1643316     336   1642980   1% /run
none                5120       0      5120   0% /run/lock
none             8216572       0   8216572   0% /run/shm
none              102400       0    102400   0% /run/user
/dev/xvdf      103081248   61044  97760940   1% /home/ubuntu/Cariaco
```
	g. Make directories within Cariaco	
	
	```
	cd Cariaco/
	mkdir metaG
	mkdir metaT
	```	

6. Everything with permissions is working fine. Start uploading files into metaG and metaT using Filezilla. (First download locally from my GoogleDrive folder: CloudStorageOnly. Try zipping folders so this is faster)
	a. Needed to update system and install zip. Then was able to start unzipping after file transfers.
	
	```
	sudo apt-get update
	sudo apt-get install unzip
	sudo unzip May_fasta-20180710T132044Z-002.zip 
	sudo unzip May_fasta-20180710T132044Z-008.zip 
	```

	- ```May_fasta-20180710T132044Z-002.zip``` contained ```R2b295A_merged.fa.gz``` and ```R3b103B_merged.fa.gz```
	- ```May_fasta-20180710T132044Z-008.zip ``` contained ```R2a103A_merged.fa.gz``` and ```R2b295A_merged.fa.gz```
	- Now up to 8% of storage space on EBS
	- While waiting for more files to download, I unzipped the actual sequence files in the instance in order to check them ```sudo gunzip * ```
		- Takes a while
		- **Now I'm up to 20% of my storage!!!! I was right in the beginning- I probably need 1000GiB not 100GiB!!**
			- Modifed storage space in EC2 console using [this tutorial](https://docs.aws.amazon.com/AWSEC2/latest/UserGuide/ebs-modify-volume.html?icmpid=docs_ec2_console)
			-  Took a few hours in the *Modifying state*. Then needed to do third step from , extend the file system. 
				-  First I stopped and restarted Instance then re mounted the EBS ```sudo mount /dev/xvdf Cariaco```
				-  Ran command ```df -h``` which shows the size of the mounted directory is ~100G while the ```lsblk``` command shows the disc space is actually 1000G
				-  Expand the disc using growpart ```sudo resize2fs /dev/xvdf```
				-  Check ```lsblk``` and ```df -h``` again to see that they both show 1000G in /home/ubuntu/Cariaco
				-  Now my volume is only 3% full!!
			



7. While waiting for that, start curating files on instance to suit my needs.  

- Need most up to date gene resources and gene resource folders that are not already in instance. 

```
cd tools/RDPTools/Xander_assembler/gene_resource/
```
 - These folders were last modified Jun 22 2015. 

- Check [RDP fungene database](http://fungene.cme.msu.edu/) and Taylor's [tutorial](https://github.com/edamame-course/Xander/blob/master/Xander.md) with some insight on this process.
	- According to tutorial, this is explanation of the files:
		- **gene.seeds:** a small set of full length, high quality protein sequences in FASTA format, used to build gene.hmm, forward and reverse HMMs. Can be downloaded from RDP's FunGene database.
		- **gene.hmm**: this is the HMM built from gene.seeds using original HMMER3. This is used to build for_enone.hmm and align contigs after assembly. Can be downloaded from RDP's FunGene database.
		- **framebot.fa**: a large near full length known protein set for identifying start kmers and FrameBot nearest matching. More diversity is better, more sequences means more starting points (more computational time) but less susceptible to noise than model creation. Prefer near full-length and well-annotated sequences. Filter with Minimum HMM Coverage at least 80 (%).
		- **nucl.fa**: a large near full length known set used by UCHIME chimera check.
	- To download the gene.seeds and gene.hmm files, click on the links at the top left of a page. 
	- To obtain the framebot.fa and nucl.fa files from RDP, click on the show/hide filter options link at the top right of the page. Here you can limit the the score, minimum size aa, and minimum HMM coverage.
		- The score will vary by gene. You may choose to leave this area blank. You can also manually look though the sequences to see if there is a large drop off in score at a specific number and use this number as a cutoff.
		- The minimum size aa will depend on the actual protein size. You'll want to find a balance between near-full length sequences and diversity.
		- The minimum HMM coverage must be greater than 80%. What percent cutoff you choose will require biological insight. Lower % HMM coverage will increase diversity but may lower the quality of your search. 
	
	
	- Once you have set your parameters, click filter. Then you will need to select sequences. It is not recommended to blindly select all sequences. Instead, manually go through and make sure no oddballs are included. For reference, the Xander paper used over 700 near full-length sequences for their .fa files. Once your sequences are selected, click begin analysis at the top right of the page. This will take you to a page where you can download your framebot.fa (protein) and nucl.fa (nucleotide) sequences. Then change the names to framebot.fa and nucl.fa respectively. Move all of these files to originaldata.

##Build new gene_resource for nirS##

- downloaded nirS.seeds and nirS. hmm files
	- Set minimum HMM coverage to 80%
		- Reduced number of sequences from 15,523 to 2,437
	- Sorted by size. Ranges from 798 to 493
		- Highest size (798) is environmental Chloroflexi strain with low score (726) compared to others.
		- Most of the models in the nirS.seeds file have a size of 596. Minimum length is 551
		- Sorted remaining sequences by size. from max size down to size of 494, there is a lot of diversity, including many environmental sequences. Starting at size 493, there is mostly Pseudomonas aeruginosa strains (except for last couple of sequences, including a potentially interesting Deltaproteobacteria environmental strain)
		- Filtered with minimum size 494. Remaining = 1,527 sequences.
	- Sorted by HMM_coverage again. Best = 99.8%, worst = 88%
		- There are a few things at the lower % coverage of this list that seem relevant to Cariaco. (eg. Dechloromonas aromatica RCB, Azoarcus sp. CIB, Rhodobacter sp. SW2) so I think keeping the HMM_coverage as low as 88% is OK.
	- Sorted by score. Best = 1097.2, worst = 59.2
		- There is a huge drop off from 531 to 181.9. I will cut off anything with score less than 530. This will get rid of some interesting environmental sequences but score is just too low :(
			- Reduce number of sequences from 1,527 to 1,251
	- **Final parameters** = Minimum Score: 530,  Minimum Size aa: 494,  Minimum HMM Coverage: 88
	- Now go through all 1,251 sequences one-by-one to remove "odd balls"
		- This is largely Pseudomonas strains
		- There are some relevant ones (Eg Cupriavidus, Marinobacter, Thioclava, Thiobacillus denitrificans, Azoarcus, Rhodobacteraceae, marine metagenome hypothetical proteins, uncultured Gammaprotebacteria, endosymbiont of Ridgeia piscesae, endosymbiont of Riftia pachyptila , Colwellia, Anaerolinea)
		- I don't see anything to remove. Checked a few using blastp to make sure they made sense (KFL36640: Arenimonas donghaensis DSM 18148 = HO3-R19 hypothetical protein and BAO30669: Sulfuritalea hydrogenivorans sk43H hydroxylamine reductase, EEW26927: Rhodobacter sp. SW2 Hydroxylamine reductase, CRP32869: Pseudomonas aeruginosa Nitrite reductase precursor).
- Selected all remaining sequences and, clicked "begin analysis" at the top right of the page. Downloaded protein and nucleotide fastas. 
	- Chose accno, fasta, not aligned [to match what I see in the tutorial gene resources)
	- Change name to framebot.fa (protein) and nucl.fa (nucleotide) respectively. Navigate to ```nirs/originaldata```. Move all 4 files to originaldata.
		- Check to make sure the files contain 1251 sequences

		```
grep ">" nucl.fa | wc
grep ">" framebot.fa | wc
```

- From these 4 files, we need to make three new files
	- ```for_enone.hmm``` and ```rev_enone.hmm``` for the forward and reverse HMMs respectively. This is used to assemble gene contigs.
	- ```ref_aligned.faa``` file containing a set of protein reference sequences aligned with for_`enone.hmm. This is used to identify starting kmers. 
	- Then just followed tutorial with following modification to prepare_gene script:

	```
	cd /home/ubuntu/tools/RDPTools/Xander_assembler/bin
nano prepare_gene_ref.sh
	```	
	- line 4: change gene to your gene of interest
	- line 9: change the jar directory to /home/ubuntu/tools/RDPTools
	- line 10: change the reference directory to /home/ubuntu/tools/RDPTools/Xander_assembler
	- Run ```./prepare_gene_ref.sh nirS```
		
	- Next, according to tutorial, need to manually examine the alignment in ```ref_aligned.faa```using Jalview or alignment viewing tool to spot any badly aligned sequences. 
		- "If found, it is likely there are no sequences in gene.seeds. Close these sequences. You need to validate these problem sequences to see if they are from the gene you are interested in, then either remove them or add some representative sequences to gene.seeds and repeat the prepartion steps."
		- Looked through it and there is a lot of diversity but there were at least 3 highly conserved regions that were ~consistent for the whole alignment


### Try running xander with new nirS generesource on same sample


In ```~/tools/RDPTools/Xander_assembler``` make new directory ```mkdir nirS_Cariaco_practice```  

Copy scripts 

```
cp /home/ubuntu/tools/RDPTools/Xander_assembler/bin/xander_setenv.sh /home/ubuntu/tools/RDPTools/Xander_assembler/nirS_Cariaco_practice
cp /home/ubuntu/tools/RDPTools/Xander_assembler/bin/run_xander_skel.sh /home/ubuntu/tools/RDPTools/Xander_assembler/nirS_Cariaco_practice
cp /home/ubuntu/tools/RDPTools/Xander_assembler/bin/get_OTUabundance.sh /home/ubuntu/tools/RDPTools/Xander_assembler/nirS_Cariaco_practice
```

Use same file as used in original try for nirS [D2b267A_merged.fa.gz]. It in in the mounted storage ~/Cariaco/metaG.  

Modify ```nano xander_setenv.sh```

```
SEQFILE=~/Cariaco/metaG/D2b267A_merged.fa
WORKDIR=/home/ubuntu/tools/RDPTools/Xander_assembler/nirS_Cariaco_practice
REF_DIR=/home/ubuntu/tools/RDPTools/Xander_assembler
JAR_DIR=/home/ubuntu/tools/RDPTools
UCHIME=/home/ubuntu/tools/third_party_tools/uchime4.2.40_i86linux32
HMMALIGN=/usr/local/bin/hmmalign
```

Change security of the file ```chmod 755 xander_setenv.sh```  

Make a tmux session and run.  
```
tmux new -s Cariaco_nirSpractice
./run_xander_skel.sh xander_setenv.sh "build find search" "nirS"
```

Saved some output files in ```Cariaco_nirS_test2_D2b267A_merged_Xanderoutput.xlsx```.  
Interesting stuff (when compared to same file but trained using the old, general nirS gene_resource files from run a few days ago: ```Cariaco_nirS_test1_D2b267A_merged_Xanderoutput```:  

- Still clustered into 12 clusters
- Same amount of reads (~220)
- The contribution from Betaproteobacteria, Anaerolinaceae, and hypothetical protein from environmental sample is still about the same
- Originally, ~90% were from gamma proteobacteria. Now half gammas (Oceanospirillales) and half alphas (Rhodobacteraceae)



###Next steps###

- Left some files uploading to EBS overnight and files backing up to external hard drive
- Continue making gene_resources for nirK, etc
- Read what Taylor has to say about picking parameters
- Finish uploading all files


##**July 11th 2018** 
### Build new gene_resource for nirK
*My thought is that there has been a lot of new info on nirK in the last ~ 5 years so it's possible that the previous gene_resource in the tutorial did not have sufficient coverage to pull out nirK from our dataset*

1. Download nirK.seeds and nirK.hmm files from fungene
2. Set minimum HMM coverage to 80%.  
	a. reduced number of sequences from 9899 to 2043.  
	b. Sorting by HMM coverage shows range is 99.2% to 80.9%.  
3. 	Sort by size
	a. Ranges 559 to 312.  
	b. Models in nirK.seeds have size range of 361-379.  
	c. Second-to-last smallest sequence is Nitrobacter, which I would like to leave in the dataset if possible.  
	d. Min size (312) is not that different than seed size (360) so did no filtering based on size
4. Sort by score.  
	a. Max 653.8, Min ~54.  
	b. End of score list contains interesting sequences for us such as Acidivorax and Bacillus species.  
	c. The biggest drop off in score I can find is from 300.5 to 226.6. Cut off anything with score less than 300
	- This will cut off some interesting things like Thiomonas, Bradyrhizobium. The "good" scores are overwhelmed by Nesseria. But there is still some diversity.

	d. Reduced number of sequences from 2043 to 1875
	
5. **Final parameters** = Minimum Score: 300, Minimum Size aa: no limit, Minimum HMM Coverage: 80

6. Manually went through looking for weird ones. Used blastp to check some:

	a. Checked ACF98152: uncultured bacterium 1116 putative nitrite reductase copper-containing NirK seems liks nirK from Bradyrhizobium.  
	b. Checked EGE13573: Moraxella catarrhalis 46P47B1 nitrite reductase AniA/Msp78. Came back with lots of hits as nirK from Moraxella and Psychrobacter.  
	c. Checked CAD32058: Moraxella catarrhalis Sequence 1 from Patent WO0218595. Came back with lots of hits as nirK from Moraxella and Psychrobacter.  
	d. Rest look OK.    

7. Make gene_resource folder and test for the one metagenome sample.  
	a. Selected all remaining sequences and, clicked "begin analysis" at the top right of the page.  
	b. Downloaded protein and nucleotide fastas. Chose accno, fasta, not aligned.  
	c. Changed names to framebot.fa (protein) and nucl.fa (nucleotide) respectively. Navigate to nirK/originaldata in virtual machine. Move all 4 files to originaldata.   
	
	


	- Log onto EC2 instance and mount EBS volume: 

	```
ssh -i amazon-key.pem ubuntu@ec2-35-172-134-113.compute-1.amazonaws.com
sudo mount /dev/xvdf Cariaco
```

	- Check that volume is mounted with ```lsblk``` and ```df```
	- Navigate to gene_resource/nirK, delete old files, and upload 4 new files. Check to make sure the files contain 1875 sequences:  

	```
grep ">" nucl.fa | wc
grep ">" framebot.fa | wc
```
	**The framebot file has 1875 sequences and the nucl file has 1858 and I can't figure out why.**

8. From these 4 files, need to make three new files: ```for_enone.hmm```, ```rev_enone.hmm```, and ```ref_aligned.faa```
	
	- Changed ```prepare_gene_ref.sh``` script:
		- line 4: change gene to your gene of interest (nirK)
		- line 9: change the jar directory to /home/ubuntu/tools/RDPTools [already did in previous for nirS]
	- line 10: change the reference directory to /home/ubuntu/tools/RDPTools/Xander_assembler [already did in previous for nirS]
	- Run ```./prepare_gene_ref.sh nirK```
	- Check ref_aligned.faa file in MEGA or Seaview to see quality of alignment
		- Nothing stands out as super weird. There are definitely some conserved regions

9.  Navigate to Cariaco_nirSpractice folder and run script for nirK in tmux. [shouldn't take so long since it's not rewriting contigs]

```
tmux new -s Cariaco_nirKpractice 
./run_xander_skel.sh xander_setenv.sh "build find search" "nirK"
```

- Didn't really work. It said it assembled 5 contigs but then didn't merge them.
- Try again with free-living community sample. I am still uploading files (it's taking days) but there are some already in there. Try file D3b234B_merged.fa
- Need to modify xander_setenv.sh to ```SEQFILE=~/Cariaco/metaG/D3b234B_merged.fa ``` first. It's OK at this point to overwrite previous results.
- Clear "nirK" folder in k45
- Run:

```
tmux new -s Cariaco_nirKpractice 
./run_xander_skel.sh xander_setenv.sh "build find search" "nirK"
```

	
**I REALIZED I DIDN'T WRITE NEW ```for_enone.hmm```,  ```rev_enone.hmm``` AND ```ref_aligned.faa``` FILES FOR NIRK OR NIRS GENE_RESOURCE. CHECK**

- OK rewrote those files for nirS, checked alignment (added notes above in section where I originally built nirS gene_resource), deleted k45 bloom filter, and running xander for nirS again on D2b267A_merged.fa  (in tmux)
- Result is very similar (~1/2 gamma and ~1/2 alpha nirS). There are 2 more organisms that are new that are in very low abundance but they seem like real candidates for Cariaco. 

```
Taxon	Abundance	Fraction Abundance
coded_by=70358..71974,organism=Bellilinea caldifistulae,definition=nitrite reductase	5.335	0.0230
coded_by=2999064..3000812,organism=Sulfuritalea hydrogenivorans sk43H,definition=hydroxylamine reductase	2.525	0.0109
coded_by=20116..21744,organism=Oleiphilus sp. HI0043,definition=nitrite reductase	100.442	0.4326
coded_by=complement(95538..97643),organism=Arenimonas donghaensis DSM 18148 = HO3-R19,definition=hypothetical protein	1.628	0.0070
coded_by=complement(2455749..2457380),organism=Oleispira antarctica RB-8,definition=Nitrite reductase	2.995	0.0129
coded_by=complement(50734..52695),organism=Labrenzia alba,definition=Nitrite reductase precursor	117.077	0.5042
coded_by=2195004..2196734,organism=Dechlorosoma suillum PS,definition=cytochrome c, mono- and diheme variants family	2.204	0.0095


Lineage	MatchName	Abundance	Fraction Abundance
coded_by=complement(50734..52695),organism=Labrenzia alba,definition=Nitrite reductase precursor	CTQ53626	117.077	0.5042
coded_by=70358..71974,organism=Bellilinea caldifistulae,definition=nitrite reductase	KPL78197	5.335	0.0230
coded_by=complement(95538..97643),organism=Arenimonas donghaensis DSM 18148 = HO3-R19,definition=hypothetical protein	KFL36640	1.628	0.0070
coded_by=2195004..2196734,organism=Dechlorosoma suillum PS,definition=cytochrome c, mono- and diheme variants family	AEV26447	2.204	0.0095
coded_by=complement(2455749..2457380),organism=Oleispira antarctica RB-8,definition=Nitrite reductase	CCK76416	2.995	0.0129
coded_by=2999064..3000812,organism=Sulfuritalea hydrogenivorans sk43H,definition=hydroxylamine reductase	BAO30669	2.525	0.0109
coded_by=20116..21744,organism=Oleiphilus sp. HI0043,definition=nitrite reductase	KZY30905	100.442	0.4326
```



###Next Steps###
  
- Modify nirK gene resource and check if can get any contigs
- Keep uploading files to passport and EBS


##**July 12th 2018** 
### Finish building new gene_resource for nirK
- Fixed gene_resource so that ./prepare_gene_ref.sh nirK references nirK directy. Re-ran it, checked nirK alignments, all looks good. 	
	- Re-ran xander for nirK on ```D2b267A_merged.fa```. No results.
	- Try again on ```D2b267B_merged.fa``` (first change '''xander_setenv.sh''' file in SEQFILE line and delete k45 folder). 
		- Still didn't work.
	- Try again on different depth: D2b247B_merged.fa. (This time ran for nirS and nirK)
		- Worked for nirS and still not nirK!
	- Try one more!: D2a237B_merged.fa
		- Finally got a nirK!!! Only one contig built. It's related to nirK from Rhodopseudomonas palustris WOOHOO!  

		
		
- Also fixed nirS gene resource by changing ```prepare_gene_ref.sh``` line 4. Then re-ran as well.
		- Tried on D2a237B_merged.fa
			- Worked well, similar results to previous files that I tried. 

			


	
### Build new gene_resource for recA housekeeping gene

- On fungene, set min HMM_coverage to 80%
	- Reduced from >103,000 to 70,651 sequences
- Score ranges from 686.8 to 63.9
	- There is rapid decline in scores below ~450
	- Cut scores <450
	- 70,154 sequences remain
- Increase minimum HMM_coverage to 90
	- 68475 sequences remain
- Size ranges between 790 and 307
	- According to E coli wiki, the size is 353 aa
	- Cut anything less than 340 aa
	- 67,333 sequences remain
- This is a huge amount of seqs but I don't know if I can justify cutting any more
- Set Hmm_coverage to 95. Remaining sequences = 38,302. Use these
	- Can't! Only up to 10,000 sequences can be downloaded at a time.
	- Also as I increase the HMM coverage it reduces the diversity dramatically
- Went back to 80% HMM coverage
- Just pick and choose sequences so there are not too many repeats. Eg. we don't need the Mycobacterium sequence 500 times
- **Parameters** = Minimum score: 350, Minimum size: 340, Minimum HMM_coverage: 80%
	- Then "Select all" and download (10,000 at a time)
		- (Make sure to sort by organism name and let whole page load, all 1000 sequences, before selecting) 
	- Concatenate all 7 files into one text file **Continued until July 16th:**
	- Delete repeats as much as possible (in both framebot and nucl files) to reduce size of these files. These are highly weighed by model organisms and do not reflect true diversity expected in Cariaco:
		- Deleted lines 573-3362, which were all Acinetobacter baumanii sequences
			- Everything between ```>JMUH01000009  location=17257..18306,organism=Acinetobacter baumanii BIDMC 56,definition=protein recA
``` and ```>CP006768  location=complement(1552931..1553980),organism=Acinetobacter baumannii ZW85-1,definition=recA protein```

		- Also lots of Mycobacterium abscessus sequences. Delete everything betweeen (new) lines 50775- 53600
			- Between ```>CPW49244  coded_by=74853..76553,organism=Mycobacterium abscessus,definition=recombinase A``` and ```>EPZ21165  coded_by=complement(85986..87026),organism=Mycobacterium abscessus V06705,definition=recombinase RecA
```
	- Also lots of Mycobacterium tuberculosis sequences. Delete between (new) lines 52033-59062
		- Between ```>CBF64426  coded_by=1..1053,organism=Mycobacterium tuberculosis,definition=Sequence 12933 from Patent EP2096177``` and ```>KBS86752  coded_by=697326..699698,organism=Mycobacterium tuberculosis XTB13-290,definition=protein RecA```
	- There is lots from Bordetella pertussis. Delete between 8163-9464
		- Between ```>CPM29417  coded_by=complement(9501..10562),organism=Bordetella pertussis,definition=recombinase A``` and ```>CAE42821  coded_by=complement(266424..267485),organism=Bordetella pertussis Tohama I,definition=RecA protein```
	- Lots of Burkholderia pseudomallei. Delete L10877-12132.
		- Between ```>LWXF01000105  location=complement(15728..16798),organism=Burkholderia pseudomallei,definition=DNA recombination/repair protein RecA``` and ```>JQGY01000142  location=6327..7397,organism=Burkholderia pseudomallei TSV5,definition=protein RecA```
		- Clostridioides difficile
			- Delete L15355-15938
			- Between```>SJS77512  coded_by=complement(201313..202359),organism=Clostridioides difficile,definition=Recombinase A``` and ```>AKP42273  coded_by=1437668..1438714,organism=Clostridioides difficile ATCC 9689 = DSM 1296,definition=recombinase A```
		- Enterobacter cloacae
			- Delete L19025-20000
			- Between ```>JPPR01000030  location=complement(5409..6467),organism=Enterobacter cloacae,definition=protein RecA``` and ```>AYIG01000001  location=complement(115772..116830),organism=Enterobacter cloacae UCICRE 9,definition=protein recA``` 
		- Enterococcus faecalis
			- Delete L19893-20702
			- Between ```>CBU86263  coded_by=1..1047,organism=Enterococcus faecalis,definition=Sequence 6644 from Patent EP2199304``` and ````>EPI28476  coded_by=complement(72085..73131),organism=Enterococcus faecalis WKS-26-18-2,definition=RecA protein```
		- Enterococcus faecium 
			- Delete L19897-21080
			- Between ```>KST46523  coded_by=complement(55445..56494),organism=Enterococcus faecium,definition=DNA recombination/repair protein RecA``` and ```>EZP94860  coded_by=complement(55381..56430),organism=Enterococcus faecium VSE1036,definition=protein RecA```
		- E coli!
			- L20233-28282
			- Between ```>KLW98841  coded_by=complement(688710..689771),organism=Escherichia coli,definition=protein RecA``` and ```>AFJ30374  coded_by=complement(3474300..3475361),organism=Escherichia coli Xuzhou21,definition=recombinase A```
		- Helicobacter pylori
			- L22901-23878
			- Between ```>MILD01000041  location=31386..32429,organism=Helicobacter pylori,definition=recombinase RecA``` and ```>CP003419  location=187939..188982,organism=Helicobacter pylori XZ274,definition=recombinase A```
		- Klebsiella pneumoniae
			- L 23560-27166
			- Between ```>AP014950  location=1124645..1125703,organism=Klebsiella pneumoniae,definition=RecA protein```` and ```>APWD01000062  location=complement(4636..5694),organism=Klebsiella pneumoniae VAKPC309,definition=RecA protein```
		- Listeria monocytogenes 
			- L 28235-29704
			- Between ```>FFHQ01000002  location=complement(648857..649903),organism=Listeria monocytogenes,definition=Recombinase A``` and ```>CP007210  location=814873..815919,organism=Listeria monocytogenes WSLC1042,definition=recombinase RecA```
		- Neisseria gonorrhoeae
			- L 32261-32990
			- Between ```>CFLE01000037  location=complement(7834..8880),organism=Neisseria gonorrhoeae,definition=recombinase A``` and ```>ATPJ01000247  location=complement(185..1312),organism=Neisseria gonorrhoeae SK8976,definition=recombinase RecA```
		- Neisseria meningitidis
			- L32275-33832
			- Between ```>FFAW01000050  location=complement(190..1269),organism=Neisseria meningitidis,definition=recombinase A``` and ```>AL157959  location=complement(1574167..1575213),organism=Neisseria meningitidis Z2491,definition=RecA protein```
		- Peptoclostridium difficile
			- L 35115-35578
			- Between ```>FJUH01000024  location=65843..66913,organism=Peptoclostridium difficile,definition=recombinase A``` and ```>AVKT01000021  location=103559..104605,organism=Peptoclostridium difficile Y41,definition=protein RecA```
		- Propionibacterium acnes
			- L36341-37028
			- Between ```>AEC46661  coded_by=1..>1023,organism=Propionibacterium acnes,definition=RecA``` and ```>AEW79077  coded_by=1061691..1062737,organism=Propionibacterium acnes TypeIA2 P.acn33,definition=protein RecA```
		- Pseudomonas aeruginosa 
			- L36925-38766
			- Between ```>KSP28737  coded_by=complement(43718..44755),organism=Pseudomonas aeruginosa,definition=DNA recombination/repair protein RecA``` and ```>ETU80718  coded_by=complement(177054..178094),organism=Pseudomonas aeruginosa Z61,definition=protein recA```
		- Salmonella enterica
			- L40969-L51840
			- Between ```>KTO00828  coded_by=72134..73195,organism=Salmonella enterica,definition=recombinase A``` and ```>ESE60056  coded_by=96469..97530,organism=Salmonella enterica subsp. salamae serovar 58:l,z13,z28:z6 str. 00-0163,definition=recombinase A```
		- Serratia marcescens
			- L41083-41676
			- Between ```>KFL05907  coded_by=17177..18286,organism=Serratia marcescens,definition=protein RecA``` and ```>AGE16630  coded_by=925691..926755,organism=Serratia marcescens WW4,definition=DNA strand exchange and recombination protein with protease and nuclease activity```
		- Shigella sonnei
			- L41479-43406
			- Between ```>AKH31711  coded_by=complement(4529005..4530066),organism=Shigella sonnei,definition=recombinase A``` and ```>EJL14482  coded_by=complement(41218..42279),organism=Shigella sonnei str. Moseley,definition=protein RecA```
		-  Staphylococcus aureus	
			- L42189-55230
			- Between ```>CVOU01000018  location=420153..421196,organism=Staphylococcus aureus,definition=multifunctional SOS repair factor``` and ```>JIFC01000001  location=406354..407397,organism=Staphylococcus aureus ZTA11/03130-3ST,definition=protein recA```
		- Streptococcus agalactiae 
			- L43195-44416
			- Between ```>CP021864  location=complement(1984793..1985932),organism=Streptococcus agalactiae,definition=DNA recombination/repair protein RecA``` and ```>AKAP01000002  location=17903..19042,organism=Streptococcus agalactiae ZQ0910,definition=recombinase A```
		- Streptococcus pneumoniae
			- L444403-59044
			- Between ```>HC931216  location=26..1192,organism=Streptococcus pneumoniae,definition=Sequence 6482 from Patent EP2199304``` and ```>AE005672  location=complement(1843413..1844579),organism=Streptococcus pneumoniae TIGR4,definition=recA protein```
		- Streptococcus pyogenes 
			- L44439-45078
			- Between ```>CP014027  location=839555..840691,organism=Streptococcus pyogenes,definition=DNA recombination/repair protein RecA``` and ```>AVCG01000046  location=complement(88731..89867),organism=Streptococcus pyogenes UTSW-2,definition=RecA protein``` 
		- Streptococcus suis
			- L44629-45498
			- Between ```>FIHS01000009  location=50789..51940,organism=Streptococcus suis,definition=recombinase A``` and ```>CP006645  location=70168..71319,organism=Streptococcus suis YB51,definition=RecA protein```
		- Vibrio cholerae
			- L47971-49050
			- Between ```>HB902267  location=1..1065,organism=Vibrio cholerae,definition=Sequence 34942 from Patent EP2096177``` and ```>AAKJ02000049  location=complement(4710..5948),organism=Vibrio cholerae V52,definition=recA protein```
		- Vibrio parahaemolyticus
			- L48197-48986
			- Between ```>OAR37860  coded_by=81586..82629,organism=Vibrio parahaemolyticus,definition=DNA recombination/repair protein RecA``` and ```>EXJ49729  coded_by=complement(3641..4684),organism=Vibrio parahaemolyticus VPTS-2010_2,definition=protein RecA```
		- Yersinia enterocolitica
			- L50137-50418
			- Between ```>CFB70057  coded_by=89980..91044,organism=Yersinia enterocolitica,definition=recombinase A``` and ```>CBX73832  coded_by=2314..3378,organism=Yersinia enterocolitica W22703,definition=protein recA```
	- Yersinia pestis
		- L50277-50602
		- Between ```>CAA53084  coded_by=399..1469,organism=Yersinia pestis,definition=RecA protein``` and ```>ADE65819  coded_by=complement(3737765..3738835),organism=Yersinia pestis Z176003,definition=recombinase A```


- **Reduced to 25,232 sequences!**

- Now finish setting up gene_resource

```
mkdir recA
cd ~/tools/RDPTools/Xander_assembler/gene_resource/
cd recA
mkdir originaldata
```

- upload 4 files to his directory. Then

```
cd /home/ubuntu/tools/RDPTools/Xander_assembler/bin
nano prepare_gene_ref.sh
```

- line 4: change gene to your gene of interest
- line 9: change the jar directory to /home/ubuntu/tools/RDPTools
- line 10: change the reference directory to /home/ubuntu/tools/RDPTools/Xander_assembler 
- Run ./prepare_gene_ref.sh recA
- Manually examine the alignment in ref_aligned.faa using Seaview to spot any badly aligned sequences.
	- "If found, it is likely there are no sequences in gene.seeds. Close these sequences. You need to validate these problem sequences to see if they are from the gene you are interested in, then either remove them or add some representative sequences to gene.seeds and repeat the prepartion steps."
	- Looked through it. Conserved regions look good although there are some very short sequences. Go with it for now. (Jump to <u>Jul 16th</u> section and run recA analysis for a sample)

	
##**July 13th 2018** 
- All files uploaded to EBS volume and unzipped!!! Finally	
- Since I am satisfied with the nirS and nirK gene_resources, while I am compiling the resources for recA, try running xander for nirS, nirK on multiple metagenome samples at once

- Make new directories in Cariaco for results. Copy scripts there

	```
cd ~/Cariaco
mkdir Results	
cp /home/ubuntu/tools/RDPTools/Xander_assembler/bin/xander_setenv.sh /home/ubuntu/Cariaco/Results
cp /home/ubuntu/tools/RDPTools/Xander_assembler/bin/run_xander_skel.sh /home/ubuntu/Cariaco/Results
cp /home/ubuntu/tools/RDPTools/Xander_assembler/bin/get_OTUabundance.sh /home/ubuntu/Cariaco/Results
```
- Modify script ```nano xander_setenv.sh``` so that the seqfile finds multiple Metagenome files. For now try only 2 (using wildcard)

	```
SEQFILE=~/Cariaco/metaG/D2b267*_merged.fa
WORKDIR=/home/ubuntu/Cariaco/Results
REF_DIR=/home/ubuntu/tools/RDPTools/Xander_assembler
JAR_DIR=/home/ubuntu/tools/RDPTools
UCHIME=/home/ubuntu/tools/third_party_tools/uchime4.2.40_i86linux32
HMMALIGN=/usr/local/bin/hmmalign
```
- Also changed sample short name to ```SAMPLE_SHORTNAME=D2b267```
- Change permissions and run in tmux

	```
chmod 755 xander_setenv.sh
tmux new -s Cariaco_nir_0713 
./run_xander_skel.sh xander_setenv.sh "build find search" "nirK nirS"
```

- Results = one k45 folder with one cluster folder and 9 different contigs (Which I've seen before.. I am pretty sure these are the contigs individually from D2b267A and D2b267B)
	- I BELIEVE that the way that I ran this was that it ran all the sequences from both samples as one sample.
	- pause for now to go do field work with Mary



	
##**July 16th 2018** 
- Check out results from last week
	- Folder k45 in Cariaco/Results is based on the run for D2b267A and D2b267B
	- Run get_OTUabundance to see if it can extract counts from each sample
		- First modify the copied script:
		- ```nano get_OTUabundance.sh``` line 17 seems to be correct
		- Make directory for results ```mkdir OTUabundanceresults```
		- Run 
		- ```./get_OTUabundance.sh k45/nirS/cluster/D2b267_nirS_45_coverage.txt  OTUabundanceresults 0.1 0.1 k45/nirS/cluster/D2b267_nirS_45_final_prot_aligned.fasta```
		- That didn't work. It treated both samples as one and split all contigs from the 2 files into OTUs. Did not split by sample.
		- According to Taylor's ebooklet, this must be done on multiple final_prot_aligned.fasta files from separate samples. So Xander must be run individually on each sample first
- Try this strategy- Run xander on multiple files (not as one cluster but sequentially) in a directory
	- Practice in directory "combo_practice" which contains D2b267A and R2b267A (in separate metaG and metaT directories. Later make one directory for each sample). Add scripts to combo_practice. 
	
	```
	mkdir metaG
	cp ~/Cariaco/metaG/D2b267A_merged.fa metaG
	
	mkdir metaT
	cp ~/Cariaco/metaT/Nov/R2b267A_merged.fa metaT
	
	cp /home/ubuntu/tools/RDPTools/Xander_assembler/bin/xander_setenv.sh /home/ubuntu/Cariaco/combo_practice
	cp /home/ubuntu/tools/RDPTools/Xander_assembler/bin/run_xander_skel.sh /home/ubuntu/Cariaco/combo_practice
	cp /home/ubuntu/tools/RDPTools/Xander_assembler/bin/get_OTUabundance.sh /home/ubuntu/Cariaco/combo_practice
```

- Modify ```xander_setenv.sh``` so the path to the file is the present working directory (then will change the working directory within loop). Use ```$pwd``` which is a variable that calls the present working directory:

```
SEQFILE=$pwd/.fa
WORKDIR=$pwd  
REF_DIR=/home/ubuntu/tools/RDPTools/Xander_assembler
JAR_DIR=/home/ubuntu/tools/RDPTools
UCHIME=/home/ubuntu/tools/third_party_tools/uchime4.2.40_i86linux32
HMMALIGN=/usr/local/bin/hmmalign

## THIS SECTION NEED TO BE MODIFIED, SAMPLE_SHORTNAME WILL BE THE PREFIX OF CONTIG ID
SAMPLE_SHORTNAME=$i  
```

- and change permissions ```chmod 755 xander_setenv.sh```


- Now copy modified script into each gene directory

```
echo meta* | xargs -n 1 cp xander_setenv.sh
echo meta* | xargs -n 1 cp run_xander_skel.sh
```

- Run individually by cd'ing into each sample directory and calling the script:

```
tmux new -s Cariaco_nirS

for i in $(\ls -d meta*)
do
	cd ~/Cariaco/combo_practice/$i/
	./run_xander_skel.sh xander_setenv.sh "build find search" "nirS";
	cd ~/Cariaco/combo_practice
done

```

*This seems to have worked. Make sure it is putting output in gene directory as well.*

Output

```
### Build bloom filter
java -Xmx2G -jar /home/ubuntu/tools/RDPTools/hmmgs.jar build /home/ubuntu/Cariaco/combo_practice/metaG/*.fa k45.bloom 45 32 2 4 30 >& k45_bloom_stat.txt
### Find starting kmers for nirS
java -Xmx2G -jar /home/ubuntu/tools/RDPTools/KmerFilter.jar fast_kmer_filter -a -o temp_starts_1.txt -t 1 45 /home/ubuntu/Cariaco/combo_practice/metaG/D2b267A_merged.fa nirS=/home/ubuntu/tools/RDPTools/Xander_assembler/gene_resource/nirS/ref_aligned.faa 
Starting kmer mapping at Mon Jul 16 17:09:56 UTC 2018
*  Number of threads:       1
*  References:              [nirS]
*  Reads file:              /home/ubuntu/Cariaco/combo_practice/metaG/D2b267A_merged.fa
*  Kmer length:             15
*  Kmer Refset Size:        62865
Processed 1000001 sequences in 58679 ms
Processed 2000000 sequences in 115648 ms
Processed 2999999 sequences in 171351 ms
Processed 4000000 sequences in 230103 ms
Processed 5000001 sequences in 286341 ms
Finished Processed 5827663 sequences in 335410 ms
### Search contigs nirS
java -Xmx2G -jar /home/ubuntu/tools/RDPTools/hmmgs.jar search -p 20 1 100 ../k45.bloom /home/ubuntu/tools/RDPTools/Xander_assembler/gene_resource/nirS/for_enone.hmm /home/ubuntu/tools/RDPTools/Xander_assembler/gene_resource/nirS/rev_enone.hmm gene_starts.txt 1> stdout.txt 2> stdlog.txt
### Merge contigs
java -Xmx2G -jar /home/ubuntu/tools/RDPTools/hmmgs.jar merge -a -o merge_stdout.txt -s test -b 50 --min-length 150 /home/ubuntu/tools/RDPTools/Xander_assembler/gene_resource/nirS/for_enone.hmm stdout.txt gene_starts.txt_nucl.fasta
Read in 527 contigs, wrote out 267 merged contigs in 1.602s
Warning, no derep mode specified, using unaligned
Processing prot_merged.fasta
Total sequences: 267
Unique sequences: 17
Dereplication complete: 149
### Cluster
Merging file alignment/aligned.stk
Using #=GC_RF as mask sequence
Processing aligned.fasta
Total sequences: 16
Unique sequences: 16
Dereplication complete: 91
Reading sequences(memratio=0.0013843589809955199)...
Using distance model edu.msu.cme.rdp.alignment.pairwise.rna.IdentityDistanceModel
Read 16 Protein sequences (memratio=0.002076515421962974)
Reading ID Mapping from file /home/ubuntu/Cariaco/combo_practice/metaG/k45/nirS/cluster/ids
Read mapping for 16 sequences (memratio=0.002076515421962974)
Starting distance computations, predicted max edges=256, at=Mon Jul 16 17:16:03 UTC 2018
Dumping 44 edges to partial_matrix0 FINAL EDGES (memory ratio=0.0016976859143264212)
Matrix edges computed: 23
Maximum distance: 0.625
Splits: 1
Partition files merged: 7
Doing complete linkage clustering with step 0.009999999776482582 (realstep=100)
Clustering complete: 30
Converted 15 sequences from [complete.clust_rep_seqs.fasta] (FASTA) to fasta in 0 s
### Chimera removal
uchime v4.2.40
by Robert C. Edgar
http://drive5.com/uchime
This code is donated to the public domain.

00:00  18Mb  100.0% Reading test_nirS_45_nucl_rep_seqs.fasta
00:00  18Mb 15 sequences                                    
00:00  20Mb  100.0% Reading /home/ubuntu/tools/RDPTools/Xander_assembler/gene_resource/nirS/originaldata/nucl.fa
00:00  20Mb 1251 sequences                                                                                      
00:02 9.9Mb  100.0% 0/14 chimeras found (0.0%)
### FrameBot
java -jar /home/ubuntu/tools/RDPTools/FrameBot.jar framebot -N -l 150 -o nirS_45 /home/ubuntu/tools/RDPTools/Xander_assembler/gene_resource/nirS/originaldata/framebot.fa nucl_rep_seqs_rmchimera.fasta
### Kmer abundance
java -Xmx2g -jar /home/ubuntu/tools/RDPTools/KmerFilter.jar kmer_coverage -t 1 -m test_nirS_45_match_reads.fa 45 test_nirS_45_final_nucl.fasta test_nirS_45_coverage.txt test_nirS_45_abundance.txt /home/ubuntu/Cariaco/combo_practice/metaG/*.fa
```

- Notice here it automatically goes to second sample directory- woohoo!

```### Build bloom filter
java -Xmx2G -jar /home/ubuntu/tools/RDPTools/hmmgs.jar build /home/ubuntu/Cariaco/combo_practice/metaT/*.fa k45.bloom 45 32 2 4 30 >& k45_bloom_stat.txt
### Find starting kmers for nirS
java -Xmx2G -jar /home/ubuntu/tools/RDPTools/KmerFilter.jar fast_kmer_filter -a -o temp_starts_1.txt -t 1 45 /home/ubuntu/Cariaco/combo_practice/metaT/R2b267A_merged.fa nirS=/home/ubuntu/tools/RDPTools/Xander_assembler/gene_resource/nirS/ref_aligned.faa 
Starting kmer mapping at Mon Jul 16 17:48:29 UTC 2018
*  Number of threads:       1
*  References:              [nirS]
*  Reads file:              /home/ubuntu/Cariaco/combo_practice/metaT/R2b267A_merged.fa
*  Kmer length:             15
*  Kmer Refset Size:        62865
Processed 1000001 sequences in 39225 ms
Processed 2000000 sequences in 79432 ms
Processed 3000000 sequences in 118317 ms
Processed 3999999 sequences in 155569 ms
Processed 5000000 sequences in 197435 ms
Processed 5999999 sequences in 234605 ms
Processed 6999999 sequences in 271750 ms
Processed 7999999 sequences in 310507 ms
Processed 8999999 sequences in 352132 ms
Processed 9999999 sequences in 387384 ms
Processed 10999999 sequences in 426250 ms
Processed 11999999 sequences in 465027 ms
Processed 12999999 sequences in 503631 ms
Processed 13999999 sequences in 543355 ms
Processed 14999999 sequences in 584201 ms
Processed 15999999 sequences in 623078 ms
Processed 16999999 sequences in 663863 ms
Processed 17999999 sequences in 708867 ms
Processed 18999999 sequences in 747729 ms
Processed 19999999 sequences in 784423 ms
Processed 20999999 sequences in 825096 ms
Processed 21999999 sequences in 863098 ms
Processed 22999999 sequences in 900957 ms
Processed 23999999 sequences in 937120 ms
Processed 24999999 sequences in 979662 ms
Finished Processed 25395972 sequences in 995338 ms
### Search contigs nirS
java -Xmx2G -jar /home/ubuntu/tools/RDPTools/hmmgs.jar search -p 20 1 100 ../k45.bloom /home/ubuntu/tools/RDPTools/Xander_assembler/gene_resource/nirS/for_enone.hmm /home/ubuntu/tools/RDPTools/Xander_assembler/gene_resource/nirS/rev_enone.hmm gene_starts.txt 1> stdout.txt 2> stdlog.txt
```

- Successfully built k45 kmers but could not build nirS contigs from metaT file
- Doesn't matter for now. The point is I got the script to work in a loop
- Next step: Combine analysis for recA, nirS, and nirK
	- First finish building gene_resource for recA
	- Then try running recA on a metaT file. It must work
- Consider changing filter size in xander script- some files are bigger than the 2G limit for the default filter

### Run analysis for recA###
```
tmux new -s Cariaco_recA

for i in $(\ls -d meta*)
do
    cd ~/Cariaco/combo_practice/$i/
    ./run_xander_skel.sh xander_setenv.sh "build find search" "recA";
    cd ~/Cariaco/combo_practice
done
```

- Tried for both combo_practice metaT and metaG files
- Don't rebuild bloom filter for now
- MetaT is still not building. There are contigs but they don't merge : ```Read in 506 contigs, wrote out 0 merged contigs in 0.606s```. Look into this tomorrow


	

	
##**July 17th 2018** 

### To do
- Start organizing metaG folder and run nirS, nirK, recA analysis
- Run get_OTUabundance for recA to see if I can get it to separate and count recA OTUs by sample
- Read more and contact Taylor about problem with metaT samples (maybe set lower limit for length of merged contigs?)
- Figure out if there is an option for anammox marker gene (there is no hzo in fungene :( )


<u>Organize metaG samples into sample folders</u>

- Follow Taylor's suggestions from ebooklet on [combining samples](https://github.com/dunivint/RDP_Tutorials/blob/master/xander_combining_samples.md)
- I have been able to modify the xander_setenv script so it can go into the present working directoy and use the .fa file. Furthermore,  can write a loop so that the script is used over and over again in different, similarly named, directories. Therefore organize samples so that each sample in a directory of the same name. Do this in the metaG master folder in mounted volume: ~/Cariaco/
	- NOTE- before moving and potentially modifying the sample files, make copies of sampels in "raw_data" to leave untouched. Can zip these to save space (doing this in tmux session called "organization")

```
cd ~/Cariaco/metaG
for i in $(\ls -d *.fa)
do
	echo $i | cut -c1-7 | xargs -n 1 mkdir
done
for i in $(\ls -d *.fa | cut -c1-7)
do
	mv "$i""_merged.fa" $i
done
```

- Move scripts into metaG folder. 


```
cp /home/ubuntu/tools/RDPTools/Xander_assembler/bin/xander_setenv.sh /home/ubuntu/Cariaco/metaG
cp /home/ubuntu/tools/RDPTools/Xander_assembler/bin/run_xander_skel.sh /home/ubuntu/Cariaco/metaG
cp /home/ubuntu/tools/RDPTools/Xander_assembler/bin/get_OTUabundance.sh /home/ubuntu/Cariaco/metaG

nano xander_setenv.sh
```

- Modify script, set permissions, and then copy into each sample folder.
	- NOTE- this is a chance to optimize parameters for all samples. Follow Taylor's recomendations on [choosing parameters](https://github.com/dunivint/RDP_Tutorials/blob/master/xander_choosing_parameters.md)
	- Leave ```K_SIZE``` as 45
	- ```MIN_COUNT``` should be 1
	- Largest Cariaco metaG sample is 8.3Gb [D3a234A]. According to ebooklet, the ```FILTER_SIZE``` should change according to the file size.
		- 32 is appropriate for 2GB. 35 is appropriate for 6GB. 38 is for 70GB.
		- Choose 36 for now since largest sample is ~8Gb. Most are ~1-3Gb. ~3 samples are 6Gb or more
	- ```MAX_JVM_HEAP``` - max amount of memory for build process. 
		- Must be larger than the size of the bloom filter, which is determined by: 2^(```FILTER_SIZE```-3))/10^9 GB
		- So far a ```FILTER_SIZE``` of 36, this calculates to 8.6GB.
		- Therefore set the ```MAX_JVM_HEAP``` to 16GB, which is more than 8 but not as large at 16 (which corresponds to FILTER_SIZE of 37)
	- ```THREADS``` set to 3. According to ebooklet :"```THREADS``` is the number of computer cores to use. Only one core is used to build the bloom filter, ```THREADS``` does not impact this step.
The find and search steps may be run in parallel, one core for each gene... Set THREADS to the number of genes you are searching for, but do not exceed one less than the number of cores you have on your computer.:
		- This virtual machine (m4.xlarge) has 4 cores (or vCPUS) so setting this to 3 is OK. 


```
SEQFILE=$PWD/*.fa
WORKDIR=$PWD
REF_DIR=/home/ubuntu/tools/RDPTools/Xander_assembler 
JAR_DIR=/home/ubuntu/tools/RDPTools 
UCHIME=/home/ubuntu/tools/third_party_tools/uchime4.2.40_i86linux32 
HMMALIGN=/usr/local/bin/hmmalign

SAMPLE_SHORTNAME=metagenome
```

- Leave rest of parameters as default
- Modify permsision ```chmod 755 xander_setenv.sh```
- Copy this version of the xander script and ```run_xander_skel``` into each sample folder:

```
cd ~/Cariaco/metaG
chmod 755 xander_setenv.sh 

echo D* | xargs -n 1 cp xander_setenv.sh 
echo D* | xargs -n 1 cp run_xander_skel.sh


```


- **And run metaG samples for all 3 genes!** Do with loop similar to practice yesterday

```
tmux new -s Cariaco_metaG
cd ~/Cariaco/metaG/
for i in $(\ls -d D*) 
do
	cd ~/Cariaco/metaG/$i/
	./run_xander_skel.sh xander_setenv.sh "build find search" "recA nirK nirS"; 
	cd ~/Cariaco/metaG
done
```

- Finally got this running after awhile. Took forever to figure out- was having permission problems. All because pwd in script wasn't capitalized!
- Leave running over night. Monitor in tmux once in awhile
- **FINISHED AT ~1PM NEXT DAY :)**


##**July 18th 2018** 

- Still running
- Start organizing metaT files into their own folders and copy scripts into working directories

###Run metaT Files for all 3 genes###
- First prepare folders for each sample

```
cd ~/Cariaco/metaT/May
for i in $(\ls -d *.fa)
do
	echo $i | cut -c1-7 | xargs -n 1 mkdir
done
for i in $(\ls -d *.fa | cut -c1-7)
do
	mv "$i""_merged.fa" $i
done

cd ~/Cariaco/metaT/Nov
for i in $(\ls -d *.fa)
do
	echo $i | cut -c1-7 | xargs -n 1 mkdir
done
for i in $(\ls -d *.fa | cut -c1-7)
do
	mv "$i""_merged.fa" $i
done


cp /home/ubuntu/tools/RDPTools/Xander_assembler/bin/xander_setenv.sh /home/ubuntu/Cariaco/metaT
cp /home/ubuntu/tools/RDPTools/Xander_assembler/bin/run_xander_skel.sh /home/ubuntu/Cariaco/metaT
cp /home/ubuntu/tools/RDPTools/Xander_assembler/bin/get_OTUabundance.sh /home/ubuntu/Cariaco/metaT

nano xander_setenv.sh
```

- Change parameters in xander_Setenv.sh
	- Largest metaT file is 13GB. Some are 5-6GB. Most are 3-6GB. A few files from Nov are very small (>.01GB)
	- ```K_SIZE``` 45
	- ```MIN_COUNT``` 1
	- ```FILTER_SIZE``` 37 (this is appropriate for a 16GB file.)
	- ```MAX_JVM_HEAP``` 32GB (corresponds to FILTER_SIZE of 38) 
	- ```THREADS``` 3. 
	- Leave rest of parameters as default

```
SEQFILE=$PWD/*.fa
WORKDIR=$PWD
REF_DIR=/home/ubuntu/tools/RDPTools/Xander_assembler 
JAR_DIR=/home/ubuntu/tools/RDPTools 
UCHIME=/home/ubuntu/tools/third_party_tools/uchime4.2.40_i86linux32 
HMMALIGN=/usr/local/bin/hmmalign

SAMPLE_SHORTNAME=$name
```


- Modify permsision ```chmod 755 xander_setenv.sh```
- Copy this version of the xander script and ```run_xander_skel``` into each sample folder:

```
echo ~/Cariaco/metaT/*/R* | xargs -n 1 cp xander_setenv.sh 
echo ~/Cariaco/metaT/*/R* | xargs -n 1 cp run_xander_skel.sh
```


- Run (not sure this will be succesful- wasn't able to get too many nirS transcripts from these files) **DIDNT WORK- BUILD BLOOM FILTER FAILED FOR ALL DIRECTORIES IN MAY. SEE WHAT HAPPENED**

```
tmux new -s Cariaco_metaT

cd ~/Cariaco/metaT/May
for name in $(\ls -d R2a103*) 
do
    cd ~/Cariaco/metaT/May/$name/
    ./run_xander_skel.sh xander_setenv.sh "build find search" "recA nirK nirS"; 
    cd ~/Cariaco/metaT/May
done

cd ~/Cariaco/metaT/Nov
for i in $(\ls -d R*) 
do
    cd ~/Cariaco/metaT/Nov/$i/
    ./run_xander_skel.sh xander_setenv.sh "build find search" "recA nirK nirS"; 
    cd ~/Cariaco/metaT/Nov
done
```
- **NOT WORKING**


### Gather MetaG Results
- I didn't give each sample a unique sample_name (because ran with same script)
- Copy results into new directory and change name of results text files to correspond to directory (sample) name



- First copy results files (all of contents of k45 folder) into directory folder of same name without the copying the sample fasta file (because it's huge)

```
cd ~/Cariaco/metaG
mkdir Results
for i in $(\ls -d D*) 
do
    mkdir ~/Cariaco/metaG/Results/$i
    cp -r ~/Cariaco/metaG/$i/k45 ~/Cariaco/metaG/Results/$i
done

```

- Then I dragged all these results into a local folder
- Deleted this folder off instance when I was confident everything was transferred. Need to save space in order to run metaT


- Working locally, within each directory, replace the phrase "metagenome" in file names with the sample name (also name of directory). This has to be done to files within the gene folder (nirS, nirK, etc) and in the subfolder within that called "cluster" **CANT GET THIS TO WORK. THIS IS WHAT I'VE TRIED. LOOK AT CODE ON THIS [PAGE](https://unix.stackexchange.com/questions/181141/rename-multiple-files-with-mv-to-change-the-extension)**
	- Also moved these to my external HD and working locally. Instance is off.
	- Needed to first install Homebrew and the function "rename"

```
/usr/bin/ruby -e "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/master/install)"

brew install rename

cd /Volumes/MyPassport/CariacoNpaper/xander_results/metaG/
for i in $(\ls -d D*test) 
do	
cd /Volumes/MyPassport/CariacoNpaper/xander_results/metaG/$i/k45/nirK
rename 's/metagenome/"$i"/' *
cd cluster
rename 's/metagenome/"$i"/' *
	
cd /Volumes/MyPassport/CariacoNpaper/xander_results/metaG/$i/k45/nirS
rename 's/metagenome/"$i"/' *
cd cluster
rename 's/metagenome/"$i"/' *
	
cd /Volumes/MyPassport/CariacoNpaper/xander_results/metaG/$i/k45/recA
rename 's/metagenome/"$i"/' *
cd cluster
rename 's/metagenome/"$i"/' *
done
```

- CRAP I can't get the loop to work. Some weird PERL scripting errors. Need to do manually- cd into each /k45/gene and /k45/gene/cluster directory and run. This is still faster then doing it in by copy paste in Finder

```
cd /Volumes/MyPassport/CariacoNpaper/xander_results/metaG/D2a143A/k45/recA
rename 's/metagenome/D2a143A/' *
cd cluster
rename 's/metagenome/D2a143A/' *
cd ../../nirS
rename 's/metagenome/D2a143A/' *
cd cluster #IN THIS CASE NO DATA IN THIS FOLDER SO DIDN'T REQUIRE NEXT LINE
rename 's/metagenome/D2a143A/' *
cd ../../nirK # NO RESULTS
rename 's/metagenome/D2a143A/' *
cd cluster
rename 's/metagenome/D2a143A/' *

cd /Volumes/MyPassport/CariacoNpaper/xander_results/metaG/D2a143B/k45/recA
rename 's/metagenome/D2a143B/' *
cd cluster
rename 's/metagenome/D2a143B/' *
cd ../../nirS # NO RESULTS
rename 's/metagenome/D2a143B/' *
cd cluster 
rename 's/metagenome/D2a143B/' *
cd ../../nirK
rename 's/metagenome/D2a143B/' *
cd cluster
rename 's/metagenome/D2a143B/' *

cd /Volumes/MyPassport/CariacoNpaper/xander_results/metaG/D2a200A/k45/recA
rename 's/metagenome/D2a200A/' *
cd cluster
rename 's/metagenome/D2a200A/' *
cd ../../nirS # NO RESULTS
rename 's/metagenome/D2a200A/' *
cd cluster 
rename 's/metagenome/D2a200A/' *
cd ../../nirK # NO RESULTS
rename 's/metagenome/D2a200A/' *
cd cluster
rename 's/metagenome/D2a200A/' *

cd /Volumes/MyPassport/CariacoNpaper/xander_results/metaG/D2a200B/k45/recA
rename 's/metagenome/D2a200B/' *
cd cluster
rename 's/metagenome/D2a200B/' *
cd ../../nirS # NO RESULTS
rename 's/metagenome/D2a200B/' *
cd cluster 
rename 's/metagenome/D2a200B/' *
cd ../../nirK # NO RESULTS
rename 's/metagenome/D2a200A/' *
cd cluster
rename 's/metagenome/D2a200B/' *


cd /Volumes/MyPassport/CariacoNpaper/xander_results/metaG/D2a237A/k45/recA
rename 's/metagenome/D2a237A/' *
cd cluster
rename 's/metagenome/D2a237A/' *
cd ../../nirS 
rename 's/metagenome/D2a237A/' *
cd cluster 
rename 's/metagenome/D2a237A/' *
cd ../../nirK # NO RESULTS
rename 's/metagenome/D2a237A/' *
cd cluster
rename 's/metagenome/D2a237A/' *

cd /Volumes/MyPassport/CariacoNpaper/xander_results/metaG/D2a237B/k45/recA
rename 's/metagenome/D2a237B/' *
cd cluster
rename 's/metagenome/D2a237B/' *
cd ../../nirS 
rename 's/metagenome/D2a237B/' *
cd cluster 
rename 's/metagenome/D2a237B/' *
cd ../../nirK
rename 's/metagenome/D2a237B/' *
cd cluster
rename 's/metagenome/D2a237B/' *

cd /Volumes/MyPassport/CariacoNpaper/xander_results/metaG/D2a247A/k45/recA
rename 's/metagenome/D2a247A/' *
cd cluster
rename 's/metagenome/D2a247A/' *
cd ../../nirS 
rename 's/metagenome/D2a247A/' *
cd cluster 
rename 's/metagenome/D2a247A/' *
cd ../../nirK # NO RESULTS FOR NIRK
rename 's/metagenome/D2a247A/' *
cd cluster
rename 's/metagenome/D2a247A/' *

cd /Volumes/MyPassport/CariacoNpaper/xander_results/metaG/D2a247B/k45/recA
rename 's/metagenome/D2a247B/' *
cd cluster
rename 's/metagenome/D2a247B/' *
cd ../../nirS 
rename 's/metagenome/D2a247B/' *
cd cluster 
rename 's/metagenome/D2a247B/' *
cd ../../nirK # NO RESULTS FOR NIRK
rename 's/metagenome/D2a247B/' *
cd cluster
rename 's/metagenome/D2a247B/' *

cd /Volumes/MyPassport/CariacoNpaper/xander_results/metaG/D2a267A/k45/recA
rename 's/metagenome/D2a267A/' *
cd cluster
rename 's/metagenome/D2a267A/' *
cd ../../nirS 
rename 's/metagenome/D2a267A/' *
cd cluster 
rename 's/metagenome/D2a267A/' *
cd ../../nirK # NO RESULTS FOR NIRK
rename 's/metagenome/D2a267A/' *
cd cluster
rename 's/metagenome/D2a267A/' *

cd /Volumes/MyPassport/CariacoNpaper/xander_results/metaG/D2a267B/k45/recA
rename 's/metagenome/D2a267B/' *
cd cluster
rename 's/metagenome/D2a267B/' *
cd ../../nirS 
rename 's/metagenome/D2a267B/' *
cd cluster 
rename 's/metagenome/D2a267B/' *
cd ../../nirK # NO RESULTS FOR NIRK
rename 's/metagenome/D2a267B/' *
cd cluster
rename 's/metagenome/D2a267B/' *

cd /Volumes/MyPassport/CariacoNpaper/xander_results/metaG/D2b143A/k45/recA
rename 's/metagenome/D2b143A/' *
cd cluster
rename 's/metagenome/D2b143A/' *
cd ../../nirS # NO RESULTS FOR NIRS
rename 's/metagenome/D2b143A/' *
cd cluster 
rename 's/metagenome/D2b143A/' *
cd ../../nirK # NO RESULTS FOR NIRK
rename 's/metagenome/D2b143A/' *
cd cluster
rename 's/metagenome/D2b143A/' *

cd /Volumes/MyPassport/CariacoNpaper/xander_results/metaG/D2b143B/k45/recA
rename 's/metagenome/D2b143B/' *
cd cluster
rename 's/metagenome/D2b143B/' *
cd ../../nirS # NO RESULTS FOR NIRS
rename 's/metagenome/D2b143B/' *
cd cluster 
rename 's/metagenome/D2b143B/' *
cd ../../nirK 
rename 's/metagenome/D2b143B/' *
cd cluster
rename 's/metagenome/D2b143B/' *

cd /Volumes/MyPassport/CariacoNpaper/xander_results/metaG/D2b200A/k45/recA
rename 's/metagenome/D2b200A/' *
cd cluster
rename 's/metagenome/D2b200A/' *
cd ../../nirS 
rename 's/metagenome/D2b200A/' *
cd cluster 
rename 's/metagenome/D2b200A/' *
cd ../../nirK # NO RESULTS FOR NIRK
rename 's/metagenome/D2b200A/' *
cd cluster
rename 's/metagenome/D2b200A/' *

cd /Volumes/MyPassport/CariacoNpaper/xander_results/metaG/D2b200B/k45/recA
rename 's/metagenome/D2b200B/' *
cd cluster
rename 's/metagenome/D2b200B/' *
cd ../../nirS 
rename 's/metagenome/D2b200B/' *
cd cluster 
rename 's/metagenome/D2b200B/' *
cd ../../nirK
rename 's/metagenome/D2b200B/' *
cd cluster
rename 's/metagenome/D2b200B/' *

cd /Volumes/MyPassport/CariacoNpaper/xander_results/metaG/D2b237A/k45/recA
rename 's/metagenome/D2b237A/' *
cd cluster
rename 's/metagenome/D2b237A/' *
cd ../../nirS 
rename 's/metagenome/D2b237A/' *
cd cluster 
rename 's/metagenome/D2b237A/' *
cd ../../nirK # NO RESULTS FOR NIRK
rename 's/metagenome/D2b237A/' *
cd cluster
rename 's/metagenome/D2b237A/' *

cd /Volumes/MyPassport/CariacoNpaper/xander_results/metaG/D2b237B/k45/recA
rename 's/metagenome/D2b237B/' *
cd cluster
rename 's/metagenome/D2b237B/' *
cd ../../nirS 
rename 's/metagenome/D2b237B/' *
cd cluster 
rename 's/metagenome/D2b237B/' *
cd ../../nirK
rename 's/metagenome/D2b237B/' *
cd cluster
rename 's/metagenome/D2b237B/' *

cd /Volumes/MyPassport/CariacoNpaper/xander_results/metaG/D2b247A/k45/recA
rename 's/metagenome/D2b247A/' *
cd cluster
rename 's/metagenome/D2b247A/' *
cd ../../nirS 
rename 's/metagenome/D2b247A/' *
cd cluster 
rename 's/metagenome/D2b247A/' *
cd ../../nirK # NO RESULTS FOR NIRK
rename 's/metagenome/D2b247A/' *
cd cluster
rename 's/metagenome/D2b247A/' *

cd /Volumes/MyPassport/CariacoNpaper/xander_results/metaG/D2b247B/k45/recA
rename 's/metagenome/D2b247B/' *
cd cluster
rename 's/metagenome/D2b247B/' *
cd ../../nirS 
rename 's/metagenome/D2b247B/' *
cd cluster 
rename 's/metagenome/D2b247B/' *
cd ../../nirK # NO RESULTS FOR NIRK
rename 's/metagenome/D2b247B/' *
cd cluster
rename 's/metagenome/D2b247B/' *

cd /Volumes/MyPassport/CariacoNpaper/xander_results/metaG/D2b267A/k45/recA
rename 's/metagenome/D2b267A/' *
cd cluster
rename 's/metagenome/D2b267A/' *
cd ../../nirS 
rename 's/metagenome/D2b267A/' *
cd cluster 
rename 's/metagenome/D2b267A/' *
cd ../../nirK # NO RESULTS FOR NIRK
rename 's/metagenome/D2b267A/' *
cd cluster
rename 's/metagenome/D2b267A/' *

cd /Volumes/MyPassport/CariacoNpaper/xander_results/metaG/D2b267B/k45/recA
rename 's/metagenome/D2b267B/' *
cd cluster
rename 's/metagenome/D2b267B/' *
cd ../../nirS 
rename 's/metagenome/D2b267B/' *
cd cluster 
rename 's/metagenome/D2b267B/' *
cd ../../nirK # NO RESULTS FOR NIRK
rename 's/metagenome/D2b267B/' *
cd cluster
rename 's/metagenome/D2b267B/' *

cd /Volumes/MyPassport/CariacoNpaper/xander_results/metaG/D3a103A/k45/recA
rename 's/metagenome/D3a103A/' *
cd cluster
rename 's/metagenome/D3a103A/' *
cd ../../nirS # NO RESULTS FOR NIRS
rename 's/metagenome/D3a103A/' *
cd cluster 
rename 's/metagenome/D3a103A/' *
cd ../../nirK # NO RESULTS FOR NIRK
rename 's/metagenome/D3a103A/' *
cd cluster
rename 's/metagenome/D3a103A/' *

cd /Volumes/MyPassport/CariacoNpaper/xander_results/metaG/D3a103B/k45/recA
rename 's/metagenome/D3a103B/' *
cd cluster
rename 's/metagenome/D3a103B/' *
cd ../../nirS # NO RESULTS FOR NIRS
rename 's/metagenome/D3a103B/' *
cd cluster 
rename 's/metagenome/D3a103B/' *
cd ../../nirK # NO RESULTS FOR NIRK
rename 's/metagenome/D3a103B/' *
cd cluster
rename 's/metagenome/D3a103B/' *

cd /Volumes/MyPassport/CariacoNpaper/xander_results/metaG/D3a198A/k45/recA
rename 's/metagenome/D3a198A/' *
cd cluster
rename 's/metagenome/D3a198A/' *
cd ../../nirS 
rename 's/metagenome/D3a198A/' *
cd cluster 
rename 's/metagenome/D3a198A/' *
cd ../../nirK 
rename 's/metagenome/D3a198A/' *
cd cluster
rename 's/metagenome/D3a198A/' *

cd /Volumes/MyPassport/CariacoNpaper/xander_results/metaG/D3a198B/k45/recA
rename 's/metagenome/D3a198B/' *
cd cluster
rename 's/metagenome/D3a198B/' *
cd ../../nirS 
rename 's/metagenome/D3a198B/' *
cd cluster 
rename 's/metagenome/D3a198B/' *
cd ../../nirK 
rename 's/metagenome/D3a198B/' *
cd cluster
rename 's/metagenome/D3a198B/' *

cd /Volumes/MyPassport/CariacoNpaper/xander_results/metaG/D3a234A/k45/recA
rename 's/metagenome/D3a234A/' *
cd cluster
rename 's/metagenome/D3a234A/' *
cd ../../nirS 
rename 's/metagenome/D3a234A/' *
cd cluster 
rename 's/metagenome/D3a234A/' *
cd ../../nirK 
rename 's/metagenome/D3a234A/' *
cd cluster
rename 's/metagenome/D3a234A/' *

cd /Volumes/MyPassport/CariacoNpaper/xander_results/metaG/D3a234B/k45/recA
rename 's/metagenome/D3a234B/' *
cd cluster
rename 's/metagenome/D3a234B/' *
cd ../../nirS 
rename 's/metagenome/D3a234B/' *
cd cluster 
rename 's/metagenome/D3a234B/' *
cd ../../nirK 
rename 's/metagenome/D3a234B/' *
cd cluster
rename 's/metagenome/D3a234B/' *

cd /Volumes/MyPassport/CariacoNpaper/xander_results/metaG/D3a295A/k45/recA
rename 's/metagenome/D3a295A/' *
cd cluster
rename 's/metagenome/D3a295A/' *
cd ../../nirS 
rename 's/metagenome/D3a295A/' *
cd cluster 
rename 's/metagenome/D3a295A/' *
cd ../../nirK #NO NIRK
rename 's/metagenome/D3a295A/' *
cd cluster
rename 's/metagenome/D3a295A/' *

cd /Volumes/MyPassport/CariacoNpaper/xander_results/metaG/D3a295B/k45/recA
rename 's/metagenome/D3a295B/' *
cd cluster
rename 's/metagenome/D3a295B/' *
cd ../../nirS 
rename 's/metagenome/D3a295B/' *
cd cluster 
rename 's/metagenome/D3a295B/' *
cd ../../nirK #NO NIRK
rename 's/metagenome/D3a295B/' *
cd cluster
rename 's/metagenome/D3a295B/' *

cd /Volumes/MyPassport/CariacoNpaper/xander_results/metaG/D3a314A/k45/recA
rename 's/metagenome/D3a314A/' *
cd cluster
rename 's/metagenome/D3a314A/' *
cd ../../nirS 
rename 's/metagenome/D3a314A/' *
cd cluster 
rename 's/metagenome/D3a314A/' *
cd ../../nirK #NO NIRK
rename 's/metagenome/D3a314A/' *
cd cluster
rename 's/metagenome/D3a314A/' *

cd /Volumes/MyPassport/CariacoNpaper/xander_results/metaG/D3a314B/k45/recA
rename 's/metagenome/D3a314B/' *
cd cluster
rename 's/metagenome/D3a314B/' *
cd ../../nirS 
rename 's/metagenome/D3a314B/' *
cd cluster 
rename 's/metagenome/D3a314B/' *
cd ../../nirK #NO NIRK
rename 's/metagenome/D3a314B/' *
cd cluster
rename 's/metagenome/D3a314B/' *

cd /Volumes/MyPassport/CariacoNpaper/xander_results/metaG/D3b103A/k45/recA
rename 's/metagenome/D3b103A/' *
cd cluster
rename 's/metagenome/D3b103A/' *
cd ../../nirS #NO NIRS
rename 's/metagenome/D3b103A/' *
cd cluster 
rename 's/metagenome/D3b103A/' *
cd ../../nirK #NO NIRK
rename 's/metagenome/D3b103A/' *
cd cluster
rename 's/metagenome/D3b103A/' *

cd /Volumes/MyPassport/CariacoNpaper/xander_results/metaG/D3b103B/k45/recA
rename 's/metagenome/D3b103B/' *
cd cluster
rename 's/metagenome/D3b103B/' *
cd ../../nirS #NO NIRS
rename 's/metagenome/D3b103B/' *
cd cluster 
rename 's/metagenome/D3b103B/' *
cd ../../nirK #NO NIRK
rename 's/metagenome/D3b103B/' *
cd cluster
rename 's/metagenome/D3b103B/' *

cd /Volumes/MyPassport/CariacoNpaper/xander_results/metaG/D3b198A/k45/recA
rename 's/metagenome/D3b198A/' *
cd cluster
rename 's/metagenome/D3b198A/' *
cd ../../nirS 
rename 's/metagenome/D3b198A/' *
cd cluster 
rename 's/metagenome/D3b198A/' *
cd ../../nirK #NO NIRK
rename 's/metagenome/D3b198A/' *
cd cluster
rename 's/metagenome/D3b198A/' *

cd /Volumes/MyPassport/CariacoNpaper/xander_results/metaG/D3b198B/k45/recA
rename 's/metagenome/D3b198B/' *
cd cluster
rename 's/metagenome/D3b198B/' *
cd ../../nirS #NO NIRS
rename 's/metagenome/D3b198B/' *
cd cluster 
rename 's/metagenome/D3b198B/' *
cd ../../nirK 
rename 's/metagenome/D3b198B/' *
cd cluster
rename 's/metagenome/D3b198B/' *

cd /Volumes/MyPassport/CariacoNpaper/xander_results/metaG/D3b234A/k45/recA
rename 's/metagenome/D3b234A/' *
cd cluster
rename 's/metagenome/D3b234A/' *
cd ../../nirS
rename 's/metagenome/D3b234A/' *
cd cluster 
rename 's/metagenome/D3b234A/' *
cd ../../nirK 
rename 's/metagenome/D3b234A/' *
cd cluster
rename 's/metagenome/D3b234A/' *

cd /Volumes/MyPassport/CariacoNpaper/xander_results/metaG/D3b234B/k45/recA
rename 's/metagenome/D3b234B/' *
cd cluster
rename 's/metagenome/D3b234B/' *
cd ../../nirS
rename 's/metagenome/D3b234B/' *
cd cluster 
rename 's/metagenome/D3b234B/' *
cd ../../nirK 
rename 's/metagenome/D3b234B/' *
cd cluster
rename 's/metagenome/D3b234B/' *

cd /Volumes/MyPassport/CariacoNpaper/xander_results/metaG/D3b295A/k45/recA
rename 's/metagenome/D3b295A/' *
cd cluster
rename 's/metagenome/D3b295A/' *
cd ../../nirS
rename 's/metagenome/D3b295A/' *
cd cluster 
rename 's/metagenome/D3b295A/' *
cd ../../nirK #NO NIRK
rename 's/metagenome/D3b295A/' *
cd cluster
rename 's/metagenome/D3b295A/' *

cd /Volumes/MyPassport/CariacoNpaper/xander_results/metaG/D3b295B/k45/recA
rename 's/metagenome/D3b295B/' *
cd cluster
rename 's/metagenome/D3b295B/' *
cd ../../nirS
rename 's/metagenome/D3b295B/' *
cd cluster 
rename 's/metagenome/D3b295B/' *
cd ../../nirK #NO NIRK
rename 's/metagenome/D3b295B/' *
cd cluster
rename 's/metagenome/D3b295B/' *

cd /Volumes/MyPassport/CariacoNpaper/xander_results/metaG/D3b314A/k45/recA
rename 's/metagenome/D3b314A/' *
cd cluster
rename 's/metagenome/D3b314A/' *
cd ../../nirS
rename 's/metagenome/D3b314A/' *
cd cluster 
rename 's/metagenome/D3b314A/' *
cd ../../nirK #NO NIRK
rename 's/metagenome/D3b314A/' *
cd cluster
rename 's/metagenome/D3b314A/' *

cd /Volumes/MyPassport/CariacoNpaper/xander_results/metaG/D3b314B/k45/recA
rename 's/metagenome/D3b314B/' *
cd cluster
rename 's/metagenome/D3b314B/' *
cd ../../nirS
rename 's/metagenome/D3b314B/' *
cd cluster 
rename 's/metagenome/D3b314B/' *
cd ../../nirK #NO NIRK
rename 's/metagenome/D3b314B/' *
cd cluster
rename 's/metagenome/D3b314B/' *

```


- **NEXT STEP**- pull out taxonabun, coverage, and prot_aligned files for each sample. The first two will provide taxonomic identity and relative abundance of each functional gene taxon, plus total reads (to compare nirS:recA or nirK:recA). prot_aligned files should be checked in MEGA or something and provided to reviewers (they specifically requested the alignment- check to make sure it is ok and regions of gene are correct). framebot file tells which contigs were linked to which blast results



##**July 23rd 2018** 

### To do
- pull out taxonabun, coverage, and prot_aligned files for each sample and move into R

```
mkdir ResultsForR

find /Volumes/MyPassport/CariacoNpaper/xander_results/metaG -name "*_45_coverage.txt" -exec cp {} /Volumes/MyPassport/CariacoNpaper/xander_results/metaG/ResultsForR/ \;

find /Volumes/MyPassport/CariacoNpaper/xander_results/metaG -name "*_45_final_prot_aligned.fasta" -exec cp {} /Volumes/MyPassport/CariacoNpaper/xander_results/metaG/ResultsForR \;

find /Volumes/MyPassport/CariacoNpaper/xander_results/metaG -name "*_45_taxonabund.txt" -exec cp {} /Volumes/MyPassport/CariacoNpaper/xander_results/metaG/ResultsForR \;

find /Volumes/MyPassport/CariacoNpaper/xander_results/metaG -name "*_45_framebot.txt" -exec cp {} /Volumes/MyPassport/CariacoNpaper/xander_results/metaG/ResultsForR \;


```

- Copy this folder into new directory ("Xander_Results" in manuscript folders) locally and start a new R project
- Notes for R analysis are in script "CariacoNCycling.R"



##**July 24th 2018** 

### To do
- Get metaT analysis to work
- Complete figures:
	- PA vs FL for all nirS and nirK taxa
	- Depth profiles of nirs:RecA and nirK:recA
- Analyze alignments of nirS and nirK
- Before analyzing alignment, concatenate alignment fasta files

###Try to run metaT files for all 3 genes again###
- First prepare folders for each sample. Deleted old failed bloom filters

```
cd ~/Cariaco/metaT/May
for i in $(\ls -d *.fa)
do
	echo $i | cut -c1-7 | xargs -n 1 mkdir
done
for i in $(\ls -d *.fa | cut -c1-7)
do
	mv "$i""_merged.fa" $i
done

cd ~/Cariaco/metaT/Nov
for i in $(\ls -d *.fa)
do
	echo $i | cut -c1-7 | xargs -n 1 mkdir
done
for i in $(\ls -d *.fa | cut -c1-7)
do
	mv "$i""_merged.fa" $i
done


cp /home/ubuntu/tools/RDPTools/Xander_assembler/bin/xander_setenv.sh /home/ubuntu/Cariaco/metaT
cp /home/ubuntu/tools/RDPTools/Xander_assembler/bin/run_xander_skel.sh /home/ubuntu/Cariaco/metaT
cp /home/ubuntu/tools/RDPTools/Xander_assembler/bin/get_OTUabundance.sh /home/ubuntu/Cariaco/metaT

nano xander_setenv.sh
```

- Change parameters in xander_Setenv.sh
	- Largest metaT file is 13GB. Some are 5-6GB. Most are 3-6GB. A few files from Nov are very small (>.01GB)
	- ```K_SIZE``` 45
	- ```MIN_COUNT``` 1
	- ```FILTER_SIZE``` 37 (this is appropriate for a 16GB file.)
	- ```MAX_JVM_HEAP``` 32GB (corresponds to FILTER_SIZE of 38) 
	- ```THREADS``` 3. 
	- Leave rest of parameters as default

```
SEQFILE=$PWD/*.fa
WORKDIR=$PWD
REF_DIR=/home/ubuntu/tools/RDPTools/Xander_assembler 
JAR_DIR=/home/ubuntu/tools/RDPTools 
UCHIME=/home/ubuntu/tools/third_party_tools/uchime4.2.40_i86linux32 
HMMALIGN=/usr/local/bin/hmmalign

SAMPLE_SHORTNAME=$name
```


- Modify permsision ```chmod 755 xander_setenv.sh```
- Copy this version of the xander script and ```run_xander_skel``` into each sample folder:

```
echo ~/Cariaco/metaT/*/R* | xargs -n 1 cp xander_setenv.sh 
echo ~/Cariaco/metaT/*/R* | xargs -n 1 cp run_xander_skel.sh
```


- Run (not sure this will be succesful- wasn't able to get too many nirS transcripts from these files) **DIDNT WORK- BUILD BLOOM FILTER FAILED FOR ALL DIRECTORIES IN MAY. SEE WHAT HAPPENED**

```
tmux new -s Cariaco_metaT

cd ~/Cariaco/metaT/May
for name in $(\ls -d R2a103*) 
do
    cd ~/Cariaco/metaT/May/$name/
    ./run_xander_skel.sh xander_setenv.sh "build find search" "recA nirK nirS"; 
    cd ~/Cariaco/metaT/May
done

cd ~/Cariaco/metaT/Nov
for i in $(\ls -d R*) 
do
    cd ~/Cariaco/metaT/Nov/$i/
    ./run_xander_skel.sh xander_setenv.sh "build find search" "recA nirK nirS"; 
    cd ~/Cariaco/metaT/Nov
done
```
- this is not working with variable ```$name``` in xander_setenv.sh for sample short name. It does work when I run it individually with sample name manually input (```$PWD``` seems to be OK in script). 
- Just do it manually for now. 
	- Run ```./run_xander_skel.sh xander_setenv.sh "build find search" "recA nirK nirS"``` after modifying accompanying ```xander_setenv.sh``` to reflect current file and current directory

<u> Notes while running</u> (*Can modify individual ```xander_setenv.sh``` scripts in each sample directory while previous sample is runnin in tmux* ```tmux attach -t Cariaco_metaT```)

- May/R2a103A: Started ~3:30pm
	- Build bloom filter OK but find kmers for nirK failed and build contigs for nirS and recA failed
	- RERUN. THERE WAS A MEMORY ISSUE. CUSTOMIZE MEMORY PARAMETERS FOR THIS FILE Rerun @ ~4:50pm
		- .fa file size = 4.4GB
		- Use filter size of 35
		- bloom file size of 8GB
		- *Merge contigs still failed.
	- Copy parameters from metaG analysis and retry
		- filter size 36 and max jvm heap of 16 
		- **STILL UNABLE TO MERGE NIRS AND RECA CONTIGS** even though it recovers many recA contigs (1903)
- May/R2a103B 
	- Build bloom filter failed
	- Can't get this to work at all. Keep changing parameters, etc.
	- Maybe file is too small? only >0.2GB
	- ```grep ">" R2a103B_merged.fa | wc```
	- Contains 1,002,240 ">" whereas R2a103A contains 24,800,256
	- CAREFULLY RE PICK MEMORY PARAMETERS AND RERUN:
		- .fa file size = 0.19GB
		- Use filter_size of 32
		- bloom file size of 4GB
		- **STILL FAILED. FILE IS PROB JUST TOO SMALL**
- R2a198A Started running ~4:15pm
	- .fa file is 5.2GB
	- Use filter_size of 35
	- bloom file size of 8GB
	- **Build bloom filter worked but couldn't merge any contigs for nirS, recA.** Couldn't find any for nirK
- R2a198B 
	- .fa is 5.5GB
	- Use filter_size of 35
	- bloom file size of 8GB
	- **Only merged contigs for recA but at least it worked.**
- R2a234A
	- .fa is 2.7GB
	- Use ```FILTER_SIZE``` of 33
	- ```MAX_JVM_HEAP``` of 4GB
	- Left running in tmux at ~8:10pm and went home. Check from home
	- **No contigs merged again**
- R2a234B
	- .fa is 3.5GB
	- Use ```FILTER_SIZE``` of 34
	- ```MAX_JVM_HEAP``` of 4GB
	- **Only recA merged again**
	- Includes contig that matches to SAR202
- R2a295A
	- .fa file is 3.5GB
	- Use ```FILTER_SIZE``` of 34
	- ```MAX_JVM_HEAP``` of 4GB
	- **NOTHING MERGED**
- R2a295B
	- .fa file is 3.7GB
	- Use ```FILTER_SIZE``` of 34
	- ```MAX_JVM_HEAP``` of 4GB
	- **only recA- some diverse organisms including Nitrospina and oil-degrading strains** 
- R2a314A
	- .fa file is 5.1GB
	- Use ```FILTER_SIZE``` of 35
	- ```MAX_JVM_HEAP``` of 8GB
	- **Got 1 contig (~6 reads) or nirS! But recA didn't merge! :(**
- R2a314B
	- .fa file is 5.1GB
	- Use ```FILTER_SIZE``` of 35
	- ```MAX_JVM_HEAP``` of 8GB
	- **merged all 3!**
	- nirK from Nitrosococcus 

** I finished modifying every xander_setenv.sh file so each one if customized to it's respective sample (According to parameters and sample shortname below). Therefore I can run rest in a loop [see loop script after this list of parameters]

- R2b198A
	- .fa file is 5.4GB
	- Use ```FILTER_SIZE``` of 35
	- ```MAX_JVM_HEAP``` of 8GB
	- **only 1 recA contig- marine metagenome- blast match to Thioglobus** 
- R2b198B
	- .fa file is 4.3GB
	- Use ```FILTER_SIZE``` of 35
	- ```MAX_JVM_HEAP``` of 8GB
	- **only recA contigs** 
- R2b234A
	- .fa file is 5.9GB
	- Use ```FILTER_SIZE``` of 35
	- ```MAX_JVM_HEAP``` of 8GB
	- **no contigs merged**
- R2b234B
	- .fa file is 3.6GB
	- Use ```FILTER_SIZE``` of 34
	- ```MAX_JVM_HEAP``` of 4GB
	- **only recA contig** 
- R2b295A
	- .fa file is 4.4GB
	- Use ```FILTER_SIZE``` of 35
	- ```MAX_JVM_HEAP``` of 8GB
	- **no contigs merged**
- R2b295B
	- .fa file is 5.0GB
	- Use ```FILTER_SIZE``` of 35
	- ```MAX_JVM_HEAP``` of 8GB
	- **only recA contig** 
- R2b314A
	- .fa file is 3.8GB
	- Use ```FILTER_SIZE``` of 35
	- ```MAX_JVM_HEAP``` of 8GB
	- **only nirS contigs**
- R2b314B
	- .fa file is 4.1GB
	- Use ```FILTER_SIZE``` of 35
	- ```MAX_JVM_HEAP``` of 8GB
	- **recA and nirS contigs**
		- most diverse nirS so far and includes expression by Scalindua
		- recA in this sample and in many samples tends to be Methylophaga (methylotrophic gammaproteobacteria)- are these CH4 oxidizers (using NO3?)
		- Many of the hits to recA come back as methylotrophic/ oil-degrading gammaproteobacteria that can use alternate e- acceptors (NO3, DMSO)
- R3b103A
	- .fa file is 3.7GB
	- Use ```FILTER_SIZE``` of 35
	- ```MAX_JVM_HEAP``` of 8GB
	- **only recA contig**
- R3b103B
	- .fa file is 2.3GB
	- Use ```FILTER_SIZE``` of 34
	- ```MAX_JVM_HEAP``` of 4GB
	- **no contigs merged**
- Nov files: R2a148A
	- .fa file is 7.3GB
	- Use ```FILTER_SIZE``` of 36
	- ```MAX_JVM_HEAP``` of 16GB
	- **nothing merged**
- R2a148B
	- .fa file is 4.8GB
	- Use ```FILTER_SIZE``` of 35
	- ```MAX_JVM_HEAP``` of 8GB
	- **recA contigs were merged.**
- R2a200A
	- .fa file is 13.6GB
	- Use ```FILTER_SIZE``` of 37
	- ```MAX_JVM_HEAP``` of 32GB
	- **error in building bloom filter. Space issue. Use on step down for parameters**
		- Deleted old k45 folder and rerun with parameters:
		- ```FILTER_SIZE``` of 36
		- ```MAX_JVM_HEAP``` of 16GB
		- **worked. got recA and nirS contigs**
- R2a200B
	- .fa file is 6.6GB
	- Use ```FILTER_SIZE``` of 36
	- ```MAX_JVM_HEAP``` of 16GB
	- **nirS and recA contigs**
		- Only a Scalindua nirS
		- recA transcripts include recA from Brocadia!
- R2a237A
	- .fa file is 0.08GB *-->very small*
	- Use ```FILTER_SIZE``` of 32
	- ```MAX_JVM_HEAP``` of 4GB
	- **nothing merged**
- R2a237B
	- .fa file is 6.1GB 
	- Use ```FILTER_SIZE``` of 35
	- ```MAX_JVM_HEAP``` of 8GB
	- **nirS and recA contigs**
		- nirS trnascripts include one from Scalindua
		- recA also include Brocadia
- R2a247A
	- .fa file is 3.1GB 
	- Use ```FILTER_SIZE``` of 34
	- ```MAX_JVM_HEAP``` of 4GB
	- **only recA contigs**
		- more methylophaga
- R2a247B
	- .fa file is 1.4GB 
	- Use ```FILTER_SIZE``` of 33
	- ```MAX_JVM_HEAP``` of 4GB
	- **nirS and recA contigs**
		- Looks like mostly gammas in both cases
- R2a267A
	- .fa file is 5.5GB 
	- Use ```FILTER_SIZE``` of 35
	- ```MAX_JVM_HEAP``` of 8GB
	- **only recA contigs**
		- Burkholderia
- R2a267B
	- .fa file is 6.4GB 
	- Use ```FILTER_SIZE``` of 35
	- ```MAX_JVM_HEAP``` of 8GB 
		- **only recA contigs**
- R2b148B
	- .fa file is 5.9GB 
	- Use ```FILTER_SIZE``` of 35
	- ```MAX_JVM_HEAP``` of 8GB 
	- **nothing merged**
- R2b200A
	- .fa file is 0.006GB *-->very small*
	- Use ```FILTER_SIZE``` of 32
	- ```MAX_JVM_HEAP``` of 4GB 
	- **nothing**
- R2b200B
	- .fa file is 0.005GB *-->very small*
	- Use ```FILTER_SIZE``` of 32
	- ```MAX_JVM_HEAP``` of 4GB 
	- **nothing**
- R2b237A
	- .fa file is 0.015GB *-->very small*
	- Use ```FILTER_SIZE``` of 32
	- ```MAX_JVM_HEAP``` of 4GB 
	- **nothing**
- R2b237B
	- .fa file is 0.24GB *-->very small*
	- Use ```FILTER_SIZE``` of 32
	- ```MAX_JVM_HEAP``` of 4GB 
	- **produced nirS contigs**
		- proteobacterial
- R2b247A
	- .fa file is 1.3GB
	- Use ```FILTER_SIZE``` of 32
	- ```MAX_JVM_HEAP``` of 4GB
	- **nothing**
- R2b247B
	- .fa file is 0.003GB *-->very small*
	- Use ```FILTER_SIZE``` of 32
	- ```MAX_JVM_HEAP``` of 4GB
	- **nothing** 
- R2b267A
	- .fa file is 5.0GB
	- Use ```FILTER_SIZE``` of 35
	- ```MAX_JVM_HEAP``` of 8GB 
	- **recA only**
		- Burkholderia
- R2b267B
	- .fa file is 0.004GB *-->very small*
	- Use ```FILTER_SIZE``` of 32
	- ```MAX_JVM_HEAP``` of 4GB 
		- **nothing**



	

```
tmux new -s Cariaco_metaT_2

cd ~/Cariaco/metaT/May
for name in $(\ls -d R*b*) # already got through the R2a and R3a files
do
    cd ~/Cariaco/metaT/May/$name/
    ./run_xander_skel.sh xander_setenv.sh "build find search" "recA nirK nirS"; 
    cd ~/Cariaco/metaT/May
done

# still need to run the following
cd ~/Cariaco/metaT/Nov
for i in $(\ls -d R*) 
do
    cd ~/Cariaco/metaT/Nov/$i/
    ./run_xander_skel.sh xander_setenv.sh "build find search" "recA nirK nirS"; 
    cd ~/Cariaco/metaT/Nov
done
```




##**July 25th 2018** 

- Continue running Cariaco Nov metaT files (above) in tmux session
- Started formatting methods in manuscript
- Check out alignment files. Need to concatenate for each gene first:
	- I realized that the resulting fastas with the contigs from different samples are aligned with each themselves but not aligned across samples. (ie. They are only aligned within the same sample)
	
```
cat *_recA_45_final_prot_aligned.fasta > recA_contigs_prot_aligned.fasta
cat *_nirS_45_final_prot_aligned.fasta > nirS_contigs_prot_aligned.fasta
cat *_nirK_45_final_prot_aligned.fasta > nirK_contigs_prot_aligned.fasta
```

- First glance in Seaview, recA alignment looks pretty good. Can't find anything suspicious
- nirK alignment- much shorter (fewer contigs) but looks good. Can't find anything suspicious
- nirS alignment looks OK. There are large regions that are not covered by every contig. Still, for the regions that are covered, the alignment is pretty consistent. These are the 3 suspicious nirS contigs:

	>```metagenome_nirS_contig_440_contig_441```
PTLSTEDFERAKTLYFQRCAGCHGVLRKGATGSSLEPADTREKGQKRLERIIQLGTEGGMNNFDDIFSSEEISLLATYIQMEPPIPPEMSLAMMRERTKVYVQPQDYPTKPLHGRNWKNFFVVIERDAGQVAIIDGDTHEIVAHIDTGYAVH

	- [top blastp hit](https://www.ncbi.nlm.nih.gov/protein/WP_106671085.1?report=genbank&log$=prottop&blast_rank=1&RID=NGY9B178015) comes back as nitrite reductase from Marinobacter halophilus but I am still suspicious. It is only 80% identity (with 100% query cover) 


	>```metagenome_nirS_contig_34_contig_35```
FLVAVKELGQMWQVDYSDIENLDITKIDSAKYLHDGFFDPTGRYFQIAANASNKMVVVDTKTRKLEAMIDTTgKKPHPGPGANWDDPKHGPVGATVHLGTGMVTVWGNDPAGTPSKAWKIVREIETDGAGAFIRTHPNSKYVWADQVKNPEPEIQQSVQV

	- [top blastp hit](https://www.ncbi.nlm.nih.gov/protein/AFO85745.1?report=genbank&log$=prottop&blast_rank=1&RID=NGZ3H9NR01R) is multipsecies nitrite reductase with 99% identity (100% query cover) but second hit drops down to 78% identity (also a nitrite reductase). This hit is from Kirkpatrick denitrification/ anammox paper from Black Sea
	- in alignment, it looks like 1st suspicious one but their alignments don't look very much like other the contigs

	>```metagenome_nirS_contig_4_contig_5```
FLIAVKELGQMWQVDYKDLDNLKISQIDSAKFLHDGFFDPTGRYFQIAANASDKMVVVDTQDWNLEAMIDVDSKPHPGPGANWNDPKCGPVAGTTHLGIGTVTVWGNDPQGHPDNAWKICYEVETDGAGLFIRTHPASDYVWADQTKHPEPEVQQSVQV


	- [top blastp hit](https://www.ncbi.nlm.nih.gov/protein/AEM36244.1?report=genbank&log$=prottop&blast_rank=1&RID=NGZT21ES015) is "putative dissimilatory reductase" (84% identity and 100% coverage). Titles of study it comes from is specific to nirS. Other hits have similar stats and also say nitrite reductase (for example [2nd](https://www.ncbi.nlm.nih.gov/protein/WP_055676367.1?report=genbank&log$=prottop&blast_rank=2&RID=NGZT21ES015) hit is Labrenzia nitrite reductase)
	
- ACTUALLY once I add gaps in front of these contigs, they aligned pretty well. They were just by themselves in their own file so didn't align well with others in concatenated file. I am beginning to trust these contigs.


- TO DO- rename "metagenome" in each contig name with name of sample name so I can track where they came from once file is concatenated. 
	- Couldn't find easy way to do this so did it manually with find and replace in text files. Also deleted the "GC_RF" line which was full of X's showing where there was concensus I think

- Re-concatenated in terminal then aligned in seaview (parsimony) **CRAP I REALIZED SEAVIEW IS CUTTING THE NAMES ARGH TRY DOING ALIGNMENT IN MEGA SO IT RETAINS FULL NAME OF SEQ**
	- ```FileS1_nirS_contigs_aligned.out.fasta```
	- ```FileS2_nirK_contigs_aligned.out.fasta```
	- ```FileS3_recA_contigs_aligned.out.fasta```
	- All 3 files are stored in "FiguresTables" directory for manuscript as they are supplemental files for the submission
- OK even though Seaview cut the names, and I don't like MEGA's aligner so I am going to keep the alignment from Seaview
- Seaview also appended a number infront of each sequence name, in order, in front of every contig. So I actually just "grep"ed the contig names from the original fasta file, then put them in a list with the number in a text file as a "key". Can mess with this later if I want to change the actual name in the alignment fasta file. This is good enough for now.
	- ```FileS1_key.txt```
	- ```FileS2_key.txt```
	- ```FileS3_key.txt```

	

### MetaT file completed running through Xander
- Move results files to hard drive and start collecting results in folders

```
cd ~/Cariaco/metaT
mkdir Results
cd Results
mkdir May

cd /home/ubuntu/Cariaco/metaT/May
for i in $(\ls -d R*) 
do
    mkdir ~/Cariaco/metaT/Results/May/$i
    cp -r ~/Cariaco/metaT/May/$i/k45 ~/Cariaco/metaT/Results/May/$i
done

cd ~/Cariaco/metaT/Results
mkdir Nov

cd /home/ubuntu/Cariaco/metaT/Nov
for i in $(\ls -d R*) 
do
    mkdir ~/Cariaco/metaT/Results/Nov/$i
    cp -r ~/Cariaco/metaT/Nov/$i/k45 ~/Cariaco/metaT/Results/Nov/$i
done
```

- Transfer "results" folder from Instance to hard drive using FileZilla. Then delete from Instance to save space (there is still a copy of each k45 folder in the original sample folder)

- Working locally. Move files to ```Cariaco_NcyclingProject/xander_results``` 
	- DNA files are here too but they all start with D and metaT files start with R
	- For DNA files, I first made a "resultsforR" folder on Mypassport, copied files into there, renamed them, then moved folder onto Macbook. This time I can move them directly because they are already appropriately named

```
find /Volumes/MyPassport/CariacoNpaper/xander_results/metaT -name "*_45_coverage.txt" -exec cp {} /Users/admin/Google\ Drive/Wagner/Research/Manuscripts/Cariaco_NitrogenCycling/FEMSMicrobiologyEcology/SecondSubmission/CariacoNCycling_Rproject/Xander_results \;

find /Volumes/MyPassport/CariacoNpaper/xander_results/metaT -name "*_45_final_prot_aligned.fasta" -exec cp {} /Users/admin/Google\ Drive/Wagner/Research/Manuscripts/Cariaco_NitrogenCycling/FEMSMicrobiologyEcology/SecondSubmission/CariacoNCycling_Rproject/Xander_results \;

find /Volumes/MyPassport/CariacoNpaper/xander_results/metaT -name "*_45_taxonabund.txt" -exec cp {} /Users/admin/Google\ Drive/Wagner/Research/Manuscripts/Cariaco_NitrogenCycling/FEMSMicrobiologyEcology/SecondSubmission/CariacoNCycling_Rproject/Xander_results \;

find /Volumes/MyPassport/CariacoNpaper/xander_results/metaT -name "*_45_framebot.txt" -exec cp {} /Users/admin/Google\ Drive/Wagner/Research/Manuscripts/Cariaco_NitrogenCycling/FEMSMicrobiologyEcology/SecondSubmission/CariacoNCycling_Rproject/Xander_results \;

```




##**Aug 8th 2018** 

- Almost done formatting manuscripts
- For response to reviews I want to report how primers compare to contigs
- Need to pull out ```nucl_rep_seqs.fasta``` files which are the nucleotide files for the contigs.
- go into removeable hard drive and put files into new directory

```
cd /Volumes/MyPassport/CariacoNpaper/xander_results/metaG/ResultsForR


find /Volumes/MyPassport/CariacoNpaper/xander_results/metaG -name "*45_nucl_rep_seqs.fasta" -exec cp {} /Volumes/MyPassport/CariacoNpaper/xander_results/metaG/ResultsForR/ \;
```

- Concatenate into one file (don't care where contig is come from for now)

```
cat *_recA_45_nucl_rep_seqs.fasta > recA_nucl_rep_seqs.fasta
cat *_nirS_45_nucl_rep_seqs.fasta > nirS_nucl_rep_seqs.fasta
cat *_nirK_45_nucl_rep_seqs.fasta > nirK_nucl_rep_seqs.fasta
```

- Copy these 3 new files into local computer and check them out Arb
- Did probe matching. But only one match to forward primer and 18 to reverse.
- This is not fair to do because contig min length is 450bp (but gene is ~1500?) so there will be large regions of each contig that are just not covered by the primer. Explained this in repsonse to editor.