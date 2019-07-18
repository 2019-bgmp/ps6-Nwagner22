DAY OF CREATION: Thursday July, 11 2019         20190711

----------------------------------------------------------------------------------
Project Set 6  -- PS6
----------------------------------------------------------------------------------

(20190711)

Our goal with this assignment is to generate some measures of accuracy for whole genome assemblies and then use these measures to assay our success in several Velvet assemblies.

## PART 1 ##

-- ON MY MACHINE --

- colaborated with: Jared Galloway

parsing a FASTA file: contigs.fa

Using Python regular expressions, extract k-mer length of each contig (in { curly braces} below). In
addition, extract the k-mer coverage for the contig (in [brackets]). Assume a k-mer length of 49.
>\>NODE_11_length_{3717}\_cov\_[19.315845]
```
temp = re.search(r'.*_([0-9]{3,5}).*_([0-9]+\.[0-9]+)',line)
KMER_LENGTH.append(temp.group(1))
KMER_COVERAGE.append(temp.group(2))
```

How to calculte N50:
Finding N50 is the same as finding the weighted median. This is where %50 of the total information lies. *Not the middle value of the array*
```
CONTIG_LENGTH_ARRAY.sort()
temp_n50_sum = 0
for element in CONTIG_LENGTH_ARRAY:
    temp_n50_sum += element
    if(temp_n50_sum > (TOTAL_NT_BY_FORMULA/2)):
        N50 = element
        break
```
For calculating the distribution i used //100 integer division to put the different contig lengths into buckets

(20190717)

Unit testing:

I grabbed the first 21 records from contigs_converted.fa file and put them into the Unit_test.fa file. I then hand calculated TOTAL LENGTH OF GENOME ASSEMBLY, NUM CONTIG, MAX CONTIG LENGTH, MEAN CONTIG LENGTH, MEAN COVERAGE DEPTH and N50 as well as the buckets that each one should go into. I stored these hand calculated results into expected_results.txt
```
(base) [Nicholass-MacBook-Pro] [~/Desktop/BGMP/Repositories/ps6-Nwagner22] [Wed Jul 17 10:46:34] 
:: ./PS6_part1.py -f Unit_test.fa > unit_test_check.txt
(base) [Nicholass-MacBook-Pro] [~/Desktop/BGMP/Repositories/ps6-Nwagner22] [Wed Jul 17 10:46:59] 
:: diff unit_test_check.txt expected_results.txt 
20c20
< 6800	1
---
> 6800	1
\ No newline at end of file
```
Results were what I was expecting. One file just had a new line at the end where the other file did not.

## PART 2 ##

-- ON TALAPAS --

In this part, we’ll be assembling fosmids from goldeye fish (Hiodon alosoides) from quality trimmed Illumina PE100 reads.

Ran the commands below to install velvet using the bioconda channel
```
% conda activate bgmp_py3
% conda install velvet -c bioconda
% velvetg                           # displays the velvet usage
```
velvetg - de Bruijn graph construction, error removal and repeat resolution 

(20190712)

Files located here:
/projects/bgmp/shared/Bi621/800_3_PE5_interleaved.fq_1 
/projects/bgmp/shared/Bi621/800_3_PE5_interleaved.fq_2
/projects/bgmp/shared/Bi621/800_3_PE5_interleaved.fq.unmatched

a.  Our fosmid library is comprised of 50 fosmids, each (approximately) 40 Kb long. How many total nt should be in this fosmid library?
    50 * 40kb = 2000000 to NTs are expected in this library

Ran the command below multiple times to work out the bugs. I initially did not provide the $ to allow my variables to work. It took some trial and error to figure out how the sbatch was working.
```
/usr/bin/time -v velveth output_directory/ 31 -shortPaired -fastq -separate $fq_1 $fq_2 -short2 -fastq $fq_unmatched

Command being timed: "velveth output_directory/ 31 -shortPaired -fastq -separate /projects/bgmp/shared/Bi621/800_3_PE5_interleaved.fq_1 /projects/bgmp/shared/Bi621/800_3_PE5_interleaved.fq_2 -short2 -fastq /projects/bgmp/shared/Bi621/800_3_PE5_interleaved.fq.unmatched"
	User time (seconds): 2090.74
	System time (seconds): 1.30
	Percent of CPU this job got: 99%
	Elapsed (wall clock) time (h:mm:ss or m:ss): 34:52.36
	Maximum resident set size (kbytes): 732208
	Exit status: 0

```

(20190713)

Next, initially,  I added this velvetg command to the script. Without giving any flags past the output directory, velvet will output a statistics file 
```
/usr/bin/time -v velvetg output_directory/

Command being timed: "velvetg output_directory/"
User time (seconds): 79.80
System time (seconds): 1.21
Percent of CPU this job got: 99%
Elapsed (wall clock) time (h:mm:ss or m:ss): 1:21.16
Maximum resident set size (kbytes): 498000
Exit status: 0
```
This outputed 5 files to the output_directory/: contigs.fa, Graph, LastGraph, PreGraph, stats.txt
I converted contigs.fa to being two lines per record using my fa_converter.py script.         created file: contigs_converted.fa

b.  Running the PS6_part1.py script on the contigs_converted.fa file I got these results:
```
./PS6_part1.py -f 31_output_directory/contigs_converted.fa -k 31

KMER SIZE: 31
TOTAL LENGTH OF GENOME ASSEMBLY: 2202636
NUM CONTIG: 12051
MAX CONTIG LENGTH: 4059
MEAN CONTIG LENGTH: 200.057423874867
MEAN COVERAGE DEPTH: 37.68471118031698
N50:  299
```
c.  I plan to feed these results into my velvetg command using the flags: -ins_length 200 -exp_cov 37.68   (mean contig length and mean coverage depth respectively)
What I understood from the manual is that if you are trying to figure out values to feed into the velvetg flags, just run velvetg on the output_directory and
analyze the statistics that are outputed. This should allow me to refine the values i am giving for the different flags in order to get bettere data out of 
velvetg

The goal is to re-run velveth into a new output directory and preform the analysis again from above but with the new flags on velvetg to hopefully improve output.

slurm-9612248-PS6.errs
```
/usr/bin/time -v velveth new_output_directory/ 31 -shortPaired -fastq -separate $fq_1 $fq_2 -short2 -fastq $fq_unmatched
/usr/bin/time -v velvetg new_output_directory/ -ins_length 200 -exp_cov 37.68 
./fa_converter.py -i new_output_directory/contigs.fa contigs_converted.fa
./PS6_part1.py -f new_output_directory/contigs_converted.fa -k 31

Command being timed: "velveth new_output_directory/ 31 -shortPaired -fastq -separate /projects/bgmp/shared/Bi621/800_3_PE5_interleaved.fq_1 /projects/bgmp/shared/Bi621/800_3_PE5_interleaved.fq_2 -short2 -fastq /projects/bgmp/shared/Bi621/800_3_PE5_interleaved.fq.unmatched"
User time (seconds): 2096.88
System time (seconds): 1.34
Percent of CPU this job got: 99%
Elapsed (wall clock) time (h:mm:ss or m:ss): 34:58.46
Maximum resident set size (kbytes): 732140
Exit status: 0

Command being timed: "velvetg new_output_directory/ -ins_length 200 -exp_cov 37.68"
User time (seconds): 101.24
System time (seconds): 1.39
Percent of CPU this job got: 99%
Elapsed (wall clock) time (h:mm:ss or m:ss): 1:42.74
Maximum resident set size (kbytes): 553148
Exit status: 0

slurm-9612268-PS6.out

Running part1 script:
KMER SIZE: 31
TOTAL LENGTH OF GENOME ASSEMBLY: 2202981
NUM CONTIG: 10824
MAX CONTIG LENGTH: 35377
MEAN CONTIG LENGTH: 203.52743902439025
MEAN COVERAGE DEPTH: 38.5063806568736
N50:  406
```

The next step was to run velveth and velvetg with kmer sizes of 41 and 49 so I edited my script to run though both amd do the part 1 analysis as well:
I bumped up my script run time to 2 hours because each one of the kmer sizes below will take over 30 minutes each (This was not needed becasue I gave it 8 cpus and it finished in under 3 minutes)

slurm-9612328-PS6.out
```
/usr/bin/time -v velveth 41_output_directory/ 41 -shortPaired -fastq -separate $fq_1 $fq_2 -short2 -fastq $fq_unmatched
/usr/bin/time -v velvetg 41_output_directory/ -ins_length 200 -exp_cov 37.68 
./fa_converter.py -i 41_output_directory/contigs.fa -o 41_output_directory/contigs_converted.fa
./PS6_part1.py -f 41_output_directory/contigs_converted.fa -k 41

Command being timed: "velveth 41_output_directory/ 41 -shortPaired -fastq -separate /projects/bgmp/shared/Bi621/800_3_PE5_interleaved.fq_1 /projects/bgmp/shared/Bi621/800_3_PE5_interleaved.fq_2 -short2 -fastq /projects/bgmp/shared/Bi621/800_3_PE5_interleaved.fq.unmatched"
User time (seconds): 53.02
System time (seconds): 1.66
Percent of CPU this job got: 336%       (FOR FUTURE ONLY GIVE 3 CPUs)
Elapsed (wall clock) time (h:mm:ss or m:ss): 0:16.26
Maximum resident set size (kbytes): 718984
Exit status: 0

Running part1 script:
KMER SIZE: 41
TOTAL LENGTH OF GENOME ASSEMBLY: 1755889
NUM CONTIG: 4179
MAX CONTIG LENGTH: 36884
MEAN CONTIG LENGTH: 420.1696578128739
MEAN COVERAGE DEPTH: 38.93400042546067
N50:  2543


/usr/bin/time -v velveth 49_output_directory/ 49 -shortPaired -fastq -separate $fq_1 $fq_2 -short2 -fastq $fq_unmatched
/usr/bin/time -v velvetg 49_output_directory/ -ins_length 200 -exp_cov 37.68 
./fa_converter.py -i 49_output_directory/contigs.fa -o 49_output_directory/contigs_converted.fa
./PS6_part1.py -f 49_output_directory/contigs_converted.fa -k 49

Command being timed: "velveth 49_output_directory/ 49 -shortPaired -fastq -separate /projects/bgmp/shared/Bi621/800_3_PE5_interleaved.fq_1 /projects/bgmp/shared/Bi621/800_3_PE5_interleaved.fq_2 -short2 -fastq /projects/bgmp/shared/Bi621/800_3_PE5_interleaved.fq.unmatched"
User time (seconds): 45.09
System time (seconds): 1.74
Percent of CPU this job got: 311%
Elapsed (wall clock) time (h:mm:ss or m:ss): 0:15.02
Maximum resident set size (kbytes): 688656
Exit status: 0

Running part1 script:
KMER SIZE: 49
TOTAL LENGTH OF GENOME ASSEMBLY: 1606161
NUM CONTIG: 2027
MAX CONTIG LENGTH: 41137
MEAN CONTIG LENGTH: 792.3833251110015
MEAN COVERAGE DEPTH: 49.084100214109434
N50:  4395
```

With a k-mer size of 49, adjust the coverage cutoff to 20x, 60x, and ‘auto’. Again, assay
your results with your code.
changed back to 1 cpu here:

ALL 3 OF THE NEXT SECTIONS ARE INCORRECT. (20X,60X,AUTO) I DID NOT USE THE CORRECT FLAG FOR COVERAGE CUTOFF. I INITALLY JUST REPLACED THE EXP-COV FLAG
REDOING BELOW the tabbed out sections

	*-- 20x --

	*slurm-9612545-PS6.out
	```
	<!-- Final graph has 5141 nodes and n50 of 1767, max 36218, total 1522688, using 1976351/2567781 reads
	Command being timed: "velvetg 49_output_directory/ -ins_length 200 -exp_cov 20"
	User time (seconds): 7.56
	System time (seconds): 0.86
	Percent of CPU this job got: 99%
	Elapsed (wall clock) time (h:mm:ss or m:ss): 0:08.44
	Maximum resident set size (kbytes): 512312
	Exit status: 0

	Running part1 script:
	TOTAL LENGTH OF GENOME ASSEMBLY: 1607125
	NUM CONTIG: 2535
	MAX CONTIG LENGTH: 36266
	MEAN CONTIG LENGTH: 633.974358974359
	MEAN COVERAGE DEPTH: 46.447045155424036
	N50:  1689 -->
	```
	*-- 60x --

	slurm-9612995-PS6.out
	```
	<!-- Final graph has 4082 nodes and n50 of 7285, max 41089, total 1558786, using 1980825/2567781 reads
	Command being timed: "velvetg 49_output_directory/ -ins_length 200 -exp_cov 60"
	User time (seconds): 7.75
	System time (seconds): 0.88
	Percent of CPU this job got: 98%
	Elapsed (wall clock) time (h:mm:ss or m:ss): 0:08.77
	Maximum resident set size (kbytes): 512312
	Exit status: 0

	Running part1 script:
	TOTAL LENGTH OF GENOME ASSEMBLY: 1604914
	NUM CONTIG: 1694
	MAX CONTIG LENGTH: 41137
	MEAN CONTIG LENGTH: 947.4108618654074
	MEAN COVERAGE DEPTH: 46.41071139197155
	N50:  6404 -->
	```
	*-- auto --

	slurm-9612997-PS6.out
	```
	<!-- Final graph has 1632 nodes and n50 of 4430, max 36250, total 1167088, using 1903968/2567781 reads
	Command being timed: "velvetg 49_output_directory/ -ins_length 200 -exp_cov auto"
	User time (seconds): 6.97
	System time (seconds): 0.88
	Percent of CPU this job got: 98%
	Elapsed (wall clock) time (h:mm:ss or m:ss): 0:07.94
	Maximum resident set size (kbytes): 512340
	Exit status: 0

	Running part1 script:
	TOTAL LENGTH OF GENOME ASSEMBLY: 1196284
	NUM CONTIG: 835
	MAX CONTIG LENGTH: 36298
	MEAN CONTIG LENGTH: 1432.6754491017964
	MEAN COVERAGE DEPTH: 74.86730145748498
	N50:  4313 -->
	```

(20190714)

I ran it in one bash script this time but I will seperate the results below.

-- 20x --

slurm-9613298-PS6.out
```
/usr/bin/time -v velvetg 49_output_directory/ -ins_length 200 -exp_cov 37.68 -cov_cutoff 20
./fa_converter.py -i 49_output_directory/contigs.fa -o 49_output_directory/20x_contigs_converted.fa
./PS6_part1.py -f 49_output_directory/20x_contigs_converted.fa -k 49
Final graph has 1461 nodes and n50 of 4397, max 36164, total 1090492, using 1868957/2567781 reads
	Command being timed: "velvetg 49_output_directory/ -ins_length 200 -exp_cov 37.68 -cov_cutoff 20"
	User time (seconds): 7.30
	System time (seconds): 0.84
	Percent of CPU this job got: 98%
	Elapsed (wall clock) time (h:mm:ss or m:ss): 0:08.30
	Maximum resident set size (kbytes): 512360
	Exit status: 0

Running part1 script:
KMER SIZE: 49
TOTAL LENGTH OF GENOME ASSEMBLY: 1116365
NUM CONTIG: 748
MAX CONTIG LENGTH: 36212
MEAN CONTIG LENGTH: 1492.466577540107
MEAN COVERAGE DEPTH: 81.0097707740642
N50:  4318
```

-- 60x --

slurm-9613298-PS6.out
```
/usr/bin/time -v velvetg 49_output_directory/ -ins_length 200 -exp_cov 37.68 -cov_cutoff 60
./fa_converter.py -i 49_output_directory/contigs.fa -o 49_output_directory/60x_contigs_converted.fa
./PS6_part1.py -f 49_output_directory/60x_contigs_converted.fa -k49

Final graph has 706 nodes and n50 of 2035, max 6503, total 385487, using 1139678/2567781 reads
	Command being timed: "velvetg 49_output_directory/ -ins_length 200 -exp_cov 37.68 -cov_cutoff 60"
	User time (seconds): 5.47
	System time (seconds): 0.79
	Percent of CPU this job got: 99%
	Elapsed (wall clock) time (h:mm:ss or m:ss): 0:06.30
	Maximum resident set size (kbytes): 512352
	Exit status: 0

Running part1 script:
KMER SIZE: 49
TOTAL LENGTH OF GENOME ASSEMBLY: 401918
NUM CONTIG: 426
MAX CONTIG LENGTH: 6551
MEAN CONTIG LENGTH: 943.4694835680751
MEAN COVERAGE DEPTH: 122.88452067136153
N50:  2044
```

-- auto --

slurm-9613298-PS6.out
```
/usr/bin/time -v velvetg 49_output_directory/ -ins_length 200 -exp_cov 37.68 -cov_cutoff auto
./fa_converter.py -i 49_output_directory/contigs.fa -o 49_output_directory/auto_contigs_converted.fa
./PS6_part1.py -f 49_output_directory/auto_contigs_converted.fa -k 49

Final graph has 1514 nodes and n50 of 4845, max 36250, total 1169500, using 1905554/2567781 reads
	Command being timed: "velvetg 49_output_directory/ -ins_length 200 -exp_cov 37.68 -cov_cutoff auto"
	User time (seconds): 6.47
	System time (seconds): 0.77
	Percent of CPU this job got: 99%
	Elapsed (wall clock) time (h:mm:ss or m:ss): 0:07.27
	Maximum resident set size (kbytes): 512420
	Exit status: 0

Running part1 script:
KMER SIZE: 49
TOTAL LENGTH OF GENOME ASSEMBLY: 1195978
NUM CONTIG: 769
MAX CONTIG LENGTH: 36298
MEAN CONTIG LENGTH: 1555.2379713914174
MEAN COVERAGE DEPTH: 77.7072203485045
N50:  4564
```

--  -ins_length 500 -exp_cov 37.68 -cov_cutoff auto --

slurm-9613300-PS6.out
```

Final graph has 1502 nodes and n50 of 5991, max 41765, total 1293124, using 1905609/2567781 reads
	Command being timed: "velvetg 49_output_directory/ -ins_length 500 -exp_cov 37.68 -cov_cutoff auto"
	User time (seconds): 7.20
	System time (seconds): 0.85
	Percent of CPU this job got: 99%
	Elapsed (wall clock) time (h:mm:ss or m:ss): 0:08.07
	Maximum resident set size (kbytes): 512356
	Exit status: 0

Running part1 script:
KMER SIZE: 49
TOTAL LENGTH OF GENOME ASSEMBLY: 1319148
NUM CONTIG: 759
MAX CONTIG LENGTH: 41813
MEAN CONTIG LENGTH: 1738.00790513834
MEAN COVERAGE DEPTH: 77.40971555731231
N50:  5810
```

slurm-9613301-PS6.out, slurm-9613302-PS6.out, slurm-9613303-PS6.out were re-runs of the three velvetg commands for each of the coverage cutoff limits to be able to store the individual graph information in seperate folder

(20190717)

Graphing: PS6_graphing.py
```
K = 31, 41, 49:
K31_distribution.jpg
K41_distribution.jpg
K49_distribution.jpg

coverage cutoff = 20x, 60x auto
20x_distribution.jpg
60x_distribution.jpg
auto_distribution.jpg

ins_length = 500 cov_cutoff = auto
ins_500_cov_auto_distribution.jpg
```


## Part 3 ##

Part 3 – Questions
1. Describe how the assembly changes with different k-mer values using the assembly statistics you have collected. How does the contig length distribution change?

	- The assembly becomes more contiguous as the kmer size increases. For k=31, 41, and 49 the total number of contigs dropped from 10824, 4179, 2027 respectively.
	MAX CONTIG LENGTH, MEAN CONTIG LENGT and MEAN COVERAGE DEPTH all increase as the kmer size increases. Becasue of all of this, N50 also increased as the kmer size 
	increased. The distribution for all three graphs is very heavy to the side of the smaller contig lengths. the k41 and k49 graphs start to become more elongated as the max contig size goes up, which means that there are a lot less contigs in the 0-99 and 100-199 buckets.

2. How does an increased coverage cutoff affect the assembly? What is happening to the de Bruijin graph when you change the value of this parameter? How does velvet calculate its value for ‘auto’?

	- The increased coverage cutoff affects the assembly in a positive manner up until it exceeds the expected coverage that we also feed into velvetg. When the coverage
	cutoff exceeds the expected coverage value you give it, the total length of genome, number of contigs, max contig length, mean contig length, and N50 all drop. The 
	mean coverage depth skyrockets to 122.88 though which does not make sense. With a higher coverage cutoff, the de Bruijin graph becomes way more cluttered with edges
	and nodes. In a de Bruijin graph, the nodes are the k-1mers (the overlap) and the edges are the kmers themselves. With a higher coverage cutoff, there will be way more
	overlap between the kmers and therefore a lot more nodes and connections between them. It would be a lot harder to traverse the graph and figure out what connections
	are supposed to be there with all of the conveluted data. Therefore a lower coverage cutoff means a simpler and easier to put together graph, but too low of a cutoff
	and you will be losing too much information to put together a correct graph. It is all about finding the happy medium which is better calculated than using auto. velvet 
	calculates its auto cov_cutoff by taking half the length weighted median contig coverage depth.

3. How does increasing minimum contig length affect your contig length distribution and N50?

	- Increasing minimum contig length produced the best distribution and longest N50 values. This must mean that the mean insert length I originally fed into the 
	velvetg was too small for the given data.