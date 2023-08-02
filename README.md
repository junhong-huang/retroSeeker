# retroSeeker

retroSeeker: A computational software for identifying retrotransposons from pairwise alignment data(.net).

## Input data:<BR>

Pairwise alignment data(.net) can be download from UCSC: https://hgdownload.soe.ucsc.edu/downloads.html#human
or customized by lastz software described in https://hgdownload.soe.ucsc.edu/goldenPath/hg38/vsMm10/

Usage:<BR>
---------

Usage:  retroSeeker [options] --fa <fasta file> --fai <fai file> --net <net file><BR>

-f/--fa <file>      : genome file, fasta format [Required]<BR>
-F/--fai <file>     : fai file of genome, you can use samtools faidx to generate [Required]<BR>
-n/--net <file>     : Net file, Net Format [Required]<BR>
-b/--bed [file]     : bed file, bed Format, scan the bed content instead of the net file [option]<BR>
-o/--output <file>  : main retroSeeker output file<BR>
-O/--output2 <file> : output the gap and fill region (used to identified the retrogene) in the net file [option]<BR>
-v/--verbose        : verbose information<BR>
-V/--version        : retroSeeker version<BR>
-h/--help           : help informations<BR>
-S/--self           : self-chain net type. default is different normal-chain net<BR>
-x/--RNA5endUp      : extend length of upStream of 5'-end of the RNA direction to search tsd [default=30]<BR>
-y/--RNA5endDown    : extend length of dpStream of 3'-end of the RNA direction to search tsd [default=30]<BR>
-z/--RNA3endUp      : extend length of upStream of 5'-end of the RNA direction to search tsd [default=30]<BR>
-w/--RNA3endDown    : extend length of dpStream of 3'-end of the RNA direction to search tsd [default=30]<BR>
-m/--minMatchLen    : minimum TSD length [default>=5]<BR>
-M/--minPolyaLen    : minimum polya length [default>=5]<BR>
-P/--maxPolyaLen    : maximum polya length [default<=50]<BR>
-g/--minGenebodyLen : minimum genebody length [default>=20]<BR>
-l/--minLen         : minimum length for retrogene [default>=50]<BR>
-L/--maxLen         : maximum length for retrogene [default<=100000]<BR>
-s/--score          : minimum score for retrogene [default>=10]<BR>


Installation:<BR>
---------

Download retroSeeker-1.0.tar.gz from https://github.com/junhong-huang/retroSeeker/releases ; unpack it, and make:<BR>
tar -xzvf retroSeeker-1.0.tar.gz<BR>
cd retroSeeker-1.0<BR>
make<BR>

System requirements:<BR>
---------

Operating system: retroSeeker is designed to run on POSIX-compatible platforms, including UNIX, Linux and Mac OS/X. We have tested  most extensively on Linux and MacOS/X because these are the machines we develop on.<BR>
Compiler: The source code is compiled with  the C++ compiler g++. We test the code using the g++ compilers.<BR>


Run retroSeeker:<BR>
---------

retroSeeker --fa hg38.fa --fai hg38.fa.fai --bed retroACA.demo.bed --minMatchLen 5 --minPolyaLen 5 --minLen 60 --maxLen 100000 --RNA5endUp 30 --RNA5endDown 30 --RNA3endUp 30 --RNA3endDown 30<BR>

Output:<BR>
---------

#chrom	chromStart	chromEnd	name	score	strand	tsdScore	tsdLength	retroLength	retroSeq	tsd5Seq	tsdPair	tsd3Seq	polyaStart	polyaEnd	polyaLength	polyaScore	polyaSeq
chr1	53770991	53771163	retroSeeker_1	27	-	15	10	172	GCAAGGAGAAGGGCATACCCGTAGACCTTGCCTGACTGTGCTCATGTCCAGGCAGGGGGGACATTGTATTCGAGATTAATTTGAAGTTCCTGCCAGCTTTATCCAGCTTAATCAGTGGCTGGATAAATAGCAGGACTGTAACATTCCCCTGGGGGAAAAAAGGCAAGAAGAA	GCAAGGAGAA	|||||.||||	GCAAGAAGAA	156	161	6	12	AAAAAA

## How to cite:<BR>

Huang J et.al, New retrotransposon classes and specificity revealed by retroSeeker, Advanced Biotechnology, 2023


Contact :<BR>
---------

Jun-Hong Huang (huangjh278@mail2.sysu.edu.cn)
RNA Information Center, State Key Laboratory for Biocontrol, Sun Yat-sen University, Guangzhou 510275, P. R. China