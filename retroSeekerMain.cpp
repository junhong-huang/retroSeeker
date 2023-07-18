/******************************************************************************
 * retroSeeker: A computational software for seeking                          *
 * the Retrotransposon/Retroelements in the net File                          *
 * Author: Jianhua Yang                                                       *
 * Email: yangjh7@mail.sysu.edu.cn                                            *
 * Copyright: School of Life Sciences, Sun Yat-sen University                 *
 * Create Date: 2023/03/18                                                    *
 *****************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <assert.h>
#include <math.h>
#include <stdint.h>
#include <getopt.h>
#include <map>
#include <algorithm>
#include <ios>
#include <iostream>
#include <string>
#include <ostream>
#include <fstream>
#include <iomanip>
#include <locale>
#include <sstream>
#include <vector>

using namespace std;

#include "bioUtils.h"
#include "faiFile.h"
#include "bedFile.h"

#include "retroSeeker.h"

void usage(void);
char version[] = "retroSeeker version 2 [2023/6/16]\n";

int main(int argc, char **argv) {
  int i            = 0;
  int c            = 0;
  char *genomeFile = NULL;
  FILE *genomefp   = NULL;
  char *faiFile    = NULL;
  FILE *faifp      = NULL;
  char *netFile    = NULL;
  FILE *netfp      = NULL;
  char *bedFile    = NULL;
  FILE *bedfp      = NULL;
  FILE *outfp      = NULL;
  FILE *outfp2     = NULL;
  char *outfile    = NULL;
  char *outfile2    = NULL;
  int showVersion  = 0;
  int showHelp     = 0;
  parameter paraInfo;

  if (argc == 1) {
    usage();
  }

  const char *shortOptions = "vhVSoO:i:a:m:f:F:l:n:M:L:s:b:e:x:y:z:w:P:g:";

  const struct option longOptions[] = {
      {"verbose",   no_argument,       NULL, 'v'},
      {"help",      no_argument,       NULL, 'h'},
      {"version",   no_argument,       NULL, 'V'},
      {"self",      no_argument,       NULL, 'S'},
      {"output",    required_argument, NULL, 'o'},
      {"output2",    required_argument, NULL, 'O'},
      {"fa",        required_argument, NULL, 'f'},
      {"fai",       required_argument, NULL, 'F'},
      {"net",       required_argument, NULL, 'n'},
      {"bed",       required_argument, NULL, 'b'},
      {"extend",    required_argument, NULL, 'e'},
      {"RNA5endUp",    required_argument, NULL, 'x'},
      {"RNA5endDown",    required_argument, NULL, 'y'},
      {"RNA3endUp",    required_argument, NULL, 'z'},
      {"RNA3endDown",    required_argument, NULL, 'w'},
      {"minMatchLen",   required_argument, NULL, 'm'},
      {"minPolyaLen", required_argument, NULL, 'M'},
      {"maxPolyaLen", required_argument, NULL, 'P'},
      {"minGenebodyLen", required_argument, NULL, 'g'},
      {"minLen",   required_argument, NULL, 'l'},
      {"maxLen",   required_argument, NULL, 'L'},
      {"score",     required_argument, NULL, 's'},
      {NULL,        0,                 NULL, 0  }, /* Required at end of array. */
  };

  paraInfo.verbose     = 0;
  paraInfo.extend   = 30;
  paraInfo.RNA5endUp   = 30;
  paraInfo.RNA5endDown   = 30;
  paraInfo.RNA3endUp   = 30;
  paraInfo.RNA3endDown   = 30;
  paraInfo.minMatchLen = 5;/*min-tsd*/
  paraInfo.minPolyaLen = 5;
  paraInfo.maxPolyaLen = 50;
  paraInfo.minGenebodyLen = 20;
  paraInfo.minLen      = 50;
  paraInfo.maxLen      = 100000;
  paraInfo.netType     = 0;
  paraInfo.minScore    = 10;
  paraInfo.bedExist    = 0;
  paraInfo.netExist    = 0;

  while ((c = getopt_long(argc, argv, shortOptions, longOptions, NULL)) >= 0) {
    switch (c) {
    case 'v':
      paraInfo.verbose = 1;
      break;
    case 'h':
      showHelp = 1;
      break;
    case 'V':
      showVersion = 1;
      break;
    case 'S':
      paraInfo.netType = 1;
      break;
    case 'o':
      outfile = optarg;
      break;
    case 'O':
      outfile2 = optarg;
      break;
    case 'f':
      genomeFile = optarg;
      break;
    case 'F':
      faiFile = optarg;
      break;
    case 'n':
      netFile = optarg;
      paraInfo.netExist = 1;
      break;
    case 'b':
      bedFile = optarg;
      paraInfo.bedExist = 1;
      break;
    case 'e':
      paraInfo.extend = atoi(optarg);
      break;
    case 'x':
      paraInfo.RNA5endUp = atoi(optarg);
      break;
    case 'y':
      paraInfo.RNA5endDown  = atoi(optarg);
      break;
    case 'z':
      paraInfo.RNA3endUp  = atoi(optarg);
      break;
    case 'w':
      paraInfo.RNA3endDown  = atoi(optarg);
      break;     
    case 'm':
      paraInfo.minMatchLen = atoi(optarg);
      break;
    case 'M':
      paraInfo.minPolyaLen = atoi(optarg);
      break;
    case 'P':
      paraInfo.maxPolyaLen = atoi(optarg);
      break;
    case 'g':
      paraInfo.minGenebodyLen = atoi(optarg);
      break;
    case 'l':
      paraInfo.minLen = atoi(optarg);
      break;
    case 'L':
      paraInfo.maxLen = atoi(optarg);
      break;
    case 's':
      paraInfo.minScore = atoi(optarg);
      break;
    case '?':
      showHelp = 1;
      break;
    default:
      usage();
    }
  }

  if (showHelp) {
    usage();
    exit(1);
  }
  if (showVersion) {
    fprintf(stderr, "%s", version);
    exit(1);
  }

  if (genomeFile == NULL) {
    fprintf(stderr, "Please set the --fa option\n");
    exit(1);
  }
  genomefp = (FILE *)fopen(genomeFile, "r");
  if (genomefp == NULL) {
    fprintf(stderr, "ERROR: Can't open %s\n", genomeFile);
    fprintf(stderr, "Please set the --fa option\n");
    exit(1);
  }

  if (faiFile == NULL) {
    fprintf(stderr, "Please set the --fai option\n");
    exit(1);
  }
  faifp = (FILE *)fopen(faiFile, "r");
  if (faifp == NULL) {
    fprintf(stderr, "ERROR: Can't open %s\n", faiFile);
    fprintf(stderr, "Please set the --fai option\n");
    exit(1);
  }

  if (bedFile != NULL){
    fprintf(stderr, "#Scan bed file...\n");
    bedfp = (FILE *)fopen(bedFile, "r");
  }
  if (netFile != NULL){
    fprintf(stderr, "#Scan net file...\n");
    netfp = (FILE *)fopen(netFile, "r");
  }
  if (netFile == NULL && bedFile == NULL) {
    fprintf(stderr, "ERROR: Can't open %s and %s\n", netFile, bedFile);
    fprintf(stderr, "Please set the --net or --bed option\n");
    exit(1);
  }

  if (outfile == NULL) {
    outfp = stdout;
  } 
  else {
    outfp = (FILE *)fopen(outfile, "w");
    if (outfp == NULL) {
      fprintf(stderr, "ERROR: Can't open the output file %s\n", outfile);
      usage();
    }
  }

  if (outfile2 == NULL) {
    outfp2 = stdout;
  } 
  else {
    outfp2 = (FILE *)fopen(outfile2, "w");
    if (outfp2 == NULL) {
      fprintf(stderr, "ERROR: Can't open the output file2 %s\n", outfile2);
      usage();
    }
  }

  searchRetroGeneInNetFile(&paraInfo, genomefp, faifp, netfp, outfp, bedfp, outfp2);

  fflush(outfp);
  fflush(outfp2);
  fclose(genomefp);
  fclose(faifp);
  if (netFile != NULL){
    fclose(netfp);
  }
  if (bedfp != NULL){
    fclose(bedfp);
  }

  return 0;
}

void usage(void) {
  fprintf(
      stderr, "%s",
      "Usage:  retroSeeker2 [options] --fa <fasta file> --fai <fai file> --net <net file> --bed [bed file]\n\
[options]\n\
-f/--fa <file>      : genome file, fasta format [Required]\n\
-F/--fai <file>     : fai file of genome, you can use samtools faidx to generate [Required]\n\
-n/--net <file>     : Net file, Net Format [Required]\n\
-b/--bed [file]     : bed file, bed Format, scan the bed content instead of the net file [option]\n\
-o/--output <file>  : main retroSeeker output file\n\
-O/--output2 <file> : output the gap and fill region (used to identified the retrogene) in the net file [option]\n\
-v/--verbose        : verbose information\n\
-V/--version        : retroSeeker version\n\
-h/--help           : help informations\n\
-S/--self           : self-chain net type. default is different normal-chain net\n\
-x/--RNA5endUp      : extend length of upStream of 5`-end of the RNA direction to search tsd [default=30]\n\
-y/--RNA5endDown    : extend length of dpStream of 3`-end of the RNA direction to search tsd [default=30]\n\
-z/--RNA3endUp      : extend length of upStream of 5`-end of the RNA direction to search tsd [default=30]\n\
-w/--RNA3endDown    : extend length of dpStream of 3`-end of the RNA direction to search tsd [default=30]\n\
-m/--minMatchLen    : minimum TSD length [default>=5]\n\
-M/--minPolyaLen    : minimum polya length [default>=5]\n\
-P/--maxPolyaLen    : maximum polya length [default<=50]\n\
-g/--minGenebodyLen : minimum genebody length [default>=20]\n\
-l/--minLen         : minimum length for retrogene [default>=50]\n\
-L/--maxLen         : maximum length for retrogene [default<=100000]\n\
-s/--score          : minimum score for retrogene [default>=10]\n\
");
  exit(1);
}
