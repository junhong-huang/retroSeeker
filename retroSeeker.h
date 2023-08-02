/******************************************************************************
 * retroSeeker: A computational software for seeking                          *
 * the Retrotransposon/Retroelements in the net File                          *
 * Author: Jianhua Yang                                                       *
 * Email: yangjh7@mail.sysu.edu.cn                                            *
 * Copyright: School of Life Sciences, Sun Yat-sen University                 *
 * Create Date: 2023/03/18                                                    *
 *****************************************************************************/

#ifndef retroSeeker_HEAD_H
#define retroSeeker_HEAD_H
#include "fasta.h"

#ifndef MIN
#define MIN(x, y) (x) < (y) ? (x) : (y)
#endif

#ifndef MAX
#define MAX(x, y) (x) > (y) ? (x) : (y)
#endif


#define GAP_CONT -3.0
#define MATCH    2.0
#define GAP_OPEN -3.0
#define MISMATCH -3.0



#define POLYA_MATCH    2.0
#define POLYA_MISMATCH -3.0
#define POLYA_BONUS    2.0

#define MIN_GENEBODY_LEN 20
#define MAX_POLYA_LEN 50

#define MIN_TSD_LEN_FOR_GREEDY_POLYA_SEARCH   7
#define TSD_EXTEND_LEN 30/* 30*2: --extend*/

struct parameterInfo {
  int verbose;
  int extend;
  int RNA5endUp;
  int RNA5endDown;
  int RNA3endUp;
  int RNA3endDown;
  int minMatchLen;
  int minLen;
  int maxLen;
  int minPolyaLen;
  int maxPolyaLen;
  int minGenebodyLen;
  int netType;
  int bedExist;
  int netExist;
  int minScore;
};

typedef struct parameterInfo parameter;

struct polyaInfo {
  int polyaScore;
  int polyaLen;
  int start;
  int end;
  char *polyaSeq;
};

typedef struct polyaInfo polyaStruct;

void searchRetroGeneInNetFile(parameter *paraInfo, FILE *genomefp, FILE *faifp,
                              FILE *netfp, FILE *outfp, FILE *bedfp, FILE *outfp2);

void identifyRetroGenes(parameter *paraInfo, FILE *genomefp, faidxMap &faiHash,
                        FILE *netfp, FILE *outfp, FILE *bedfp, FILE *outfp2);

int identifyOneRetroGene(parameter *paraInfo, FILE *genomefp, faidxMap &faiHash,
                         CBed6 *fillBed, int chromLen, FILE *outfp);

int tsdSeeker(parameter *paraInfo, char *retroSeq, char *tsd5Seq, char *tsd3Seq,
              int chromStart, int chromEnd, CBed6 *fillBed, FILE *outfp);

int get5pTrimNum(char *seq, int tsdOffset);

int getTsdStart(char *retroSeq, char *tSeq, int polyaEnd);

/*polyaStruct *findPolyA(parameter *paraInfo, char *inputSeq, int tsdLen, char strand);*/
polyaStruct *findPolyA(parameter *paraInfo, char *inputSeq, int tsd3Len, int tsd5Len, char strand);/*new*/

void freeFloatMatrix(double **matrix);

void freeIntMatrix(int **matrix);

double max(double m[], int len);

int argmax(double m[], int len);

double scorePairs(char a, char b);

double scoreTsdPairs(char *tsd1, char *tsd2);

int encodeInt(char ch);

CBed6 *parseFillLine(char *line, char *chrom);

CBed6 *parseGapLine(char *line, char *chrom);

char *parseNetLine(char *line, int *chromLen);

void printHeader(parameter *paraInfo, FILE *outfp);

void freePolyaItem(polyaStruct *polya);

#endif /* End retroSeeker_HEAD_H */
