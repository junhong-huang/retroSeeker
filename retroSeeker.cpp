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

int retroIdx = 1;

void searchRetroGeneInNetFile(parameter *paraInfo, FILE *genomefp, FILE *faifp, FILE *netfp, FILE *outfp, FILE *bedfp, FILE *outfp2) {
  faidxMap faiHash;
  fprintf(stderr, "#read genome fai file\n");
  readFai(faifp, faiHash);

  printHeader(paraInfo, outfp);
  fprintf(stderr, "#identify Retrogenes\n");
  identifyRetroGenes(paraInfo, genomefp, faiHash, netfp, outfp, bedfp, outfp2);
}

void printHeader(parameter *paraInfo, FILE *outfp) {
  if (paraInfo->verbose)
    fprintf(stderr, "print header information\n");

  fprintf(outfp, "#chrom\tchromStart\tchromEnd\tname\tscore\tstrand\t");
  fprintf(outfp, "tsdScore\ttsdLength\tretroLength\t");
  fprintf(outfp, "retroSeq\ttsd5Seq\ttsdPair\ttsd3Seq\t");
  fprintf(outfp, "polyaStart\tpolyaEnd\tpolyaLength\tpolyaScore\tpolyaSeq\n");
}

void identifyRetroGenes(parameter *paraInfo, FILE *genomefp, faidxMap &faiHash,
                        FILE *netfp, FILE *outfp, FILE *bedfp, FILE *outfp2) {
  char *line         = NULL;
  char *chrom        = NULL;
  int chromLen       = 0;
  long int rowNum    = 0;
  long int gapRowNum = 0;
  char *gapLine      = NULL;
  if (paraInfo->bedExist){
    while (line = getLine(bedfp)) {
      if (feof(bedfp) || line == NULL) {
        if (gapLine != NULL)
          safeFree(gapLine);
        safeFree(line);
        break;
      }    
      CBed6 *fillBed = parseBed6Line(line);
      identifyOneRetroGene(paraInfo, genomefp, faiHash, fillBed, chromLen, outfp);
      freeBed6Item(fillBed);
      if (line != NULL)
        safeFree(line);
      rowNum++;
    }/*while*/
  }
  else{
    while (line = getLine(netfp)) {
      if (feof(netfp) || line == NULL) {
        if (gapLine != NULL)
          safeFree(gapLine);
        safeFree(line);
        break;
      }

      if (startStr(line, "#")) {
        continue;
      }  

      // net line
      if (startStr(line, "net")) {
        if (chrom != NULL) {
          safeFree(chrom);
          chromLen = 0;
        }
        chrom = parseNetLine(line, &chromLen);
      }                                 

      // for self-chain net file
      if (paraInfo->netType && startStr(line + 1, "fill")) { 
        CBed6 *fillBed = parseFillLine(line + 1, chrom);
        if (!(outfp2 == NULL)){
            fprintf(outfp2, "%s\t%d\t%d\tfill_%s_%d_%d\t%.0f\t%c\n", fillBed->chrom,
            fillBed->chromStart, fillBed->chromEnd,
            fillBed->name, fillBed->otherIntInfo_1, fillBed->otherIntInfo_2, fillBed->score, fillBed->strand);
        }
        int insertLen  = fillBed->chromEnd - fillBed->chromStart;
        if (insertLen >= paraInfo->minLen && insertLen <= paraInfo->maxLen)
          identifyOneRetroGene(paraInfo, genomefp, faiHash, fillBed, chromLen, outfp);
        freeBed6Item(fillBed);
      } 

     /*gap line*/
      if (!(outfp2 == NULL)){
          if (startStr(line + 2, "gap")) {
            gapLine   = strClone(line);
            CBed6 *gapBed  = parseGapLine(gapLine, chrom);
            if( ((gapBed->chromEnd)-(gapBed->chromStart))>30 ){          
                fprintf(outfp2, "%s\t%d\t%d\tgap_%s_%d_%d\t%.0f\t%c\n", gapBed->chrom,
                gapBed->chromStart, gapBed->chromEnd,
                gapBed->name, gapBed->otherIntInfo_1, gapBed->otherIntInfo_2, gapBed->score, gapBed->strand);
            }
            freeBed6Item(gapBed);
          } 
      }

      /*maybe level 2*/
      if (startStr(line + 3, "fill")) {
        /*"line + 3"表示指针s向后移动3个位置。line 是一个字符数组或者字符串，通过 "line + 3" 可以获得指向 line 数组中第4个元素的指针。*/
        // CBed6 *gapBed  = parseGapLine(gapLine + 2, chrom);
        CBed6 *fillBed = parseFillLine(line + 3, chrom);
        if (!(outfp2 == NULL)){
            fprintf(outfp2, "%s\t%d\t%d\tfill_%s_%d_%d\t%.0f\t%c\n", fillBed->chrom,
            fillBed->chromStart, fillBed->chromEnd,
            fillBed->name, fillBed->otherIntInfo_1, fillBed->otherIntInfo_2, fillBed->score, fillBed->strand);
        }
        int insertLen  = fillBed->chromEnd - fillBed->chromStart;
        if (insertLen >= paraInfo->minLen && insertLen <= paraInfo->maxLen)
          identifyOneRetroGene(paraInfo, genomefp, faiHash, fillBed, chromLen,outfp);

        freeBed6Item(fillBed);
      }

      /*maybe level 1*/
      if ( startStr(line + 1, "fill")) {
          CBed6 *fillBed = parseFillLine(line + 1, chrom);
          if (!(outfp2 == NULL)){
              fprintf(outfp2, "%s\t%d\t%d\tfill_%s_%d_%d\t%.0f\t%c\n", fillBed->chrom,
              fillBed->chromStart, fillBed->chromEnd,
              fillBed->name, fillBed->otherIntInfo_1, fillBed->otherIntInfo_2, fillBed->score, fillBed->strand);
          }
          freeBed6Item(fillBed);
      }

      // if (gapLine != NULL)
      //   safeFree(gapLine);
      //} // gap line
      if (line != NULL)
        safeFree(line);
      rowNum++;
    }    
  }/*else*/

  if (chrom != NULL)
    safeFree(chrom);
}

int identifyOneRetroGene(parameter *paraInfo, FILE *genomefp, faidxMap &faiHash,
                         CBed6 *fillBed, int chromLen, FILE *outfp) {

/*  int externLen = paraInfo->extend;
  int insertLen = fillBed->chromEnd - fillBed->chromStart;
  if (insertLen < externLen * 2)
    externLen = int(insertLen / 2.0);
  int tsdLen  = externLen * 2;*/

  /*junhong add*/
  int externLen = paraInfo->extend;
  int RNA5endUp = paraInfo->RNA5endUp;
  int RNA5endDown = paraInfo->RNA5endDown;
  int RNA3endUp = paraInfo->RNA3endUp;
  int RNA3endDown = paraInfo->RNA3endDown;
  if (fillBed->strand == '-') {
    RNA5endUp =  paraInfo->RNA3endDown;
    RNA5endDown =  paraInfo->RNA3endUp;
    RNA3endUp =  paraInfo->RNA5endDown;
    RNA3endDown =  paraInfo->RNA5endUp;
  } 

  int insertLen = fillBed->chromEnd - fillBed->chromStart;
  if ( insertLen < RNA5endDown+RNA3endUp )
    RNA5endDown = RNA3endUp = int(insertLen / 2.0);
  int tsdLen  = externLen * 2;
  int tsd5Len  = RNA5endUp+RNA5endDown;
  int tsd3Len  = RNA3endUp+RNA3endDown;

  if (paraInfo->verbose)
    fprintf(stderr, "\n\ninsertLen, extendLen, tsdLen, tsd5Len, tsd3Len:\n%d\n%d\n%d\n%d\n%d\n", insertLen, externLen, tsdLen, tsd5Len, tsd3Len);

  int i       = 0;
  char *chrom = fillBed->chrom;
  string chromStr(chrom);
  if (faiHash.find(chromStr) == faiHash.end()) {
    fprintf(stderr, "can't not find the chromosome %s, skip it.\n", chrom);
    return 0;
  }

  faidx *fai     = faiHash[chromStr];

  /*junhong added*/
  if (paraInfo->verbose)
    fprintf(stderr, "\n\nchromLen before modified: %d\n", chromLen);
  if (chromLen == 0){
    chromLen = fai->len;
    if (paraInfo->verbose)
      fprintf(stderr, "chromLen after modified: %d\n", chromLen);
  }
  /*junhong added*/

  int chromStart = fillBed->chromStart - RNA5endUp;
  if (chromStart < 0)
    chromStart = 0;
  int chromEnd = fillBed->chromEnd + RNA3endDown;
  if (chromEnd > chromLen)
    chromEnd = chromLen;
  char *fillExtendSeq = faidxFetchSeq(genomefp, fai, chromStart, chromEnd, fillBed->strand);/*需要考虑fillBed的RNA方向*/
  int fillExtendSeqLen = strlen(fillExtendSeq);

  /*TSD也需要和bed的fill同向！！！*/
  int allocSpace   = 0;
  allocSpace = tsd5Len > tsd3Len ? tsd5Len : tsd3Len;
  char *tsd5Seq   = (char *)safeMalloc(sizeof(char) * (allocSpace + 1));
  char *tsd3Seq   = (char *)safeMalloc(sizeof(char) * (allocSpace + 1));

 if (fillBed->strand == '+') {
    strncpy(tsd5Seq, fillExtendSeq, tsd5Len);/*fillExtendSeq第一个碱基开始取，取tsdLen的长度，即5`端的2e*/
    strncpy(tsd3Seq, fillExtendSeq + (fillExtendSeqLen - tsd3Len), tsd3Len);/*3`端的2e*/
    tsd5Seq[tsd5Len] = '\0';
    tsd3Seq[tsd3Len] = '\0';
 }
 else{/*junhong add because previous 没考虑fillBed的RNA方向*/
      strncpy(tsd5Seq, fillExtendSeq, tsd3Len);
      strncpy(tsd3Seq, fillExtendSeq + (fillExtendSeqLen - tsd5Len), tsd5Len);
      tsd5Seq[tsd3Len] = '\0';
      tsd3Seq[tsd5Len] = '\0';
  }


  int tsd5SeqLen = strlen(tsd5Seq);
  int tsd3SeqLen = strlen(tsd3Seq);
  if (paraInfo->verbose)
    fprintf(stderr, "\n\ntsd5Seq, tsd5SeqLen, fillExtendSeq, fillExtendSeqLen, tsd3Seq, tsd3SeqLen: \n%s\n%d\n%s\n%d\n%s\n%d\n", tsd5Seq, tsd5SeqLen, fillExtendSeq, fillExtendSeqLen, tsd3Seq, tsd3SeqLen);

  tsdSeeker(paraInfo, fillExtendSeq, tsd5Seq, tsd3Seq, chromStart, chromEnd, fillBed, outfp);
  safeFree(tsd5Seq);
  safeFree(tsd3Seq);

  return 1;
}

int tsdSeeker (parameter *paraInfo, char *fillExtendSeq, char *tsd5Seq, char *tsd3Seq,
              int chromStart, int chromEnd, CBed6 *fillBed, FILE *outfp) {

  /*step1: scan TSD*/
  if (paraInfo->verbose)
    fprintf(stderr, "\n\n######step1: scan TSD by tsdSeeker");
  double **scoreMatrix;
  double tmp[4];
  int i, j, k;
  int tsd5Len_to_seek = strlen(tsd5Seq);
  int tsd3Len_to_seek = strlen(tsd3Seq);
  int fillExtendSeqLen = strlen(fillExtendSeq);

  int M       = tsd5Len_to_seek + 1;
  int N       = tsd3Len_to_seek + 1;

  if (paraInfo->verbose)
    fprintf(stderr, "\n\ntsdSeeker: tsd5Seq, tsd5Len_to_seek, M, fillExtendSeq, fillExtendSeqLen, tsd3Seq, tsd3Len_to_seek, N: \n%s\n%d\n%d\n%s\n%d\n%s\n%d\n%d\n", tsd5Seq, tsd5Len_to_seek, M, fillExtendSeq, fillExtendSeqLen, tsd3Seq, tsd3Len_to_seek, N);

  char *seq1  = tsd5Seq;
  char *seq2  = tsd3Seq;
  char *qSeq;
  char *mSeq;
  char *tSeq;
  int allocSpace   = 0;
  int bestI        = 0;
  int bestJ        = 0;
  int traceI       = 0;
  int traceJ       = 0;
  int bestScore    = 0;
  int matNum       = 0;
  int misNum       = 0;
  int insNum       = 0;
  int delNum       = 0;
  int maxPrintRetroLen  = 500;
  double errorRate = 100;

  /*needleman: global alignment*/

  // allocate memory: M表示行数, N表示列数
  scoreMatrix    = (double **)safeMalloc(sizeof(double *) * M);
  scoreMatrix[0] = (double *)safeMalloc(sizeof(double) * (M * N));

  for (i = 1; i < M; i++)
    scoreMatrix[i] = scoreMatrix[0] + N * i;

  /*动态规划算法：寻找两序列的最佳匹配*/
  // initialize scoreMatrix
  for (i = 0; i < N; i++)
    scoreMatrix[0][i] = 0;
  for (i = 0; i < M; i++)
    scoreMatrix[i][0] = 0;

  // matrix score
  for (i = 1; i < M; i++) {/*fill the matrix row by row*/
    for (j = 1; j < N; j++) {
      /*当下最优解：当下碱基+上一个碱基。 当下的选择：要么gap要么不gap*/
      tmp[0] = scoreMatrix[i - 1][j - 1] + scorePairs(seq1[i - 1], seq2[j - 1]); /*不gap：要么配，要不错配*//* sequence index is i-1 and j-1*/
      tmp[1] = scoreMatrix[i - 1][j] + GAP_OPEN;/*gap*/
      tmp[2] = scoreMatrix[i][j - 1] + GAP_OPEN;/*gap*/
      tmp[3] = 0;
      scoreMatrix[i][j] = max(tmp, 4);
      if (scoreMatrix[i][j] >= bestScore) {
        bestI     = i;
        bestJ     = j;
        bestScore = scoreMatrix[i][j];
      }
    }
  }

  // trace back
  allocSpace = M > N ? M : N;
  qSeq       = (char *)safeMalloc(sizeof(char) * (allocSpace + 1));
  mSeq       = (char *)safeMalloc(sizeof(char) * (allocSpace + 1));
  tSeq       = (char *)safeMalloc(sizeof(char) * (allocSpace + 1));

  i = bestI;
  j = bestJ;
  k = 0;
  for (;;) {/*类似while(1)*/
    if (i == 0 && j == 0)
      break;
    if (scoreMatrix[i][j] <= 0) {
      traceI = i;
      traceJ = j;
      break;
    }
    if (k >= allocSpace) {
      qSeq       = (char *)safeRealloc(qSeq, k + 2);
      mSeq       = (char *)safeRealloc(mSeq, k + 2);
      tSeq       = (char *)safeRealloc(tSeq, k + 2);
      allocSpace = k;
    }
    if (i >= 1 && j >= 1 &&
        scoreMatrix[i][j] == scoreMatrix[i - 1][j - 1] + scorePairs(seq1[i - 1], seq2[j - 1])) {
      qSeq[k] = seq1[i - 1];
      if (seq1[i - 1] == seq2[j - 1]) {
        mSeq[k] = '|';
        matNum++;
      } else {
        mSeq[k] = '.';
        misNum++;
      }
      tSeq[k] = seq2[j - 1];
      i--;
      j--;
    } else if (i >= 1 && j >= 0 &&
               scoreMatrix[i][j] == scoreMatrix[i - 1][j] + GAP_OPEN) {
      qSeq[k] = seq1[i - 1];
      mSeq[k] = '-';
      tSeq[k] = '-';
      i--;
      insNum++;
    } else if (j >= 1 && i >= 0 &&
               scoreMatrix[i][j] == scoreMatrix[i][j - 1] + GAP_OPEN) {
      qSeq[k] = '-';
      mSeq[k] = '-';
      tSeq[k] = seq2[j - 1];
      j--;
      delNum++;
    }
    k++;
  }

  qSeq[k] = '\0';
  mSeq[k] = '\0';
  tSeq[k] = '\0';

  reverseBytes(qSeq, strlen(qSeq));
  reverseBytes(mSeq, strlen(mSeq));
  reverseBytes(tSeq, strlen(tSeq));

  int totalLen = matNum + misNum + insNum + delNum;
  int errorLen = misNum + insNum + delNum;
  int diffMat  = matNum - misNum - insNum - delNum;
  errorRate    = (double)errorLen / (double)totalLen;

  if (paraInfo->verbose)
      fprintf(stderr, "\n\nDynamic programming algorithm, trace back：bestI, traceI, bestJ, traceJ: \n%d\n%d\n%d\n%d\n\n", bestI, traceI, bestJ, traceJ);

  if (matNum >= paraInfo->minMatchLen) {
    /*满足了TSD匹配得分的，fillExtend就变成了retroEleSeq了*/
    /*retroEleSeq: fill+2tsd*/
    int strand        = '+';
    int fillExtendLen      = strlen(fillExtendSeq);
    
    /*新的结束：bestJ*/
    int endPos        = fillExtendLen - (tsd3Len_to_seek - bestJ);
    int endBase       = fillExtendSeq[endPos];
    fillExtendSeq[endPos]  = '\0';/*新endBase后面的字符串全截掉*/

    /*新的开始：traceI*/
    int startPos      = traceI;/*因为用的是整个extend-fill去扫最佳匹配，似乎没有限制TSD5不能落在fill上*/
    char *retroEleSeq = strClone(fillExtendSeq + startPos);/*复制字符串fillExtendSeq从startPos位置开始的子字符串*/
   
    int retroEleLen   = strlen(retroEleSeq);
    int matchLen      = strlen(mSeq);

  if (paraInfo->verbose)
      fprintf(stderr, "\n\nretroEleSeq, retroEleLen, mSeq, matNum, matchLen, insNum, tsdLen: \n%s\n%d\n%s\n%d\n%d\n%d\n%d\n\n", retroEleSeq, retroEleLen, mSeq, matNum, matchLen, insNum, (matchLen-insNum) );





    /*step2 :scan polyA*/
  if (paraInfo->verbose)
      fprintf(stderr, "\n\n######step2 :scan polyA");
  if (retroEleLen >= paraInfo->minLen && retroEleLen <= paraInfo->maxLen) {
      int retroChromStart = chromStart + startPos;
      int retroChromEnd   = retroChromStart + retroEleLen;
      if (fillBed->strand == '-') {
        retroChromEnd   = chromEnd - startPos;
        retroChromStart = retroChromEnd - retroEleLen;
      }

      polyaStruct *pPolyA = NULL; // plus strand
      polyaStruct *mPolyA = NULL; // minus strand

      pPolyA = findPolyA(paraInfo, retroEleSeq, (matchLen - insNum), (matchLen - delNum), '+');
      int plusScore          = bestScore + pPolyA->polyaScore;
      int totalScore         = plusScore;     
      polyaStruct *polyaItem = pPolyA;

      if (paraInfo->netExist){/*junhong add*/
          mPolyA = findPolyA(paraInfo, retroEleSeq, (matchLen - delNum), (matchLen - insNum), '-');
          int minusScore         = bestScore + mPolyA->polyaScore;
          if (minusScore > plusScore) {
            strand = '-';
            if (fillBed->strand == '+') {
              fillBed->strand = '-';
            } else if (fillBed->strand == '-') {
              fillBed->strand = '+';
            }
            totalScore = minusScore;/*refresh*/
            polyaItem  = mPolyA;/*refresh*/
            reverseComp(retroEleSeq);
            reverseComp(qSeq);
            reverseComp(tSeq);
            reverseBytes(mSeq, strlen(mSeq));
            char *tmpSeq = qSeq;
            qSeq         = tSeq;
            tSeq         = tmpSeq;
          }
      }
      
      int tsdStart = 0;
      int tsdScore = bestScore;
      int trimNum  = 0;

      if (polyaItem->end > retroEleLen - matchLen) {
        tsdStart   = getTsdStart(retroEleSeq, tSeq, polyaItem->end);
        tsdScore   = scoreTsdPairs(qSeq + tsdStart, tSeq + tsdStart);
        totalScore = polyaItem->polyaScore + tsdScore;
        trimNum    = get5pTrimNum(qSeq, tsdStart);/*留给TSD 7nt，剩下的尽可能留给PolyA, 此处的作用就是刷新TSD的区域*/
      }
      int tsdLen = strlen(qSeq + tsdStart);


      // fprintf(stderr, "%d\t%d\n", tsdStart, trimNum);
      if (fillBed->strand == '-') {
        retroChromEnd = retroChromEnd - trimNum;
      } else {
        retroChromStart = retroChromStart + trimNum;
      }

      if (totalScore >= paraInfo->minScore) {
        char *newRetroEleSeq = retroEleSeq + trimNum;
        retroEleLen          = strlen(newRetroEleSeq);
        fprintf(outfp, "%s\t%d\t%d\tretroSeeker_%d\t%d\t%c\t", fillBed->chrom,
                retroChromStart, retroChromEnd, retroIdx, totalScore, fillBed->strand);
        fprintf(outfp, "%d\t%d\t%d\t", tsdScore, tsdLen, retroEleLen);
        if (retroEleLen > maxPrintRetroLen) {
          int cutLen             = int(maxPrintRetroLen / 2.0);
          char retroBase         = newRetroEleSeq[cutLen];
          newRetroEleSeq[cutLen] = '\0';
          fprintf(outfp, "%sXXXXXX", newRetroEleSeq);
          newRetroEleSeq[cutLen] = retroBase;
          fprintf(outfp, "%s\t", newRetroEleSeq + (retroEleLen - cutLen));
        } else {
          fprintf(outfp, "%s\t", newRetroEleSeq);
        }
        fprintf(outfp, "%s\t", qSeq + tsdStart);
        fprintf(outfp, "%s\t", mSeq + tsdStart);
        fprintf(outfp, "%s\t", tSeq + tsdStart);

        if (polyaItem->polyaScore > 0) {
          /*polyA_TSD3_dist = retroEleLen -  polyaItem->end - tsd3length;*/
          fprintf(outfp, "%d\t%d\t%d\t%d\t%s", polyaItem->start - trimNum,
                  polyaItem->end - trimNum, polyaItem->polyaLen,
                  polyaItem->polyaScore, polyaItem->polyaSeq);
        } else {
          fprintf(outfp, "0\t0\t0\t0\tnoPolyA");
        }
        fprintf(outfp, "\n");
        retroIdx++;
      }
      freePolyaItem(pPolyA);
      if (paraInfo->netExist)
        freePolyaItem(mPolyA);
    }
    fillExtendSeq[endPos] = endBase;
    safeFree(retroEleSeq);

  } // TSD

  freeFloatMatrix(scoreMatrix);
  safeFree(qSeq);
  safeFree(mSeq);
  safeFree(tSeq);
  return matNum;
}


int getTSDnum(char *seq) {
  int i       = 0;
  int baseNum = 0;
  int seqLen  = strlen(seq);
  for (i = 0; i < seqLen; i++) {
    if (seq[i] != '-' && seq[i] != '.')
      baseNum += 1;
  }
  return baseNum;
}

int get5pTrimNum(char *seq, int tsdOffset) {
  int i       = 0;
  int trimNum = 0;
  int seqLen  = strlen(seq);
  for (i = 0; i < tsdOffset && i < seqLen; i++) {
    if (seq[i] != '-')
      trimNum += 1;
  }
  return trimNum;
}

int getTsdStart(char *retroSeq, char *tSeq, int polyaEnd) {
  int i        = 0;
  int j        = 0;
  int offset   = 0;
  int tsdLen   = strlen(tSeq);
  int retroLen = strlen(retroSeq);
  for (i = tsdLen - 1; i >= 0; i--) {
    j += 1;
    if (tSeq[i] != '-') {
      retroLen -= 1;
    }
    if (retroLen == polyaEnd) {
      break;
    }
  }
  offset = tsdLen - j;
  if (offset < 0)
    offset = 0;
  return offset;
}



/*简而言之，就是任何一段序列的区域均由i index +j index构成, 现在统计任意j开始和任意i结束的A含量*/
/*polyAoffset=15决定了i能运动的范围是TSD前15nt内，而由于maxPolyaLen=5, j也只能在i前面50nt内运动*/
/*特别注意一下，搜索是从后往前搜的，searchStart反而是从polyAend开始的*/
polyaStruct *findPolyA(parameter *paraInfo, char *inputSeq, int tsd3Len, int tsd5Len, char strand) {
  int i           = 0;
  int j           = 0;
  int polyAoffset = 15; // maximum distance between polyA and Tsd: 15
  char *retroSeq  = strClone(inputSeq);
  if (strand == '-')
    reverseComp(retroSeq);
  int retroLen = strlen(retroSeq);
  int polyAend_confineRegion_start = retroLen - tsd3Len - polyAoffset;/*离3`-TSD 15 nt*/
  int polyAend_confineRegion_end   = retroLen - (MIN_TSD_LEN_FOR_GREEDY_POLYA_SEARCH) - 1;/*给TSD3 留出最小的长度7，剩下的范围都可以用来搜索A*/
/*  int polyAend_confineRegion_end   = retroLen - (paraInfo->minMatchLen) - 1;*/
  // int end                = retroLen - tsdLen + MIN_TSD_LEN - 1;
  if (paraInfo->verbose){
      fprintf(stderr, "\n\ntsd3Len, polyAoffset: \n%d\n%d", tsd3Len, polyAoffset);
      fprintf(stderr, "\n\nsearch poly A in candaidate retroSeq (already with tsd5 and tsd3),  polyAend_confineRegion_start, polyAend_confineRegion_end, tsd3Len: \n%s\n%d\n%d\n%d\n\n", retroSeq, polyAend_confineRegion_start, polyAend_confineRegion_end, tsd3Len);    
  }


  int bestScore          = 0;
  int score              = 0;
  int bestStart          = 0;
  int bestEnd            = 0;
  polyaStruct *polyaItem = (polyaStruct *)safeMalloc(sizeof(polyaStruct));
  polyaItem->start       = 0;
  polyaItem->end         = 0;
  polyaItem->polyaScore  = 0;
  polyaItem->polyaLen    = 0;
  polyaItem->polyaSeq    = NULL;

  if (polyAend_confineRegion_start < (tsd3Len + paraInfo->minPolyaLen))/*5*/
      polyAend_confineRegion_start = tsd3Len + paraInfo->minPolyaLen;/*什么意思?*/

  if (paraInfo->verbose){
      fprintf(stderr, "\n\npolyAend_confineRegion_start, polyAend_confineRegion_end: \n%d\n%d\n\n", polyAend_confineRegion_start,polyAend_confineRegion_end);    
  }

  for (i = polyAend_confineRegion_end; i >= polyAend_confineRegion_start; i--) {/*end->start(离tsd3:15nt) 倒回去搜索*/
    score = 0;

    /*redefined the start and end*/
    int searchStart = i;
    int polyAstart_confineRegion = i - (paraInfo->maxPolyaLen);
    int searchEnd   = polyAstart_confineRegion;/*searchStart-50, SearchRegion=50*/
    if (searchEnd < ( tsd5Len + (paraInfo->minGenebodyLen) ) )/*junhong add minGenebodyLen*/
      searchEnd = tsd5Len + paraInfo->minGenebodyLen;

    if (paraInfo->verbose){
        fprintf(stderr, "\n\npolyAend_confineRegion_end, polyAend_confineRegion_start, searchStart, searchEnd: \n%d\n%d\n%d\n%d\n\n", polyAend_confineRegion_end, polyAend_confineRegion_start, searchStart, searchEnd);    
    }

    for (j = searchStart; j >= searchEnd; j--) {/*确定i后，统计以每个i为末端，前面50nt内,以任意一个j开头的丰度*/
      if (retroSeq[j] == 'A' || retroSeq[j] == 'a') {
        score += POLYA_MATCH;
      } else {
        score += POLYA_MISMATCH;
      }
      if (score > bestScore) {
        bestScore = score;
        bestStart = j;
        bestEnd   = i;
      }
    }
  }

  int polyaLen = bestEnd - bestStart + 1;
  if (paraInfo->verbose){
      fprintf(stderr, "\n\npolyaLen, bestEnd(i), bestStart(j): \n%d\n%d\n%d\n\n", polyaLen, bestEnd, bestStart);    
  }

  /*safe the polya result*/
  if (polyaLen >= paraInfo->minPolyaLen) {
    char *polyaSeq = (char *)safeMalloc(sizeof(char) * (polyaLen + 1));
    strncpy(polyaSeq, retroSeq + bestStart, polyaLen);
    polyaSeq[polyaLen]    = '\0';
    polyaItem->polyaSeq   = polyaSeq;
    polyaItem->polyaLen   = polyaLen;
    polyaItem->start      = bestStart + 1;
    polyaItem->end        = bestEnd + 1;
    polyaItem->polyaScore = bestScore;
  }
  safeFree(retroSeq);
  return polyaItem;
}

void freePolyaItem(polyaStruct *polya)
/* free float matrix */
{
  if (polya->polyaSeq != NULL)
    safeFree(polya->polyaSeq);
  safeFree(polya);
}

void freeFloatMatrix(double **matrix)
/* free float matrix */
{
  safeFree(matrix[0]);
  safeFree(matrix);
}

void freeIntMatrix(int **matrix)
/*free int matrix */
{
  safeFree(matrix[0]);
  safeFree(matrix);
}

double max(double m[], int len)
/* return maxinum value */
{
  int i           = 0;
  double maxScore = m[0];

  for (i = 1; i < len; i++) {
    if (m[i] > maxScore) {
      maxScore = m[i];
    }
  }
  return maxScore;
}

int argmax(double m[], int len)
/* return the index for maximum value */
{
  int i           = 0;
  int best_i      = 0;
  double maxScore = m[0];
  for (i = 1; i < len; i++) {
    if (m[i] > maxScore) {
      maxScore = m[i];
      best_i   = i;
    }
  }
  return best_i;
}

int encodeInt(char ch) {
  /* translate character to number */

  ch = toupper(ch);

  if (ch == 'A')
    return 0;
  if (ch == 'C')
    return 1;
  if (ch == 'G')
    return 2;
  if (ch == 'U' || ch == 'T')
    return 3;

  return 4;
}

double scorePairs(char a, char b)
/* score match and mismatch */
{
  double pairMatrix[5][5] = {
  /* A  C  G  T  N*/
      {2,  -3, -3, -3, -3}, // A
      {-3, 2,  -3, -3, -3}, // C
      {-3, -3, 2,  -3, -3}, // G
      {-3, -3, -3, 2,  -3}, // T
      {-3, -3, -3, -3, -3}  // N
  };
  return pairMatrix[(int)(encodeInt(a))][(int)(encodeInt(b))];
}

double scoreTsdPairs(char *tsd1, char *tsd2)
/* score match and mismatch */
{
  int i        = 0;
  int tsdLen   = strlen(tsd1);
  double score = 0;
  for (i = 0; i < tsdLen; i++) {
    score += scorePairs(tsd1[i], tsd2[i]);
  }
  return score;
}

CBed6 *parseFillLine(char *line, char *chrom) {
  CBed6 *bed6   = NULL;
  int fieldNum  = 0;
  char **fields = NULL;
  fields        = splitWhitespace(line, &fieldNum); // with three spaces
  if (fieldNum >= 5) {
    bed6             = (CBed6 *)safeMalloc(sizeof(CBed6));
    bed6->chrom      = strClone(chrom);
    bed6->name     = strClone(fields[3]);
    bed6->chromStart = atoi(fields[1]);
    bed6->score      = atoi(fields[2]);
    bed6->chromEnd   = bed6->chromStart + bed6->score;
    bed6->strand     = fields[4][0];
    bed6->otherIntInfo_1  = atoi(fields[5]);
    bed6->otherIntInfo_2  = atoi(fields[6]);
  }
  freeWords(fields, fieldNum);
  return bed6;
}

CBed6 *parseGapLine(char *line, char *chrom) {
  CBed6 *bed6   = NULL;
  int fieldNum  = 0;
  char **fields = NULL;
  fields        = splitWhitespace(line, &fieldNum); // with two spaces
  if (fieldNum >= 5) {
    bed6             = (CBed6 *)safeMalloc(sizeof(CBed6));
    bed6->chrom      = strClone(chrom);
    bed6->name       = strClone(fields[3]);
    bed6->chromStart = atoi(fields[1]);
    bed6->score      = atoi(fields[2]);
    bed6->chromEnd   = bed6->chromStart + bed6->score;
    bed6->strand     = fields[4][0];
    bed6->otherIntInfo_1  = atoi(fields[5]);
    bed6->otherIntInfo_2  = atoi(fields[6]);
  }
  freeWords(fields, fieldNum);
  return bed6;
}

char *parseNetLine(char *line, int *chromLen) {
  char *chrom   = NULL;
  int fieldNum  = 0;
  char **fields = NULL;
  fields        = splitWhitespace(line, &fieldNum);
  if (fieldNum == 3) {
    if (chrom != NULL)
      safeFree(chrom);
    chrom     = strClone(fields[1]);
    *chromLen = atoi(fields[2]);
  } else {
    fprintf(stderr, "Error Net line %s with %d elements\n", line, fieldNum);
    freeWords(fields, fieldNum);
    safeFree(line);
    exit(1);
  }
  freeWords(fields, fieldNum);
  return chrom;
}


