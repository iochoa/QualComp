#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "randistrs.h"
#include "mtwist.h"
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_blas.h>

#define ARRAYLEN 20

// Arguments:
// 1: File where quality value strings are stored
// 2: File where statistics are stored
// 3: File where svd information and v is stored
// 4: File where rho is stored
// 5: Number of reads  
// 6: Size of each read
// 7: Output file


// ###### Functions ######

// This function adds number to array, writes it to file and updates position of markers
void writeFile(FILE *f, int *buffer, int tmpbits, int number, int *markerstart, int *markerend, int willend){
  //add bits to buffer
  int ctr1,ctr,numchars=0,tmpnumber,lenrecord,posctr1;
  for (ctr=1;ctr<=tmpbits;ctr++){
    buffer[((*markerend)+ctr)%ARRAYLEN] = number%2;
    number = number/2 ;
  }
  *markerend = ((*markerend) + tmpbits)%ARRAYLEN;

  //write out bits from markerstart to markerend (both inclusive)
  lenrecord = *markerend-*markerstart+1;
  if(lenrecord<=0){
    lenrecord = *markerend + ARRAYLEN - *markerstart + 1;
  }
  //check if lenrecord is >= 8 bis
  if(lenrecord>=8){
    numchars = lenrecord/8;
  }
  if(willend && 8*numchars<lenrecord){
    numchars = numchars+1;
  }
  
//  printf("numchars is %d \n",numchars);
  for(ctr=0;ctr<numchars;ctr++){
    tmpnumber = 0;
    for(ctr1=0; ctr1<8;ctr1++){
      posctr1 = ctr1 + *markerstart;
      tmpnumber = 2*tmpnumber + buffer[posctr1%ARRAYLEN];
    }
    *markerstart = ((*markerstart)+8)%ARRAYLEN;
    fputc(tmpnumber,f);
  }
}

//This generates ineger bits such that mean is equal to rho
void computeNumBits(gsl_vector *rho,gsl_vector *rhoint, int szReads){  
  int ctr;
  double intpart;
  double frac;
  for(ctr=0;ctr<szReads;ctr++){
    frac = modf(gsl_vector_get(rho,ctr),&intpart);
    gsl_vector_set(rhoint,ctr, intpart+(frac>mt_drand()));
  }
}

//This generates position of value in the decision region 
int computeRegion(int tmpint, double* decisionRegions,double value){
  int region  = 0;
  while (value > decisionRegions[region] && region < (pow(2,tmpint) - 1)){
    region = region + 1;
  }
  return region;
}


// ###### main #####

int main(int argc, char *argv[])
{
  if(argc<5){
    printf("Error!! \nArguments :\n 1: File where quality value strings are stored\n 2: File where statistics are stored\n 3: File where svd information and v is stored \n 4: File where rho is stored \n 5: Number of reads \n 6: Size of each reads \n 7: Output file");
    return 1;
  }
  srand(1);
  FILE *fin;
  FILE *fout;
  int szReads, numReads,tmpint,positionMap = 0;
  double rate,loglambda;
  gsl_matrix *v;
  gsl_vector *tmpread1;
  gsl_vector *s,*rho,*mean,*tmpread,*sqrts,*rhoint;
  gsl_matrix *covariance;
  int ctr,ctr1,ctr2,continuevar=1;
  double tmp;
  double **decisionRegions;
  int markerstart=0,markerend=-1;
  int buffer[ARRAYLEN] = {0};
  
  int numRegions = 20;
  decisionRegions = (double **) calloc(numRegions,sizeof(double*));
  for (ctr = 0; ctr<numRegions; ctr++){
    decisionRegions[ctr] = (double *) calloc(pow(2,ctr+1) - 1, sizeof(double));
  }

  numReads = atof(argv[5]);
  szReads = atoi(argv[6]);
  
  //Read Mean and Covariance matrix from file
  mean = gsl_vector_alloc(szReads);
  covariance = gsl_matrix_alloc(szReads,szReads);
  fin = fopen(argv[2],"r");
  gsl_matrix_fscanf(fin,covariance);
  gsl_vector_fscanf(fin,mean);
  //for (ctr=0;ctr<szReads;ctr++){
  //  printf("mean is %d \n", (int)gsl_vector_get(mean,ctr));
  //}
  fclose(fin);

  // Read info for SVD
  s = gsl_vector_alloc(szReads);
  v = gsl_matrix_alloc(szReads,szReads);
  fin = fopen(argv[3],"r");
  gsl_vector_fscanf(fin,s);
  gsl_matrix_fscanf(fin,v);
  fclose(fin); 

  // Compute sigma (square root of variance)
  sqrts = gsl_vector_alloc(szReads);
  for (ctr = 0 ; ctr < szReads; ctr++){
    tmp = sqrt(gsl_vector_get(s,ctr));
    gsl_vector_set(sqrts,ctr,tmp);
  }

  // Read rho, optimal number of bits per position
  rho = gsl_vector_alloc(szReads);
  rhoint = gsl_vector_alloc(szReads);
  fin = fopen(argv[4],"r");
  gsl_vector_fscanf(fin,rho);
  fclose(fin);

  // Open file containing decision regions
  fin = fopen("decision_value.dat", "r");
  for (ctr = 0; ctr < numRegions ; ctr ++){
    for (ctr1=0; ctr1<pow(2,ctr+1) - 1; ctr1++){
      fscanf(fin,"%lf", &decisionRegions[ctr][ctr1]);
    }
  }
  fclose(fin);

  // Open binary file
  fout = fopen(argv[7],"w");

  // int numbitstotal = 0;
  // Read Quality Values to be compressed
  fin = fopen(argv[1],"r");
  tmpread = gsl_vector_alloc(szReads);
  tmpread1 = gsl_vector_alloc(szReads);
  for(ctr = 0; ctr < numReads; ctr++){
    for(ctr1 = 0; ctr1 < szReads; ctr1++){
      gsl_vector_set(tmpread,ctr1,(int) fgetc(fin));
    }
    fgetc(fin);
    
    // Normalize the quality values, and convert to a N(0,1)
    gsl_vector_sub(tmpread,mean);
    gsl_vector_memcpy(tmpread1,tmpread);
    gsl_blas_dgemv (CblasTrans, 1.0,  v, tmpread1, 0.0, tmpread); 
    gsl_vector_div(tmpread,sqrts);
    
    // Randomly determine how many bits use at each position of the read
    computeNumBits(rho,rhoint,szReads);
    
    for(ctr1 = 0; ctr1< szReads; ctr1++){
      tmpint = (int) gsl_vector_get(rhoint,ctr1);
      // printf("number of bits is %d ", tmpint);
      // numbitstotal = numbitstotal + tmpint;
      if (tmpint != 0){
        // Compute the decision region at each position based on the number of bits
        positionMap = computeRegion(tmpint,decisionRegions[tmpint-1],gsl_vector_get(tmpread,ctr1) );
        // printf("and the region is %d \n", positionMap);
        // Write to file with tmpint bits
        writeFile(fout, buffer, tmpint, positionMap, &markerstart, &markerend,0);       
      }
      // If it is the last element of the last read, write the buffer to file
      if (ctr1 == (szReads-1) && ctr == (numReads -1)){
        writeFile(fout, buffer, 0, 0, &markerstart, &markerend,1);       
      }
    }
  }
  fclose(fout);    
  fclose(fin);

  gsl_matrix_free(v);
  gsl_matrix_free(covariance);
  gsl_vector_free(s);
  gsl_vector_free(mean);
  gsl_vector_free(tmpread);
  gsl_vector_free(tmpread1);
  gsl_vector_free(sqrts);
  gsl_vector_free(rhoint);
  gsl_vector_free(rho);
  // printf("numbitstotal is %d \n", numbitstotal);
  return 0;
}
