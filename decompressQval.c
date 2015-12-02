#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <randistrs.h>
#include <mtwist.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_blas.h>

#define ARRAYLEN 20

// Arguments:
// 1: File where binary compressed quality values are stored
// 2: File where statistics are stored
// 3: File where svd information and v is stored
// 4: File where rho is stored
// 5: Number of reads  
// 6: Size of each read
// 7: Output file
// 8: minQV
// 9: maxQV

// ###### Functions ######

// This function reads bits from file
int readFile(FILE *f, int *buffer, int tmpbits, int *markerstart, int *markerend, int *lenrecord){

  int ctr1,ctr,numchars=0,tmpnumber,posctr1;

  // Read from file if necessary
  if (tmpbits > *lenrecord){
    numchars = (tmpbits-*lenrecord + 8 - 1)/8;
    for (ctr=0; ctr<numchars; ctr++){
      tmpnumber = fgetc(f);
      for (ctr1=0; ctr1<8; ctr1++){
        buffer[((*markerend)+7-ctr1)%ARRAYLEN] = tmpnumber%2;
        tmpnumber = tmpnumber/2 ;
      }
    }
  }
  *markerend = ((*markerend) + numchars*8)%ARRAYLEN;
  *lenrecord = *lenrecord + numchars*8 - tmpbits;

  // Read tmpint bits from buffer
  tmpnumber = 0;
  for (ctr=0; ctr < tmpbits; ctr++){
    tmpnumber = tmpnumber + pow(2,ctr)*buffer[((*markerstart)+ctr)%ARRAYLEN];
  }
  *markerstart = (*markerstart+tmpbits)%ARRAYLEN;

  return tmpnumber;
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


// ###### main #####

int main(int argc, char *argv[])
{
  if(argc<5){
    printf("Error!! \nArguments :\n 1: File where binary compressed quality values are stored\n 2: File where statistics are stored\n 3: File where svd information and v is stored \n 4: File where rho is stored \n 5: Number of reads \n 6: Size of each reads \n 7: Output file");
    return 1;
  }
  srand(1);
  FILE *fin;
  FILE *fout;
  int szReads, numReads,tmpint,positionMap = 0;
  double rate,loglambda;
  int minQV, maxQV;
  gsl_matrix *v;
  gsl_vector *tmpread1;
  gsl_vector *s,*rho,*mean,*tmpread,*sqrts,*rhoint;
  gsl_matrix *covariance;
  int ctr,ctr1,ctr2,continuevar=1;
  double tmp;
  double **decisionRegions;
  int markerstartread = 0, markerendread = 0;
  int bufferread[ARRAYLEN] = {0};
  int lentotal=0;  
  int numRegions = 20;
  decisionRegions = (double **) calloc(numRegions,sizeof(double*));
  for (ctr = 0; ctr<numRegions; ctr++){
    decisionRegions[ctr] = (double *) calloc(pow(2,ctr+1), sizeof(double));
  }

  numReads = atof(argv[5]);
  szReads = atoi(argv[6]);
  minQV = atoi(argv[8]);
  maxQV = atoi(argv[9]);

  //Read Mean and Covariance matrix from file
  mean = gsl_vector_alloc(szReads);
  covariance = gsl_matrix_alloc(szReads,szReads);
  fin = fopen(argv[2],"r");
  gsl_matrix_fscanf(fin,covariance);
  gsl_vector_fscanf(fin,mean);
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
  fin = fopen("decision_region.dat", "r");
  for (ctr = 0; ctr < numRegions ; ctr ++){
    for (ctr1=0; ctr1<pow(2,ctr+1); ctr1++){
      fscanf(fin,"%lf", &decisionRegions[ctr][ctr1]);
    }
  }
  fclose(fin);

  // Open binary file to be read
  fin = fopen(argv[1],"r");

  // Open file to write the uncompressed Quality Values
  fout = fopen(argv[7],"w");

  tmpread = gsl_vector_alloc(szReads);
  tmpread1 = gsl_vector_alloc(szReads);

  for(ctr = 0; ctr < numReads; ctr++){
    // Randomly determine how many bits to use at each position of the read
    computeNumBits(rho,rhoint,szReads);
    // Read from the binary file as many bits as needed to reconstruct the read
    for(ctr1 = 0; ctr1< szReads; ctr1++){ 
      tmpint = (int) gsl_vector_get(rhoint,ctr1);
      if (tmpint != 0){
        // Read tmpint bits from file and return positionMap
        positionMap = readFile(fin, bufferread, tmpint, &markerstartread, &markerendread, &lentotal);
        // Compute the quality value at each position
	gsl_vector_set(tmpread, ctr1, decisionRegions[tmpint-1][positionMap]);
      }else{
        gsl_vector_set(tmpread, ctr1, 0);
      }
    }
    // Convert to the correspongind normal dristribution
    gsl_vector_mul(tmpread,sqrts);
    gsl_vector_memcpy(tmpread1,tmpread);
    gsl_blas_dgemv(CblasNoTrans, 1.0, v, tmpread1, 0.0, tmpread);
    gsl_vector_add(tmpread,mean);
    
    // Write to file
    for (ctr1=0; ctr1<szReads; ctr1++){
      if ((int)gsl_vector_get(tmpread,ctr1) < minQV){
        fputc(minQV,fout);
      }else if ( (int)gsl_vector_get(tmpread,ctr1) > maxQV){
        fputc(maxQV,fout);
      }else{
        fputc((int)gsl_vector_get(tmpread,ctr1),fout);
      }
    }
    fputc('\n',fout);
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
  return 0;
}
