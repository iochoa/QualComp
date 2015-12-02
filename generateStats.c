#include <stdio.h>
#include <stdlib.h>
#include <string.h>

// Arguments :
// 1: file name with ONLY quality values
// 2: Length of the reads
// 3: Number of reads
// 4: Output file

void updateStatistics(double *read, double *mean, double **covariance, int numReads, int szReads){
  int ctr,ctr1 ;
  for (ctr = 0; ctr < szReads; ctr++){
    mean[ctr] = (numReads * mean[ctr] + read[ctr])/(numReads+1);
      for (ctr1 = 0; ctr1 < szReads; ctr1++){
        covariance[ctr][ctr1] = (numReads * covariance[ctr][ctr1] + read[ctr]*read[ctr1])/(numReads+1);
      }
  }
}
int main(int argc, char *argv[])
{
  if(argc<5){
    printf("Error!! \nArguments :\n 1: file name with ONLY quality values\n 2: Length of the reads \n 3: Number of reads \n 4: Output file\n");
    return 1;
  }
  FILE *f;
  FILE *f1;
  int szReads;
  double **covariance;
  double *mean;
  double *tmpread;
  int ctr,ctr1;
  int numTotalReads;

  szReads = atoi(argv[2]);
  numTotalReads = atoi(argv[3]);
  //szReads = (int)argv[2];  
  f = fopen(argv[1],"r");
  f1 = fopen(argv[4],"w");

  
  covariance = (double **) calloc(szReads, sizeof(double*));

  for (ctr = 0; ctr < szReads; ctr++){
    covariance[ctr] = (double *) calloc(szReads, sizeof(double));
  }

  mean = (double *) calloc(szReads, sizeof(double));
  tmpread = (double *) calloc(szReads, sizeof(double));

  for (ctr = 0; ctr < numTotalReads; ctr++){
    for (ctr1 = 0; ctr1 < szReads; ctr1++){
      tmpread[ctr1] = (int) fgetc(f);
    }
    fgetc(f);
    updateStatistics(tmpread,mean,covariance,ctr, szReads);
  }
  for (ctr = 0; ctr < szReads; ctr++){
    for (ctr1 = 0; ctr1 < szReads; ctr1++){
      covariance[ctr][ctr1] = covariance[ctr][ctr1]-mean[ctr]*mean[ctr1];
    }
  }
  for (ctr = 0; ctr < szReads; ctr++){
    for (ctr1 = 0; ctr1 < szReads; ctr1++){
      fprintf(f1,"%f ",covariance[ctr][ctr1]);
    }
      fprintf(f1,"\n");
  }
  for (ctr = 0; ctr < szReads; ctr++){
    fprintf(f1,"%f ",mean[ctr]);
  }
  fprintf(f1,"\n");
  
  fclose(f1);
  fclose(f);
  return 0;
}
