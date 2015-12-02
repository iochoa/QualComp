#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <fcntl.h>
#include <math.h>

#define LINE_MAX 300

// Arguments:
// 1: Fastq file
// 2: Fastq file
// 3: szReads
// 4: numReads
 
int main(int argc, char *argv[])
{
  FILE *f;
  FILE *f1;
  char linequal[LINE_MAX];
  char linequal1[LINE_MAX];
  int ctr;
  int ctr1 = 0;
  double mse = 0;
  double tmp;
  int szReads;
  int numReads;
  
  f = fopen(argv[1],"r");
  f1 = fopen(argv[2],"r");
  szReads = atoi(argv[3]);
  numReads = atoi(argv[4]);

  for (ctr1=0;ctr1<numReads;ctr1++){
      fgets(linequal,LINE_MAX, f);
      fgets(linequal,LINE_MAX, f);
      fgets(linequal,LINE_MAX, f);
      fgets(linequal,LINE_MAX, f);
      fgets(linequal1,LINE_MAX,f1);
      fgets(linequal1,LINE_MAX,f1);
      fgets(linequal1,LINE_MAX,f1);
      fgets(linequal1,LINE_MAX,f1);
      for (ctr=0;ctr<szReads;ctr++){
        tmp = pow(((double)linequal[ctr] - (double)linequal1[ctr]),2);
        mse = mse + tmp; 
      }
  }
  fclose(f);
  fclose(f1);
  mse = mse/numReads/szReads;
  printf("MSE is: %f \n", mse);
  return 0;
}
