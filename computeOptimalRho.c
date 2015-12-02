#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_linalg.h>
#include <math.h>

// Arguments:
// 1: File where singular vector and v are stored
// 2: Length of the reads
// 3: Rate 
// 4: Output file


int main(int argc, char *argv[])
{
  if(argc<5){
    printf("Error!! \nArguments :\n 1: File where singular vector and v are stored\n 2: Length of the reads\n 3: Rate \n 4: Output file\n");
    return 1;
  }
  FILE *fin;
  FILE *fout;
  int szReads;
  double rate,loglambda;
  gsl_vector *s,*rho;
  int ctr,ctr1,continuevar=1;
  double tmp;

  szReads = atoi(argv[2]);
  rate = atof(argv[3]);

  s = gsl_vector_alloc(szReads);
  rho = gsl_vector_alloc(szReads);

  fin = fopen(argv[1],"r");
  gsl_vector_fscanf(fin, s);
  fclose(fin);

  // Compute optimal rho
  ctr = szReads;
  gsl_vector_set_all(rho,-1.0);
  
  if (rate == 0){
    continuevar = 0;
    for (ctr1=0; ctr1<szReads; ctr1++){
      gsl_vector_set(rho,ctr1,0);
    }
  }

  while(continuevar)
  {
    loglambda = 0;
    for(ctr1=0;ctr1<ctr;ctr1++){
      loglambda = loglambda + log2(gsl_vector_get(s,ctr1)*2/szReads*log(2));
    }
    loglambda = (loglambda - 2*rate)/(ctr);
    for(ctr1=0;ctr1<ctr;ctr1++){
      tmp = 0.5*log2(2*log(2)*gsl_vector_get(s,ctr1)/szReads)-0.5*loglambda;
      gsl_vector_set(rho,ctr1,tmp);
    }
    if(gsl_vector_get(rho,ctr-1)>0){
      continuevar = 0;
    }else{
      ctr = ctr - 1;
      gsl_vector_set(rho,ctr,0);
    }
  }


  fout = fopen(argv[4],"w");
  gsl_vector_fprintf(fout,rho,"%f");
  fclose(fout);
  
  gsl_vector_free(rho);
  gsl_vector_free(s);
  return 0;
}
