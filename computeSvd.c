#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_linalg.h>

// Arguments:
// 1: File where covariance matrix and mean are stored
// 2: Length of the reads
// 3: Output file


int main(int argc, char *argv[])
{
  if(argc<4){
    printf("Error!! \nArguments :\n 1: File where covariane matrix is stored\n 2: Length of the reads\n 3: Output file\n");
    return 1;
  }
  FILE *f;
  FILE *f1;
  int szReads;
  gsl_matrix *covariance, *v;
  gsl_vector *s,*work;

  szReads = atoi(argv[2]);
  

  covariance = gsl_matrix_alloc(szReads,szReads);
  v = gsl_matrix_alloc(szReads, szReads);
  s = gsl_vector_alloc(szReads);
  work = gsl_vector_alloc(szReads);

  f = fopen(argv[1],"r");
  gsl_matrix_fscanf(f, covariance);
  fclose(f);
  gsl_linalg_SV_decomp_jacobi (covariance, v, s);
  //gsl_linalg_SV_decomp (covariance, v, s, work);
  f1 = fopen(argv[3],"w");
  
  gsl_vector_fprintf(f1,s,"%f");
  fprintf(f1,"\n");
  gsl_matrix_fprintf(f1, v,"%f");
  
  fclose(f1);
  
  //gsl_matrix_fprintf(stdout,covariance,"%f");
  gsl_matrix_free(covariance);
  gsl_matrix_free(v);
  gsl_vector_free(s);
  gsl_vector_free(work);
  return 0;

}
