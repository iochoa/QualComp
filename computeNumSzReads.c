#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <fcntl.h>

#define LINE_MAX 300


// Arguments:
// 1: File containing the QVs
// 2: Output file

int main(int argc, char *argv[])
{
  FILE  *f;
  char line[LINE_MAX];
  int szRead = 0;
  int numReads = 0;
  int minQV = 1000;
  int maxQV = 0;
  int value, ctr;

  f = fopen(argv[1],"r");

  while (fgets(line,LINE_MAX, f)!=NULL){
      numReads++; 
      szRead = strlen(line)-1;
      for (ctr=0;ctr<szRead;ctr++){
        value = line[ctr];
        if (value < minQV){
          minQV = value;
        }
        if (value > maxQV){
          maxQV = value;
        }
      }
  }
  fclose(f);
 
  f = fopen(argv[2],"w");

  fprintf(f,"%d\n",szRead);
  fprintf(f,"%d\n",numReads);
  fprintf(f,"%d\n",minQV);
  fprintf(f,"%d\n",maxQV);
 
  fclose(f);
  return 0;
}
    


