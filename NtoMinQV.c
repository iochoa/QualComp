#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <fcntl.h>

#define LINE_MAX 300

// Arguments:
// 1: Fastq file
// 2: Output file
// 3: szReads
// 4: minQV

int main(int argc, char *argv[])
{
  FILE  *fin, *fout;

  char line[LINE_MAX];
  char lineread[LINE_MAX];
  char linequal[LINE_MAX];
  int ctr;
  int szRead, minQV;

  fin = fopen(argv[1],"r");
  fout = fopen(argv[2], "w");
  szRead = atoi(argv[3]);
  minQV = atoi(argv[4]);

  while (fgets(line,LINE_MAX, fin)!=NULL){
      fprintf(fout,"%s",line);
      
      fgets(lineread,LINE_MAX,fin);
      fprintf(fout,"%s",lineread);

      fgets(line,LINE_MAX,fin);
      fprintf(fout,"%s",line);

      fgets(linequal,LINE_MAX,fin);
      
      for (ctr=0;ctr<szRead;ctr++){
        if (lineread[ctr] == 'N'){
          linequal[ctr] = minQV;
        }
      }
      fprintf(fout,"%s",linequal);
  }
  fclose(fin);
  fclose(fout);
  return 0;
}

