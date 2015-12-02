#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#define LINE_MAX 300


int classify(double **meanvect, int* qval, int *class, int numrecord, int lengthread,int numclusters){
double *cluster;
  double value=10000000,tmpvalue=0;
  int ctr,ctr1;
  int classprev=class[numrecord];
  for (ctr=0;ctr<numclusters;ctr++){
    tmpvalue=0;
    for(ctr1=0;ctr1<lengthread;ctr1++){
      tmpvalue = tmpvalue + pow((qval[ctr1]-meanvect[ctr][ctr1]),2);
    }
    if(value>tmpvalue){
      class[numrecord]=ctr+1;
      value = tmpvalue;
    }
  }
  if(class[numrecord]==classprev){
    return 0;
  }else{
    return 1;
  }
}

void printarr(int *arr, int length){
  int ctr;
  for (ctr=0;ctr<length;ctr++){
    printf("%d ",arr[ctr]);
  }
  printf("\n");
}

void updatemeans(char* filename, double **meanvect, int *class, int numreads, int lengthread, int numclusters, int numlinestoread){
    int ctr,ctr1,ctr2;
    char lineqval[lengthread];
    int qval[lengthread];
    int haschanged=0;
    FILE *ffastq;
    int clustnum[numclusters];//this stores the number of reads of the specific class which is read
    double tmp;

    //initialize clustnum
    for(ctr=0;ctr<numclusters;ctr++){
      clustnum[ctr] = 0;
    }
    
    ffastq = fopen(filename,"r");
    for (ctr=0; ctr<numreads;ctr++){
      for (ctr2=0;ctr2<numlinestoread;ctr2++){
        fgets(lineqval, LINE_MAX, ffastq);
      }
      clustnum[class[ctr]-1] = clustnum[class[ctr]-1]+1;
      tmp = clustnum[class[ctr]-1];
      //read in quality values
      for(ctr1=0;ctr1<lengthread;ctr1++){
        qval[ctr1] = lineqval[ctr1];
        meanvect[class[ctr]-1][ctr1] = ((tmp-1)*meanvect[class[ctr]-1][ctr1]+qval[ctr1])/tmp;
      }
    }
    fclose(ffastq);
    //printarr(clustnum,numclusters);
}

int updatecluster(char* filename, double **meanvect, int *class, int numreads,int lengthread, int numclusters, int numlinestoread){
    int ctr,ctr1,ctr2;
    char lineqval[lengthread];
    int qval[lengthread];
    int haschanged=0;
    FILE *ffastq;
   
    ffastq = fopen(filename,"r");
    for (ctr=0; ctr<numreads;ctr++){
      for (ctr2=0;ctr2<numlinestoread;ctr2++){
        fgets(lineqval, LINE_MAX, ffastq);
      }
      for(ctr1=0;ctr1<lengthread;ctr1++){
        qval[ctr1] = lineqval[ctr1];
      }
      haschanged = haschanged+classify(meanvect, qval,class,ctr,lengthread,numclusters);
    }
    fclose(ffastq);
    return haschanged;
}


int main(int argc, char *argv[])
{
  unsigned int minqval ;
  unsigned int maxqval ;
  unsigned int numclusters;
  unsigned int lengthread;
  unsigned int numreads;
  double **meanvect;
  int i,ctr,ctr1;
  double step;
  int *class;
  int haschanged=1;
  int MAX_ITER = 100;
  int numlinestoread = 1;

  if (argc != 8){
    printf("Error: You are missing input arguments \n");
    printf("1: File with Qvals \n");
    printf("2: Number of clusters \n");
    printf("3: Min QV \n");
    printf("4: Max QV \n");
    printf("5: Length of the reads \n");
    printf("6: Number of reads \n");
    printf("7: File with noQvals \n");
    exit(1);
  }

  sscanf(argv[2],"%d ", &numclusters);
  sscanf(argv[3],"%d ", &minqval);
  sscanf(argv[4],"%d ", &maxqval);
  sscanf(argv[5],"%d",&lengthread);
  sscanf(argv[6],"%d",&numreads);

  meanvect = (double **) calloc( (numclusters),sizeof(double*));
  for (i = 0; i<numclusters; i++){
    meanvect[i] = (double*) calloc( (lengthread), sizeof(double));
  }
  class = (int*) calloc(numreads,sizeof(int));

  step = (maxqval-minqval)/(numclusters-1.0);
  for (i=0; i< numclusters; i++){
    for (ctr= 0 ; ctr< lengthread; ctr++){
      meanvect[i][ctr] = step*i+minqval;
    }
  }
  ctr = 0;
  printf("Number of reads in each cluster is: \n");
  while(haschanged && ctr < MAX_ITER){
    ctr = ctr+1;
    haschanged = updatecluster(argv[1],meanvect,class,numreads,lengthread,numclusters,numlinestoread);
    updatemeans(argv[1],meanvect,class,numreads,lengthread,numclusters,numlinestoread);
  }
  
  // Create as many files as clusters to write the reads on them 
  FILE **f = malloc(sizeof(FILE*)*(numclusters));
  FILE **f1 = malloc(sizeof(FILE*)*(numclusters));
  for (ctr = 1; ctr<=numclusters; ctr++){
    char* tmp0 = malloc( strlen(argv[1]+2) );
    char* tmp1;
    tmp1 = (char*)calloc(255,sizeof(char));
    sprintf(tmp1, "%d", ctr);
    strcpy(tmp0,argv[1]);
    strcat(tmp0,"_");
    strcat(tmp0,tmp1);
    f[ctr-1] = fopen(tmp0, "w");

    char* tmp00 = malloc( strlen(argv[7]+2) );
    char* tmp11;
    tmp11 = (char*)calloc(255,sizeof(char));
    sprintf(tmp11, "%d", ctr);
    strcpy(tmp00,argv[7]);
    strcat(tmp00,"_");
    strcat(tmp00,tmp11);
    f1[ctr-1] = fopen(tmp00, "w");

  }

  // Write qvs to their corresponding cluster
  FILE *fastq;
  FILE *fnoQval;
  char line[lengthread];
  fastq = fopen(argv[1],"r");
  fnoQval = fopen(argv[7],"r");
  for (ctr = 0; ctr<numreads; ctr++){
    ctr1 = class[ctr];
    for (i = 0; i<numlinestoread; i++){
      fgets(line,LINE_MAX,fastq);
      fprintf(f[ctr1-1],"%s",line);
    }
      fgets(line,LINE_MAX,fnoQval);
      fprintf(f1[ctr1-1],"%s",line);
      fgets(line,LINE_MAX,fnoQval);
      fprintf(f1[ctr1-1],"%s",line);
      fgets(line,LINE_MAX,fnoQval);
      fprintf(f1[ctr1-1],"%s",line);
  }
  free(f);
  free(f1);
  fclose(fastq);
  fclose(fnoQval);
}
