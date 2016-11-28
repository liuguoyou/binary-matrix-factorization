/********************************************************************************
* Authors     : Mehmet Koyuturk & Ananth Grama                                  *
* Institution : Purdue University, Department of Computer Sciences              *
* Last update : 03/31/2004                                                      *
* Copyright   : See README file for copyright information.                      *
********************************************************************************/

#include "../include/io.h"

ArgList *
processArguments(int argc, char *argv[]){
  unsigned int i;
  ArgList *arguments = (ArgList *) malloc (sizeof(ArgList));

  arguments->init       = D_INIT;
  arguments->epsilon    = D_EPSILON;
  arguments->write      = D_WRITE;
  arguments->mcs        = D_MCS;
  arguments->algorithm  = D_ALGORITHM;
  arguments->poolsize   = POOLSIZE;
  arguments->tmpsize    = TMPSIZE;

  if (argc <= 1){
    printf("Usage:\n");
    printf("%s <filename> <optional arguments>\n",argv[0]);
    printf(" <filename> : Name of the file containing input matrix, required.\n");
    printf("  Optional Arguments (Default values in parentheses):\n");
    printf("   -a <algorithm> : Algorithm for rank-one approximation (%d).\n",arguments->algorithm);
    printf("     (%d) Discrete\n", DISCRETE);
    printf("     (%d) Continuous\n", CONTINUOUS);
    printf("   -i <init> : Initialization strategy(%d).\n", arguments->init);
    printf("     (%d) Allones\n", ALLONES);
    printf("     (%d) Center\n", CENTER);
    printf("     (%d) Maximum\n", INITMAX);
    printf("     (%d) Partition\n", PARTITION);
    printf("     (%d) GGG\n", GGG);
    printf("     (%d) Neighbor\n", NEIGHBOR);
    printf("     (%d) Random Row\n", RANDOMROW);
    printf("     (%d) Random\n", RANDOM);
    printf("   -e <threshold> : Bound on normalized hamming radius (%d)\n", arguments->epsilon);
    printf("   -w : Write output vectors (FALSE). If used, presence and pattern\n");
    printf("        vectors will be written in files <filename>.X.out and <filename>.Y.out\n");
    printf("        respectively.\n");
    printf("   -c <minclustersize> : Minimum size of a cluster(%d)\n",arguments->mcs);
    printf("   -m <memorysize> : Log(size of static memory pool)(%d)\n",arguments->poolsize);
    printf("   -t <tmpmemsize> : Log(size of temporary memory pool)(%d)\n",arguments->tmpsize);
    exit(0);
  }
  
  strcpy(arguments->fname, argv[1]); 

  for (i=2; i<argc; i+=2){
    char next = argv[i][1];
    if (argv[i][0] != '-')
      errexit("Unrecognized input symbol: %s.\nRun without arguments to see a list of arguments\n",argv[i]);

    switch ( next ){
    case 'i' : arguments->init       = atoi(argv[i+1]);
      break;
    case 'e' : arguments->epsilon    = atoi(argv[i+1]);
      break;
    case 'w' : arguments->write      = TRUE;
      i--;
      break;
    case 'c' : arguments->mcs        = atoi(argv[i+1]);
      break;
    case 'a' : arguments->algorithm  = atoi(argv[i+1]);
      break;
    case 'm' : arguments->poolsize   = atoi(argv[i+1]);
      break;
    case 't' : arguments->tmpsize    = atoi(argv[i+1]);
      break;
    default: 
      errexit("Unrecognized input symbol: %c.\nRun without arguments to see a list of arguments.\n",
	      next);
    }
  }
  return arguments;
}

void
writeOutput(char *fname, Matrix *X, Matrix *Y){
  FILE *wfp;
  char tmpfname[60];
  Row *nrow;
    
  // Write Presence & Pattern Matrices
  strcpy(tmpfname, fname);
  strcat(tmpfname, ".X.out");
  writeMatrix(X, tmpfname); 

  strcpy(tmpfname, fname);
  strcat(tmpfname, ".Y.out");
  writeMatrix(Y, tmpfname); 

  // Write weight file
  /*strcpy(tmpfname, fname);
  strcat(tmpfname, ".wgt");
  wfp = fopen(tmpfname, "w");
  if (wfp == NULL){
    mprintf("***WARNING: Cannot write into file: %s\n",tmpfname);
    return;
  }
  for(nrow = X->rlist; nrow != (Row *) NULLROW; nrow=nrow->next)
    fprintf(wfp, "%d ", nrow->nz);
  fprintf(wfp,"0");
  fclose(wfp);*/
}

