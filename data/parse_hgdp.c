#include <stdlib.h>
#include <stdio.h>
#include <string.h>
/**
 * This program reads a genome data file formatted as those that can be downloaded from the Human
 * Genome Project, that is, a space separated plain text file where pairs of nucleotides are
 * represented as two adjacent letters as in "AG" and missing data is represented as "--".
 *
 *
 * The output consists of three PBM binary pseudo-image files (lookup "netpnm" or run "man netpnm"
 * to see what a PBM is).  Let n be the number of individuals (columns) in the input file, and m be
 * the number of genes (rows).
 *  
 * The first, mode.pbm, is a 4xn binary pseudo-image where each row corresponds to one of four
 * symbols: "A" is represented as "0 0 0 1", T as "0 0 1 0", G as "0 1 0 0" and C as "1 0 0 0".  The
 * value at each row corresponds to the symbol that occurs more often in the corresponding row in
 * the input file.
 * 
 * The second, dist.pbm, is an (2m)xn binary pseudo-image. Given the pair at row i and column j in
 * the input file, the corresponding pair of values (2*i,j) and (2*i+1,j) of this files are defined
 * as follows:
 * 0 0 ---- both nucleotides coincide with the mode of row i
 * 0 1 ---- one of the two nucleotides is different from the mode of row i
 * 1 0 ---- both nucleotides differ from the i-th row mode
 * 0 0 ---- This value is assumed when the (i,j) datum is missing
 *
 * The third matrix, mask.pbm, cointains a 1 if the (i,j) value in the input file is known
 * and a 0 if it is missing ("--" in the input file)
 */
int main(int argc, char* argv[]) {
  FILE* fin = NULL;
  FILE* fmode = NULL;
  FILE* fdist = NULL;
  FILE* fmask = NULL;

  char* line = NULL;
  char* tmp;
  size_t len;
  size_t m = 0,n = 0;
  ssize_t r;
  char* tok;
  char* modes;
  unsigned long Af = 0;
  unsigned long Tf = 0;
  unsigned long Gf = 0;
  unsigned long Cf = 0;
  
  
  fin = fopen(argv[1],"r");
  if(!fin) return 1;

  //
  // 1 scan number of rows and columns in file
  //
  r = getline(&line,&len,fin); // skip header
  //
  // 1.1 scan number of columns from header
  tok = strtok(line," \n\t\r");
  while (tok) {
    n++;
    tok = strtok(NULL," \n\t\r");
  }
  if (n == 0) {
    fprintf(stderr,"Empty or invalid file: %s\n",argv[1]);
    exit(5);
  }
  printf("Number of columns: %lu\n",n);
  //
  // 1.2 scan number of rows
  //
  while ((r = getline(&line,&len,fin)) > 0) {
    m++;
  }
  printf("Number of rows: %lu\n",m);
  //
  // compute mode
  // 
  fmode = fopen("mode.pbm","w");
  fprintf(fmode,"P4 1 %lu\n",m); 
  if (!fmode)
    return 2;
  modes = (char*) malloc(m*sizeof(char));

  printf("SECOND PASS.\n");
  rewind(fin);
  r = getline(&line,&len,fin); // skip header
  size_t i = 0;
  int c = 0;
  while ((r = getline(&line,&len,fin)) > 0) {
    size_t ni = 0; 
    tok = strtok(line," \n\t\r"); // skip gene name
    Af = Tf = Gf = Cf = 0;
    tok = strtok(NULL," \n\t\r");
    while (tok) {
      ni++;
      switch (tok[0])  { // first letter
      case 'A':
	Af++;
	break;
      case 'T':
	Tf++;
	break;
      case 'G':
	Gf++;
	break;
      case 'C':
	Cf++;
	break;
      case '-':
	break;
      }
      switch (tok[1])  { // first letter
      case 'A':
	Af++;
	break;
      case 'T':
	Tf++;
	break;
      case 'G':
	Gf++;
	break;
      case 'C':
	Cf++;
	break;
      case '-':
	break;
      }
      //
      // next 
      //
      tok = strtok(NULL," \n\t\r");
    }
    if (ni != n) {
      fprintf(stderr,"Wrong number of colums (%lu should be %lu) at line %lu\n",ni,n,m);
      exit(7);
    }
    //
    // compute and output mode to mode.pbm
    // A,T,G,C
    char mode='A';
    size_t Mf = Af;
    if (Tf > Mf) {
      mode='T';
      Mf = Tf;
    }
    if (Gf > Mf) {
      mode='G';
      Mf = Gf;
    }
    if (Cf > Mf) {
      mode='C';
      Mf = Cf;
    }
    modes[i] = mode;
    //
    // write two nibbles to  mode.pbm
    //
    if (i % 2) {
      c <<= 4;
      switch (mode) {
      case 'A':
	c |= 0x01;
	break;
      case 'T':
	c |= 0x02;
 	break;
      case 'G':
	c |= 0x04;
	break;
      case 'C':
	c |= 0x08;
	break;
      }
      fputc(c,fmode);
      c = 0;
    } else {  
      switch (mode) {
      case 'A':
	c = 0x01;
	break;
      case 'T':
	c = 0x02;
	break;
      case 'G':
	c = 0x04;
	break;
      case 'C':
	c = 0x08;
	break;
      }
    }
    i++;
  } // end second pass
  fclose(fmode);
  if (m == 0) {
    fprintf(stderr,"Empty or invalid file: %s\n",argv[1]);
    exit(6);
  }
  printf("THIRD PASS.\n");

  tmp = (char*)malloc(n*sizeof(char));
  //
  // third pass: generate dist and mask
  //
  printf("Input file has %lu rows and %lu columns.\n",m,n);
  
  fdist = fopen("dist.pbm","w");
  if (!fdist) return 3;	 
  fmask = fopen("mask.pbm","w");
  if (!fmask) return 4;
  fprintf(fdist,"P4 %lu %lu 8\n",2*m,n);
  fprintf(fmask,"P4 %lu %lu\n",2*m,n);
  // now parse each line; write to fdist and fmask
  //
  i = 0;  
  rewind(fin);
  r = getline(&line,&len,fin); // skip header
  while (getline(&line,&len,fin) > 0) {
    tok = strtok(line," \n\t\r"); // skip gene name
    tok = strtok(line," \n\t\r");
    size_t j = 0;
    char cd = 0, cm = 0;
    while (tok) {
      const char a = tok[0];
      const char b = tok[1];
      const char mode = modes[i];
      if (a == '-') {
	tmp[j++]='-';
      } else {
	if (a == mode) {
	  if (b == mode) {
	    tmp[j++] = '0';
	  } else {
	    tmp[j++] = '1';
	  }
	} else { // a != mode
	  if (b == mode) {
	    tmp[j++] = '1';
	  } else {
	    tmp[j++] = '2';
	  }
	}
      }
      tok = strtok(NULL," \n\t\r");
    }
    tmp[j] = 0;
    puts(line);
    // now we've got a string of length n with '-','0','1' or '2'
    // we traverse it twice, and write TWO lines  on mask and dist
    // the first row is the 'least significant row'
    // and can only be '1' if dist='1'
    char mask = 0x80;
    for (j = 0; j < n; j++) {      
      switch (tmp[j]) {
      case '-': // 0
	break;
      case '0': // 0 above, 0 below
	cm |= mask;
	break;  
      case '1': // 
	cd |= mask;
	cm |= mask;
	break;
      case '2':
	cm |= mask;
	break;
      }
      mask >>= 1;
      if (!mask) {
	fputc(cd,fdist);
	fputc(cm,fmask);
	mask = 0x80;
      }
    }
    if (mask) { // did not finish byte
      while (mask) {
	cd <<= 1;
	cm <<= 1;
      }
      fputc(cd,fdist);
      fputc(cm,fmask);      
    }
    // second 'most significant' row, can only be '1' if dist='2'
    for (j = 0; j < n; j++) {      
      switch (tmp[j]) {
      case '-': // 0
	break;
      case '0': // 0 above, 0 below
	cm |= mask;
	break;  
      case '1': // 
	cm |= mask;
	break;
      case '2':
	cd |= mask;
	cm |= mask;
	break;
      }
      mask >>= 1;
      if (!mask) {
	fputc(cd,fdist);
	fputc(cm,fmask);
	mask = 0x80;
      }
    }
    if (mask) { // did not finish byte
      while (mask) {
	cd <<= 1;
	cm <<= 1;
      }
      fputc(cd,fdist);
      fputc(cm,fmask);      
    } 
    i++;
    fputc('\n',fmask);
    fputc('\n',fdist);
  } // third pass loop
  free(modes);
  fclose(fin);
  fclose(fmask);
  fclose(fdist);
  free(line);
  return 0;
}
  
