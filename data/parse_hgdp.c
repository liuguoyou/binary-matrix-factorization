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
  size_t m,n = 0;
  ssize_t r;

  unsigned long Af = 0;
  unsigned long Tf = 0;
  unsigned long Gf = 0;
  unsigned long Cf = 0;
  
  
  fin = fopen(argv[1],"r");
  if(!fin) return 1;

  //
  // 1 scan number of rows and columns in file
  //
  r = getline(&line,&n,fin); // skip header
  //
  // 1.1 scan number of columns from header
  tok = strtok(line," \n\t\r");
  while (tok) {
    n++;
    tok = strtok(NULL," \n\t\r");
  }
  if (n == 0) {
    fprintf(stderr,"Empty or invalid file: %s\n,"argv[1]);
    exit(5);
  }
  //
  // 2.1 scan number of rows and compute mode
  //
  fmode = fopen("mode.pbm","w");
  if (!mode) return 2;
  while ((r = getline(&line,&n,fin)) > 0) {
    m++;
    tok = strtok(line," \n\t\r"); // skip gene name
    Af = Tf = Gf = Cf = 0;
    tok = strtok(line," \n\t\r");
    while (tok) {
      switch (*tok)  { // first letter
      case 'A':
	atgc_f
	break;
      case 'T':
	break;
      case 'G':
	break;
      case 'C':
	break;
      case '-':
	break;
      }
      tok = strtok(NULL," \n\t\r");
      if (tok) nr++;
    }
  }
  if (m == 0) {
    fprintf(stderr,"Empty or invalid file: %s\n,"argv[1]);
    exit(6);
  }
  //
  // go back to the top, compute the mode
  //
  printf("Input file has %d rows and %d columns.\n",m,n);
  rewind(fin);
  
  f = fopen("dist.pbm","w");
  if (!mode) return 3;	 
  fmode = fopen("mask.pbm","w");
  if (!mode) return 4;
  
  // copy header as is
  // it is also the longest line, so it will suffice for latter calls
  r = getline(&line,&n,fin);
  fputs(line,fout);  
  // now parse each line
  size_t n = 0; // number of columns in data file
  size_t nri = 0; // number of columns in row i, for control
  size_t m = 0; // number of rows
  while (getline(&line,&n,fin) > 0) {
    char *tok;
    nr = 0;
    tok = strtok(line," \n\t\r"); // skip gene name
    tok = strtok(line," \n\t\r");
    while (tok) {
      fprintf(fout,"%s|",tok);
      printf("%s|",tok);
      tok = strtok(NULL," \n\t\r");
      if (tok) nr++;
    }
    if (nr1 == 0) {
      nr1 = nr;
      printf("Num
    } else {
      if (__1 != nr) {
        fprintf(stderr,"ERROR: wrong number of columns at line %lu\n",ln);
        fclose(fin);
        fclose(fout);
        free(line);
        return 3;
      }
    }
    ln++;
    fputc('\n',fout);
    printf("\n");
  }
  fclose(fin);
  fclose(fout);
  free(line);
  return 0;
}

