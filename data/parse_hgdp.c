#include <stdlib.h>
#include <stdio.h>
#include <string.h>

int main(int argc, char* argv[]) {
  FILE* fin = NULL;
  FILE* fout = NULL;
  char* line = NULL;
  size_t n = 0;
  ssize_t r;
  fin = fopen(argv[1],"r");
  if(!fin) return 1;

  fout = fopen(argv[2],"w");
  if (!fout) return 2;

  // copy header as is
  // it is also the longest line, so it will suffice for latter calls
  r = getline(&line,&n,fin);
  fputs(line,fout);  
  // now parse each line
  size_t nr = 0;
  size_t nr0 = 0;
  size_t ln = 0;
  while (getline(&line,&n,fin) > 0) {
    char *tok;
    nr = 0;
    tok = strtok(line," \n\t\r");
    while (tok) {
      fprintf(fout,"%s|",tok);
      printf("%s|",tok);
      tok = strtok(NULL," \n\t\r");
      if (tok) nr++;
    }
    if (nr0 == 0) {
      nr0 = nr;
    } else {
      if (nr0 != nr) {
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

