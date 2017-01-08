%{
#include <stdio.h>
#include <errno.h>
#include <stdbool.h>
#include <string.h>
#include "mojyh_color.h"

char data[10000000];
int nblank_lines, nlines, nword, nchar, ncols;
bool is_first_ech_line = true;

char* remove_string (int length, char text[length]) {
  int c;
  char * r = malloc(length * sizeof(char));
  for (c = 0; c < length - 1; c++){
    r[c] = text[c];
  }
  return r;
}

void initialize_data () {
  int ech = 0;
  char str[Word];
  for (ech = 0; ech < ncols - 1; ech++) {
    strcat (data, "column"); sprintf(str, "%d",ech); strcat(data, str); strcat (data, " "); 
  }
  strcat (data, "column"); sprintf(str, "%d",ech); strcat(data, str); strcat (data, "\n"); 
}

void check_ncols (int length, char text[length]) {
  int i = 0, w = 0; 
  for (i=0; i<length; i++){
    if (text[i] != ' ' && text[i] != '\t'){
      w++;
      while (text[i] != ' ' && i!= length - 1 && text[i] != '\t')
        {
          i++;
        }
    }
  }
  if (is_first_ech_line) {
    is_first_ech_line = false;
    ncols = w;
    initialize_data ();
  }
  else {if (ncols != w) {printf("Bad number of samples : %s\n", text); abort();}}
}

%}
nl [\n\r]
Emin Emin.+[\n\r]
dE dE.+[\n\r]
int [\+-]?[0-9]+
float {int}\.?{int}
num {float}[eE]?{int}?
num_blank {num}[ \t]+
num_windows {num_blank}+{num}[\r]$
num_list {num_blank}+{num}$

%%
{Emin}|{dE} { ; }
{num_windows} { nlines++; char* tmp = remove_string (yyleng, yytext); check_ncols (yyleng - 1, tmp); strcat(data, tmp); strcat(data, "\n"); yymore(); }
{num_list} { nlines++; check_ncols (yyleng, yytext); strcat(data, yytext); strcat(data, "\n"); yymore(); }
^{nl} { nblank_lines++; nchar++; }
[^ \t\n\r]+ { nword++; nchar+= yyleng; }
.|{nl} { ; }
%%



int main (int argc, char* argv[]) {
// open a file handle to a particular file:

  char* data_file = argv[1];
  if (argc < 3) {
    printf("Please indicate either spectro or quanty as second argument.\n");
    return -1;
  }

  char* filetype = argv[2];
  bool is_quanty = false;
  if (strcmp(filetype,"quanty") == 0) {
    is_quanty = true;
  }
  else if (strcmp(filetype,"spectro") == 0) {
    is_quanty = false;
  }
  else {
    printf("Please indicate either spectro or quanty as second argument.\n");
    return -1;
  }

	FILE *myfile = fopen(data_file, "r");
	// make sure it's valid:
	if (!myfile) {
		printf("I can't open the file!\n");
		return -1;
	}
	yyin = myfile;
	
	// lex through the input:
	yylex();
  fclose(myfile);
  printf("%d\t%d\t%d\n", nblank_lines, nlines, ncols - 1);
  //printf("Data of interest:\n %s\n", data);
  if (is_quanty) {ncols = 3;} 
  return main_mojyh(nlines, ncols - 1, data, is_quanty);

}
