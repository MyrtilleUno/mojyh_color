%{
#include <stdio.h>
#include <errno.h>
#include <stdbool.h>
#include "yoshi_color.h"
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
  }
  else {if (ncols != w) {printf("Bad number of samples : %s\n", text); abort();}}
}

%}
nl [\n\r]
int [\+-]?[0-9]+
float {int}\.?{int}
num {float}[eE]?{int}?
num_blank {num}[ \t]+
num_windows {num_blank}+{num}[\r]$
num_list {num_blank}+{num}$

%%
{num_windows} { nlines++;  check_ncols (yyleng - 1, remove_string (yyleng, yytext)); yymore(); }
{num_list} { nlines++; check_ncols (yyleng, yytext); yymore(); }
^{nl} { nblank_lines++; nchar++; }
[^ \t\n\r]+ { nword++; nchar+= yyleng; }
.|{nl} { ; }
%%
int main (int argc, char* argv[]) {
// open a file handle to a particular file:
  char* data_file = argv[1];
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
  return main_yoshi(nlines, ncols-1, data_file);

}
