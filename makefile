CC=gcc
LEX=flex
#CFLAGS=-Wall
CFLAGS=
LDFLAGS=-lm
LDLEXFLAGS=-lm -ll -lc
yoshi_color : scanner.c
	$(CC) $(CFLAGS) -o yoshi_color scanner.c yoshi_color.o $(LDLEXFLAGS)
scanner.c : yoshi_color.o
	$(LEX) -o scanner.c scanner.lex
yoshi_color.o : yoshi_color.c
	$(CC) $(CFLAGS) -c yoshi_color.c -o yoshi_color.o $(LDFLAGS)
clean : 
	rm -rf *.o yoshi_color
