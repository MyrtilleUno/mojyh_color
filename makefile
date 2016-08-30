CC=gcc
LEX=flex
#CFLAGS=-Wall
CFLAGS=
LDFLAGS=
LDLEXFLAGS=-ll -lc
yoshi_color : scanner.c
	$(CC) $(CFLAGS) -o yoshi_color scanner.c yoshi_color.o fmemopen.o $(LDLEXFLAGS)
scanner.c : yoshi_color.o
	$(LEX) -o scanner.c scanner.lex
yoshi_color.o : fmemopen.o
	$(CC) $(CFLAGS) -c yoshi_color.c -o yoshi_color.o $(LDFLAGS)
fmemopen.o : fmemopen.c
	$(CC) $(CFLAGS) -c fmemopen.c -o fmemopen.o

clean : 
	rm -rf *.o yoshi_color
