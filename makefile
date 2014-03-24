CC=gcc
#CFLAGS=-Wall
CFLAGS=
LDFLAGS=-lm
yoshi_color : yoshi_color.o
	$(CC) $(CFLAGS) -o yoshi_color yoshi_color.o $(LDFLAGS)
yoshi_color.o : yoshi_color.c
	$(CC) $(CFLAGS) -c yoshi_color.c -o yoshi_color.o $(LDFLAGS)
clean : 
	rm -rf *.o yoshi_color
