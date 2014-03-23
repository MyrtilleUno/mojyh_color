Yoshi_color : yoshi_color.o
	gcc -o yoshi_color yoshi_color.o
yoshi_color.o : yoshi_color.c
	gcc -c yoshi_color.c -o yoshi_color.o
clean : 
	rm -rf *.o yoshi_color
