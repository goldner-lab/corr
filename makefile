CC = /usr/bin/g++ -Wall -g


all: favia monaco

f = favia.o pakdot.o Getopt.o
favia : $(f)
	$(CC) $(f) -o $@

# dependencies:
favia.o pakdot.o : pakdot.h

m = monaco.o  Getopt.o
monaco : $(m)
	$(CC) $(m) -o $@

monaco.o : pakdot.h
