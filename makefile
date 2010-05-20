CC = /usr/bin/g++ -Wall -g


all: favia


tp = favia.o pakdot.o
favia : $(tp)
	$(CC) $(tp) -o $@

# dependencies:
favia.o pakdot.o : pakdot.h
