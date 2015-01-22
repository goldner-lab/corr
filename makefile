CC = g++
CFLAGS = -Wall -g ${INCLUDES}
PREFIX ?= /usr

all: favia monaco fudge

favia : favia.o pakdot.o Getopt.o

# dependencies:
favia.o pakdot.o : pakdot.h

monaco : monaco.o  Getopt.o

monaco.o : pakdot.h

fudge: fudge.o

install : favia
	cp $+ ${PREFIX}/bin
