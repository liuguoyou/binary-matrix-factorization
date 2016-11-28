#********************************************************************************
# Authors     : Mehmet Koyuturk & Ananth Grama                                  *
# Institution : Purdue University, Department of Computer Sciences              *
# Last update : 04/01/2004                                                      *
# Copyright   : See README file for copyright information.                      *
#********************************************************************************/

CC      = gcc
CFLAGS  = -O3 -c
LIBS    = -lm

INCLUDE = include
SRC     = src
BIN     = bin

HDRS    = $(INCLUDE)/bnd.h $(INCLUDE)/binvector.h $(INCLUDE)/vector.h $(INCLUDE)/row.h $(INCLUDE)/matrix.h $(INCLUDE)/initialize.h $(INCLUDE)/bndsolve.h $(INCLUDE)/io.h $(INCLUDE)/system.h

OBJS    = $(BIN)/bnd.o $(BIN)/bndsolve.o $(BIN)/initialize.o $(BIN)/io.o $(BIN)/matrix.o $(BIN)/binvector.o $(BIN)/vector.o $(BIN)/system.o

EVAL 	= $(BIN)/matrix.o  $(BIN)/system.o	$(BIN)/evaluate.o	$(BIN)/binvector.o

all:	bnd

bnd: 	$(OBJS) Makefile
	$(CC) $(OBJS) $(LIBS) -o bnd

evaluate:	$(EVAL) Makefile
		$(CC) $(EVAL) $(LIBS) -o evaluate

$(BIN)/system.o : $(HDRS) $(SRC)/system.c
	$(CC) $(CFLAGS) $(SRC)/system.c -o $(BIN)/system.o  

$(BIN)/bnd.o : $(HDRS) $(SRC)/bnd.c
	$(CC) $(CFLAGS) $(SRC)/bnd.c -o $(BIN)/bnd.o  

$(BIN)/bndsolve.o : $(HDRS) $(SRC)/bndsolve.c
	$(CC) $(CFLAGS) $(SRC)/bndsolve.c -o $(BIN)/bndsolve.o  

$(BIN)/initialize.o : $(HDRS) $(SRC)/initialize.c
	$(CC) $(CFLAGS) $(SRC)/initialize.c -o $(BIN)/initialize.o  

$(BIN)/io.o : $(HDRS) $(SRC)/io.c
	$(CC) $(CFLAGS) $(SRC)/io.c -o $(BIN)/io.o  

$(BIN)/matrix.o : $(HDRS) $(SRC)/matrix.c
	$(CC) $(CFLAGS) $(SRC)/matrix.c -o $(BIN)/matrix.o  

$(BIN)/vector.o : $(HDRS) $(SRC)/vector.c
	$(CC) $(CFLAGS) $(SRC)/vector.c -o $(BIN)/vector.o  

$(BIN)/binvector.o : $(HDRS) $(SRC)/binvector.c
	$(CC) $(CFLAGS) $(SRC)/binvector.c -o $(BIN)/binvector.o 

$(BIN)/evaluate.o : $(HDRS) $(SRC)/evaluate.c
	$(CC) $(CFLAGS) $(SRC)/evaluate.c -o $(BIN)/evaluate.o 

clean:
	rm -f bin/*.o src/*~ include/*~ *~ core
