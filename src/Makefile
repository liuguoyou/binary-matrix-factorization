SRC=$(wildcard *.cpp)
HDR=$(wildcard *.h)
OBJ=bsvd.o \
	initialize_dictionary.o\
	update_dictionary.o\
	update_dictionary_missing_data.o\
	coefficients_update.o\
	random_number_generation.o\
	config.o\
	binmat.o\
	intmat.o\
	entropy_coding.o\
	util.o\
	pbm.o\
	binimage.o\
	pnm.o

# tests
TESTSRC=$(wildcard *_test.cpp)
TEST=$(patsubst %.cpp,%,$(TESTSRC))
# tools
TOOLSRC=$(wildcard *_tool.cpp)
TOOL=$(patsubst %.cpp,%,$(TOOLSRC))
# experiments
EXPSRC=$(wildcard *_exp.cpp)
EXP=$(patsubst %.cpp,%,$(EXPSRC))

ifdef DEBUG
CFLAGS=-DDEBUG -O -Wall -g -fopenmp -march=native -mtune=native
else
CFLAGS=-O3 -Wall -Wno-unused-variable -Wno-unused-parameter -fopenmp -msse3 -march=native -mtune=native
endif

all: $(OBJ) $(TOOL) $(TEST) $(EXP)

tests: $(TEST)

exps: $(EXP)

tools: $(TOOL)

%.o: %.cpp $(HDR) Makefile
	echo $(HDR)
	g++ $(CFLAGS) -c -o $@ $<

%_test: %_test.cpp $(HDR) $(OBJ) Makefile
	g++  $(CFLAGS)  -o $@ $< $(OBJ) -lgsl -lgslcblas -lm

%_exp: %_exp.cpp $(HDR) $(OBJ) Makefile
	g++  $(CFLAGS)  -o $@ $< $(OBJ) -lgsl -lgslcblas -lm

%_tool: %_tool.cpp $(HDR) $(OBJ) Makefile
	g++  $(CFLAGS)  -o $@ $< $(OBJ) -lgsl -lgslcblas -lm

clean:
	rm -f $(TEST) $(TOOL) $(EXP) $(OBJ)
	rm -f *.o
