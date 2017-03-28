SRC=$(wildcard *.cxx)
OBJ=$(SRC:.cxx=.ex)
CC=g++ -std=c++11
LIB=-I ${CLASS_include} -L ${CLASS_lib} -L SSENAR -lCLASSpkg -lSsenar `root-config --cflags` `root-config --libs` -fopenmp -lgomp -Wunused-result -lTMVA

all: $(OBJ)

%.ex: %.cxx
	${CC} -o $@  $< ${LIB}

clean:
	rm *.ex
                                                                    

