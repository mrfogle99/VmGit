# Author: David Holmqvist <daae19@student.bth.se>

CXX=g++
CXXFLAGS=-std=c++17 -g -Wunused -Wall -Wunused -pthread 

ifeq ($(BASIC),1)
	CXXFLAGS += -DMAKE_BASERINO
	OUTPUT = pearson
else
	OUTPUT = pearson_par
endif

all: pearson basic

pearson: vector dataset analysis pearson.cpp
	$(CXX) $(CXXFLAGS) pearson.cpp vector.o dataset.o analysis.o -o $(OUTPUT)

analysis: vector analysis.hpp analysis.cpp
	$(CXX) $(CXXFLAGS) -c analysis.cpp -o analysis.o

dataset: vector dataset.hpp dataset.cpp
	$(CXX) $(CXXFLAGS) -c dataset.cpp -o dataset.o

vector: vector.hpp vector.cpp
	$(CXX) $(CXXFLAGS) -c vector.cpp -o vector.o



analysis_basic: vector analysis.hpp analysis.cpp
	$(CXX) $(CXXFLAGS) -DMAKE_BASERINO -c analysis.cpp -o analysis_basic.o

basic: vector dataset analysis_basic pearson.cpp
	$(CXX) $(CXXFLAGS) -DMAKE_BASERINO pearson.cpp vector.o dataset.o analysis_basic.o -o pearson
	
clean:
	rm -rf pearson *.o *.dSYM 2> /dev/null
