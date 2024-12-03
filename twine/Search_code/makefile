CX = g++
CXXFLAGS = -O3 -march=native -Wall -Wextra -std=c++11 -fopenmp -Wno-deprecated-copy #-g -fsanitize=address -fno-omit-frame-pointer

IFLAGS = -I $(GUROBI_HOME)/include/ -I /usr/include/x86_64-linux-gnu/
LFLAGS = -L $(GUROBI_HOME)/lib/ -lgurobi_g++5.2 -lgurobi100 #-lgmpxx -lgmp

%.o: %.cpp
	$(CXX) $(IFLAGS) $(LFLAGS) $(CXXFLAGS) -c $< -o $@ 

all : main clean

main :  twine.o SearchKey.o Callback.o
	$(CXX) $(CXXFLAGS) $(IFLAGS) -o twine.run twine.o SearchKey.o Callback.o $(LFLAGS)

clean :
	rm -rf *.o
