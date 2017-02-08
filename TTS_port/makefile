CC = g++                  #compiler name
CFLAGS = -std=c++0x         #extra compiler flags
OUT_NAME = oop_test.x           #executable name

all: *.cpp
	$(CC) -o $(OUT_NAME) $(CFLAGS) kmc_reps_mpi.cpp
    
clean:
	rm $(OUT_NAME) *.o