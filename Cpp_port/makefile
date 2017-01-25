CC = mpiCC                  #compiler name
CFLAGS = -std=c++0x         #extra compiler flags
OUT_NAME = MSA-KMC.x           #executable name



all: *.cpp *.h 
	$(CC) -o $(OUT_NAME) $(CFLAGS) kmc_reps_mpi.cpp
    
clean:
	rm $(OUT_NAME) *.o