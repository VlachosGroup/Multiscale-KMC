CC = mpiCC                  #compiler name
CFLAGS = -std=c++0x         #extra compiler flags
OUT_NAME = MSA-KMC-MPI.x           #executable name
OBJ = Traj_stats.o run_sim_par.o STS_sim.o TTS_sim.o in_file_reader.o
SRC = ../src/

all: $(OBJ)
	$(CC) -o $(OUT_NAME) $(CFLAGS) $(OBJ)

%.o: $(SRC)%.cpp
	$(CC) -c $(CFLAGS) $< $(SRC)myHeader.h
    
clean:
	rm $(OUT_NAME) *.o