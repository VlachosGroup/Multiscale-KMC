CC = mpiCC                  #compiler name
CFLAGS = -std=c++0x         #extra compiler flags
OUT_NAME = oop_test_par.x           #executable name

all: *.cpp *.h
	$(CC) -o $(OUT_NAME) $(CFLAGS) Traj_stats.cpp run_sim_par.cpp STS_sim.cpp in_file_reader.cpp myHeader.h
    
clean:
	rm $(OUT_NAME) *.o