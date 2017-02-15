CC = g++                  #compiler name
CFLAGS = -std=c++0x         #extra compiler flags
OUT_NAME = MSA-KMC.x           #executable name
SOURCE = ../Cpp_files/Traj_stats.cpp ../Cpp_files/run_sim_ser.cpp ../Cpp_files/STS_sim.cpp ../Cpp_files/in_file_reader.cpp ../Cpp_files/myHeader.h

all: ../Cpp_files/*.cpp ../Cpp_files/*.h
	$(CC) -o $(OUT_NAME) $(CFLAGS) $(SOURCE)
    
clean:
	rm $(OUT_NAME) *.o