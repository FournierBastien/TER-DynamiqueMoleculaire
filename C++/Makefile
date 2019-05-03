DIR=$(notdir $(PWD))
MPICXX = mpic++
SRC= $(wildcard *.cpp)
OBJ= $(SRC:.cpp=.o)
EXEC= mpi openmp seq

all : $(DIR)_run_seq $(DIR)_run_mpi $(DIR)_run_openmp

sequentiel: sequentiel.cpp
	g++ -std=c++11 -Wall -Wextra $< -o $@

openmp : openmp.cpp
	g++ -std=c++11 -Wall -Wextra -fopenmp $< -o $@

mpi: mpi.cpp 
	$(MPICXX) -o $@ $^

run_mpi: mpi
	@echo "\033[32m### Building MPI ###\033[0m"
	mpirun -np 4 ./mpi

run_openmp: openmp
	@echo "\033[32m### Building OpenMP ###\033[0m"
	./openmp

run_seq: sequentiel
	@echo "\033[32m### Building Sequentiel ###\033[0m"
	./sequentiel

clean:
	@rm -rf *.o

cleanall: clean
	@rm -rf $(EXEC)