DIR=$(notdir $(PWD))
1;2802;0cSRC=$(DIR).cpp
OBJ=$(DIR)


all : $(DIR)_noopenmp $(DIR)_openmp


noopenmp : $(DIR)_noopenmp


%_openmp.o:%.cpp
	@echo "\033[34m### Building with OpenMP ###\033[0m"
	g++ -fopenmp -std=c++11 -c $<  -o $@

%_noopenmp.o:%.cpp
	@echo "\033[32m### Building without OpenMP ###\033[0m"
	g++ -std=c++11 -c $<  -o $@


$(DIR)_noopenmp:$(DIR)_noopenmp.o
	g++ -o $@  $< -lpthread


$(DIR)_openmp:$(DIR)_openmp.o
	g++ -fopenmp -o $@  $< -lpthread

clean:
	rm -f *.o *~ \#*\# $(DIR) $(DIR)_noopenmp $(DIR)_openmp
