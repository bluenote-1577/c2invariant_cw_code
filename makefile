gen: generation_c2.cpp
	c++ -g -o gen -std=c++11 generation_c2.cpp -lgiac -lgmp; 
c2: server_c2.cpp
	c++ -g -o c2cw.out server_c2.cpp -std=c++0x -fopenmp;

