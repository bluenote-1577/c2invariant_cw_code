default: mpi_cust_c2.cpp
	openmpi/bin/mpic++ -o mpi_cust -fopenmp -std=c++0x mpi_cust_c2.cpp
threads: threads_cust_c2.cpp
	c++ -o threads_cust -fopenmp -std=c++0x threads_cust_c2.cpp
