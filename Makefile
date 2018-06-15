default: int_cust_c2.cpp
	g++ -g -o cust -fopenmp -std=c++0x int_cust_c2.cpp -I/usr/include/openmpi-x86_64/ 
cust: int_cust_c2.cpp
	g++ -g -o cust -fopenmp -std=c++0x int_cust_c2.cpp ;./cust
