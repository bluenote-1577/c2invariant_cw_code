#!/bin/bash

counter=1
while [ $counter -le 10 ]
do
	let "a=counter * 3"
	echo "MPI"
	openmpi/bin/mpirun -np $a mpi_cust
	echo "MPI"
	echo "THR"
	./threads_cust $a
	echo "THR"
	((counter++))
done
