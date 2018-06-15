# c2invariant_cw_code

Some code for obtaining the c2 invariant of a graph.

Usage:

For parallel processing using MPI
make;
mpirun -np (num_process) mpi_cust

For parallel processing using threads:
make (threads);
./threads_cust (num_threads)
