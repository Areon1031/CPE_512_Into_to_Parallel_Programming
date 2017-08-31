/* 
   MPI/OpemMP - Hello World - C++ Version (utilizing C function calling Conventions)
   FILE: hybrid.cpp

   Compilation on dmc.asc.edu

   first set up environment by typing from the command line

      module load openmpi

   to compile the program type

      mpic++ -o hybrid -fopenmp hybrid.cpp 

   to run on eight processors type

      mpirun -np 4 ./hybrid
*/

// Kyle Ray
// CPE_512_Intro_to_Parallel_Programming
// August 31, 2017
// Homework #1

// Hybrid version utilizing MPI and OpenMP
// Create 4 main MPI processes
// Each process will have two executing threads

using namespace std;
#include <iostream>
#include <mpi.h>
#include <omp.h>

int main (int argc, char *argv[]) 
{
   MPI_Status status;
   int nmtsks, rank;

   MPI_Init(&argc,&argv); // Initalize MPI environment
   MPI_Comm_size(MPI_COMM_WORLD,&nmtsks); //get total number of processes
   MPI_Comm_rank(MPI_COMM_WORLD,&rank); // get process identity number

   // Create the parallel team of two threads for this MPI process
   #pragma omp parallel num_threads(2)
   {
      // Allow the output statements to finish
      #pragma omp critical
      {
         cout << "Hello World from MPI Process #" << rank << " Thread # " << omp_get_thread_num() << endl << flush;
      }
   }

   // Only root MPI process does this
   if (rank == 0)
   {
      cout << "Number of MPI Processes = " << nmtsks << endl;
   }

  /* Terminate MPI Program -- clear out all buffers */
  MPI_Finalize();

}
