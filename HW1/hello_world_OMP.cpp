/*
openMP Example - Hello World - C/C++ Version
FILE: hello_world_OMP.cpp

Compilation on dmc.asc.edu

   to set number of threads through environment variable in bash shell
   (have to do this before you run the program)
   export OMP_NUM_THREADS=8 

   to compile the program

      g++ -o hello_world_OMP -fopenmp hello_world_OMP.cpp

   to run -- set number of threads as explained above and then type

     ./hello_world_OMP  

*/

using namespace std;
#include <iostream>
#include <omp.h>

int main (int argc, char *argv[]) {

int nthreads, tid;

   // Fork a team of threads giving them their own copies of variables
   #pragma omp parallel private(nthreads, tid)
   {
      // Obtain thread number
      tid = omp_get_thread_num();

      // print out Hello World Message -- do this as a noninterruptable op
      #pragma omp critical
      {
         cout << "Hello World from thread = " << tid << endl;
      }
      // get number of threads
      nthreads = omp_get_num_threads();

      // print out No of Threads Message -- do this as a noninterruptable op
      #pragma omp critical
      {
         // Only master thread does this
         if (tid == 0) {
            cout << "Number of Threads = "  << nthreads << endl;
         }
      }
   }  // All threads join master thread and disband 

}

