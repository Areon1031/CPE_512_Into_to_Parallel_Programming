#include <iostream>
#include <vector>
#include <mpi.h>
#include <fstream>
using namespace std;

// Hashing Function Include
#include "md5.h"
#include "sha1.h"
#include "sha256.h"

bool readDictionary(string filename);

int main(int argc, char* argv[])
{
  // Variables
  int numtasks, rank, num;
  MPI_Status status;

  // Start up the MPI Processes
  MPI_Init(&argc, &argv); // initialize MPI environment
  MPI_Comm_size(MPI_COMM_WORLD, &numtasks); // get total number of MPI processes
  MPI_Comm_rank(MPI_COMM_WORLD, &rank); // get unique task id number

  // Master needs to read in and parse the dictionary file
  // Should then send the size that each process should allocate to hold their part
  // Scatter the dictionary to other tasks.

  return 0;
}

bool readDictionary(string filename)
{
  std::vector<string> passwords;
  fstream dictFile(filename.c_str());

  if(dictFile.fail())
    return false;

  while (!dictFile.eof())
  {
    // Read and store each password

  }

}
