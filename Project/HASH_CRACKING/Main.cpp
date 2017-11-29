#include <iostream>
#include <vector>
#include <fstream>
using namespace std;

// Hashing Function Include
#include "md5.h"
#include "sha1.h"
#include "sha256.h"

// MPI
#include <mpi.h>

// Constants
const int ROOT_MAX_PASS_COUNT = 5000;

int readDictionary(string filename, std::vector<string>& passwords);
void scatter(std::vector<string>& passwords, std::vector<string>& local_passwords, int num_passes, int root, int rank, int numtasks);

int main(int argc, char* argv[])
{
  // Variables
  int numtasks, rank, num;
  MPI_Status status;

  // Start up the MPI Processes
  MPI_Init(&argc, &argv); // initialize MPI environment
  MPI_Comm_size(MPI_COMM_WORLD, &numtasks); // get total number of MPI processes
  MPI_Comm_rank(MPI_COMM_WORLD, &rank); // get unique task id number

  if (rank == 0)
  {
    if (argc < 2)
    {
      cout << "Usage: mpirun -np numtasks programName dictionaryFile passwordToFind" << endl;
      MPI_Abort(MPI_COMM_WORLD, MPI_ERR_INFO);
    }
  }

  // Actual Password
  std::string actual_pass = argv[2];

  // Password Container
  std::vector<string> passwords;

  // Number of passwords each process will have
  int num_passes = 0;

  // Root needs to read in and parse the dictionary file
  // Should then send the size that each process should allocate to hold their part
  // Scatter the dictionary to other tasks.

  if (rank == 0)
  {
    cout << "Reading Dictionary" << endl;
    //passwords = new string[ROOT_MAX_PASS_COUNT];
    //int pass_count = readDictionary(argv[1], passwords);
    if (!readDictionary(argv[1], passwords))
    {
      // Something went wrong with the file
      MPI_Abort(MPI_COMM_WORLD, MPI_ERR_FILE);
    }

    // Calculate the number of passwords for each MPI process
    //num_passes = std::ceil((double)passwords.size() / numtasks);
    num_passes = passwords.size() / numtasks;
  }

  MPI_Bcast(&num_passes, 1, MPI_INT, 0, MPI_COMM_WORLD);

  // Local Set of Passwords
  std::vector<string> local_passwords(num_passes);

  if (rank == 0)
    cout << "Made it to the scatter" << endl;

  scatter(passwords, local_passwords, num_passes, 0, rank, numtasks);

  if (rank == 1)
    cout << "Made it out of the scatter" << endl;

  //MPI_Scatter(passwords, num_passes, MPI_CHAR, local_passwords, num_passes, MPI_CHAR, 0, MPI_COMM_WORLD);

  if (rank == 0)
  {
    for (int i = 0; i < num_passes; i++)
    {
      cout << "Local Password " << i << " : " << local_passwords[i] << endl;
    }
    cout << "Num Passes: " << num_passes << endl;
  }



  // Terminate MPI Program -- perform necessary MPI housekeeping
  // clear out all buffers, remove handlers, etc.
  MPI_Finalize();
}

int readDictionary(string filename, std::vector<string>& passwords)
{
  fstream dictFile(filename.c_str());
  int count = 0;

  if (dictFile.fail())
  {
    dictFile.close();
    return count;
  }

  while (!dictFile.eof())
  {
    // Read and store each password
    string line;
    getline(dictFile, line);
    passwords.push_back(line);
    ++count;
  }

  dictFile.close();
  return count;
}

void scatter(std::vector<string>& passwords, std::vector<string>& local_passwords, int num_passes, int root, int rank, int numtasks)
{
  MPI_Status status;

  if (rank == root)
  {
    int begin_element = 0;
    cout << "List size " << passwords.size() << endl;

    for (int mpitask = 0; mpitask < numtasks; mpitask++)
    {
      if (mpitask == root)
      {
        for (int i = 0; i < num_passes; i++)
        {
          cout << "Begin element: " << begin_element << endl;
          local_passwords[i] = passwords[begin_element];
          ++begin_element;
        }
      }
      else
      {
        for (int i = 0; i < num_passes; i++)
        {
          cout << "Begin element: " << begin_element << endl;
          MPI_Send(passwords[begin_element].c_str(), passwords[begin_element].length(), 
            MPI_CHAR, mpitask, mpitask, MPI_COMM_WORLD);
          ++begin_element;
        }
      }
    }
  }
  
  MPI_Barrier(MPI_COMM_WORLD);
  if (rank != root)
  {
    for (int i = 0; i < num_passes; i++)
    {
      cout << "In recv proc " << rank << " recv loop " << i << endl;
      // Probe the message
      MPI_Probe(root, rank, MPI_COMM_WORLD, &status);
      int len = 0;
      MPI_Get_count(&status, MPI_CHAR, &len);
      char* buf = new char[len];

      MPI_Recv(buf, len, MPI_CHAR,
        root, rank, MPI_COMM_WORLD, &status);

      string buf_string(buf, len);
      cout << "In recv for proc " << rank << " string is " << buf_string << endl;
      local_passwords.push_back(buf_string);
    }
  }
}
