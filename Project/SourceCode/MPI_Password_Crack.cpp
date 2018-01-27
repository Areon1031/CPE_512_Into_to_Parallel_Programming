/*
* MPI_Password_Crack.cpp
* 
* Program designed to utilize the message passing interface (MPI)
* to parallelize the process of password cracking using the dictionary 
* attack.  This will take in a text file of passwords as well as the real
* password we are looking for, this is to simulate when a hacker has already
* dumped hashes from the victim machine and is now trying to crack the hash,
* and perform the appropriate hashing technique to each word in the dictionary.
* If the hacker is lucky the victim's actual password will be in the dictionary
* file chosen and the hacker will be able to crack the hash and compromise the 
* victim's system.  Performing this attack serially takes a very long time, so
* this program faciliates a study to see how much of a speed up a hacker could 
* gain by using MPI.
*
* Note:
*   The hashing algorithms used in this study are from zedwood.com
*   Refer to the license file attached.
*
* Compile
*   mpic++ md5.cpp sha1.cpp sha256.cpp MPI_Password_Crack.cpp -o MPI_Password_Crack 
*
* Run:
*   mpirun -np NUM_PROCESSES MPI_Password_Crack DICTIONARY ACTUAL_PASSWORD
*
* Author: Kyle Ray
* CPE 512 Intro to Parallel Programming 
* Project
* December 5, 2017
*/

#include <iostream>
#include <vector>
#include <fstream>
#include <cstdlib> // atoi
using namespace std;

// Hashing Function Include
#include "md5.h"
#include "sha1.h"
#include "sha256.h"

// MPI
#include <mpi.h>

// Prototypes
int readDictionary(string filename, std::vector<string>& passwords);
void scatter(std::vector<string>& passwords, std::vector<string>& local_passwords, int num_passes, int root, int rank, int numtasks);

int main(int argc, char* argv[])
{
  // Variables
  int numtasks, rank, num;
  double tot_finish;
  double start, finish;
  int hash_mode; // 0 - MD5, 1 - SHA-1, 2 - SHA-256
  MPI_Status status;

  // Start up the MPI Processes
  MPI_Init(&argc, &argv); // initialize MPI environment
  MPI_Comm_size(MPI_COMM_WORLD, &numtasks); // get total number of MPI processes
  MPI_Comm_rank(MPI_COMM_WORLD, &rank); // get unique task id number

  if (rank == 0)
  {
    if (argc < 3)
    {
      std::cout << "Usage: mpirun -np numtasks programName hashmode dictionaryFile passwordToFind" << endl;
      MPI_Abort(MPI_COMM_WORLD, MPI_ERR_INFO);
    }
  }

  // Store the Hash Mode
  hash_mode = atoi(argv[1]);

  // Actual Password
  std::string actual_pass;
  switch (hash_mode)
  {
    case 0:
      actual_pass = md5(argv[3]);
      break;
    case 1:
      actual_pass = sha1(argv[3]);
      break;
    case 2:
      actual_pass = sha256(argv[3]);
      break;
  }

  // Password Container
  std::vector<string> passwords;

  // Done Flag
  bool pass_found = false;

  // Total number of passwords, used for sizing
  int total_num_passes = 0;

  // Root needs to read in and parse the dictionary file
  // Should then send the size that each process should allocate to hold their part
  // Scatter the dictionary to other tasks.

  if (rank == 0)
  {
    std::cout << "Reading Dictionary" << endl;
    if (!readDictionary(argv[2], passwords))
    {
      // Something went wrong with the file
      MPI_Abort(MPI_COMM_WORLD, MPI_ERR_FILE);
    }

    // Set the total number of passwords to make each process aware
    total_num_passes = passwords.size();

    // Display Size of Dictionary
    std::cout << "Finished reading Dictionary with " << total_num_passes << " passwords!" << endl;
  }

  MPI_Bcast(&total_num_passes, 1, MPI_INT, 0, MPI_COMM_WORLD);

  // Local Set of Passwords
  std::vector<string> local_passwords;

  // Setup a barrier to make sure that all of the processes enter the region at the same time
  // for timing analysis
  MPI_Barrier(MPI_COMM_WORLD);

  if (rank == 0)
  {
    std::cout << "Starting the timer. " << endl;
    std::cout << "Scattering...\n";
    start = MPI_Wtime();   
  }

  // Scatter the passwords to the processes
  scatter(passwords, local_passwords, total_num_passes , 0, rank, numtasks);

  // Flag for a process to report when done, allows other processes to know when to stop
  bool someone_done;

  int pass_count = 0;
  double loc_finish = 0;

  // Check each local pass hash and if we find the password then quit and let all the other processes know
  // Start the main loop
  // 1.) Hash the current password in the MPI process list
  // 2.) Compare the hashed string with the password hash given by the user
  // 3.) If it is a match, then quit and let the other processes know.
  // 4.) Else, process the next password.
  for (int i = 0; i < local_passwords.size(); i++)
  {
    std::string check_pass;
    switch (hash_mode)
    {
      case 0:
        check_pass = md5(local_passwords[i]);    // MD5 algorithm by zedwood.com
        break;
      case 1:
        check_pass = sha1(local_passwords[i]);   // SHA1 algorithm by zedwood.com
        break;
      case 2:
        check_pass = sha256(local_passwords[i]); // SHA256 algorithm by zedwood.com
        break;
    }

    // If we have a match then set the flag
    if (actual_pass == check_pass)
    {
      pass_found = true;
      loc_finish = MPI_Wtime();
      std::cout << "Pass found by process " << rank << endl
        << "Password is " << local_passwords[i] << endl
        << "Hash is " << check_pass << endl;
    }

    // Let other processes know if the password is found or not
    // TODO: Could probably use non blocking send and receives to do this somehow, just to make it more efficient.
    // Although, all of the loops should be somewhere around the same time complexity, so it might not matter too much.
    // The communication takes too much time because it is executing every loop.
    //MPI_Status status;
    //for (int mpitask = 0; mpitask < numtasks; mpitask++)
    //{
    //  // If not current process, send flag to other processes.
    //  if (mpitask != rank)
    //  {
    //    MPI_Send(&pass_found, 1, MPI_C_BOOL, mpitask, rank, MPI_COMM_WORLD);

    //    // Receive flag from other processes
    //    MPI_Recv(&someone_done, 1, MPI_C_BOOL, mpitask, mpitask, MPI_COMM_WORLD, &status);

    //    // If password is found by a process, then set exit condition
    //    if (someone_done)
    //      pass_found = true;
    //  }
    //}

    //pass_count++;

    // If the password has been found, no more processing necessary
    if (pass_found)
      break;
  }

  if (!pass_found)
    loc_finish = MPI_Wtime();

  MPI_Reduce(&loc_finish, &finish, 1, MPI_DOUBLE, MPI_MIN, 0, MPI_COMM_WORLD);
  std::cout << "CLEAN UP: Process " << rank << " is finished processing " << endl;//<< pass_count << " number of passwords ending on " << local_passwords[pass_count-1] << endl;

  MPI_Barrier(MPI_COMM_WORLD);
  
  // Display the results
  if (rank == 0)
  {
    tot_finish = MPI_Wtime();
    std::cout << "Time to find the password is " << (finish - start) << " seconds" << endl;
    std::cout << "Total execution time: " << (tot_finish - start) << " seconds" << endl;
  }


  // Terminate MPI Program -- perform necessary MPI housekeeping
  // clear out all buffers, remove handlers, etc.
  MPI_Finalize();
}

// Method to read in the dictionary text file and store the passwords
// into a vector.
int readDictionary(string filename, std::vector<string>& passwords)
{
  // Create file stream object and open the file
  fstream dictFile(filename.c_str());

  // Keep a count of the number of passwords, needed for equal distribution
  int count = 0;

  if (dictFile.fail())
  {
    dictFile.close();
    return count;
  }

  // Read and store each password
  while (!dictFile.eof())
  {
    string line;
    getline(dictFile, line);
    passwords.push_back(line);
    ++count;
  }

  dictFile.close();
  return count;
}

// Method to equally scatter the dictionary to each MPI process
// TODO: Find a faster way to send strings.
// These blocking send and receives take a very long time.
void scatter(std::vector<string>& passwords, std::vector<string>& local_passwords, int num_passes, int root, int rank, int numtasks)
{
  MPI_Status status;

  // Root handles the distribution of the dictionary file
  if (rank == root)
  {
    int base = num_passes / numtasks;
    int extra = num_passes % numtasks;

    int begin_element = 0;

    for (int mpitask = 0; mpitask < numtasks; mpitask++)
    {
      int num_passes_by_rank = mpitask < extra ? base + 1 : base;
      if (mpitask == root)
      {
        for (int i = 0; i < num_passes_by_rank; i++)
        {
          local_passwords.push_back(passwords[begin_element]);
          ++begin_element;
        }
      }
      else
      {
        for (int i = 0; i < num_passes_by_rank; i++)
        {
          char* pass = const_cast<char*>(passwords[begin_element].c_str());
          MPI_Send(pass, passwords[begin_element].length(), 
            MPI_CHAR, mpitask, mpitask, MPI_COMM_WORLD);
          ++begin_element;
        }
      }
    }
  }
  
  // Other processes receive their share of the dictionary
  if (rank != root)
  {
    int base = num_passes / numtasks;
    int extra = num_passes % numtasks;
    int num_passes_by_rank = rank < extra ? base + 1 : base;
    for (int i = 0; i < num_passes_by_rank; i++)
    {
      // Probe the message
      MPI_Probe(root, rank, MPI_COMM_WORLD, &status);
      int len = 0;
      MPI_Get_count(&status, MPI_CHAR, &len);
      char* buf = new char[len];

      MPI_Recv(buf, len, MPI_CHAR,
        root, rank, MPI_COMM_WORLD, &status);

      string buf_string(buf, len);
      local_passwords.push_back(buf_string);
    }
  }
}
