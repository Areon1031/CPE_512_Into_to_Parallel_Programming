/*
* OpenMP_Password_Crack.cpp
*
* Program designed to utilize OpenMP parallel library 
* to parallelize the process of password cracking using the dictionary
* attack.  This will take in a text file of passwords as well as the real
* password we are looking for, this is to simulate when a hacker has already
* dumped hashes from the victim machine and is now trying to crack the hash,
* and perform the appropriate hashing technique to each word in the dictionary.
* If the hacker is lucky the victim's actual password will be in the dictionary
* file chosen and the hacker will be able to crack the hash and compromise the
* victim's system.  Performing this attack serially takes a very long time, so
* this program faciliates a study to see how much of a speed up a hacker could
* gain by using OpenMP.
*
* Note:
*   The hashing algorithms used in this study are from zedwood.com
*   Refer to the license file attached.
*
* Compile
*   g++ md5.cpp sha1.cpp sha256.cpp MPI_Password_Crack.cpp -o MPI_Password_Crack -fopenmp
*
* Run:
*   ./MPI_Password_Crack DICTIONARY ACTUAL_PASSWORD NUM_THREADS
*
* Author: Kyle Ray
* CPE 512 Intro to Parallel Programming
* Project
* December 5, 2017
*/

#include <iostream>
#include <vector>
#include <fstream>
using namespace std;

// Hashing Function Include
#include "md5.h"
#include "sha1.h"
#include "sha256.h"

// OpenMP
#ifdef _OPENMP
  #include <omp.h>
#endif

// Prototypes
int readDictionary(string filename, std::vector<string>& passwords);

int main(int argc, char* argv[])
{
  // Variables
  double start, finish;
  const int hash_mode = 2; // 0 - MD5, 1 - SHA-1, 2 - SHA-256

  if (argc < 2)
  {
    std::cout << "Usage: programName dictionaryFile passwordToFind numberOfThreads" << endl;
    exit(EXIT_FAILURE);
  }

  // Actual Password
  std::string actual_pass;
  switch (hash_mode)
  {
  case 0:
    actual_pass = md5(argv[2]);
    break;
  case 1:
    actual_pass = sha1(argv[2]);
    break;
  case 2:
    actual_pass = sha256(argv[2]);
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
  std::cout << "Reading Dictionary" << endl;
  if (!readDictionary(argv[1], passwords))
  {
    // Something went wrong with the file
    exit(EXIT_FAILURE);
  }

  // Set the total number of passwords to make each process aware
  total_num_passes = passwords.size();
  
  // Start the timer
  start = omp_get_wtime();
  

  // Check each local pass hash and if we find the password then quit
  int num_threads;
  if (argc < 3)
    num_threads = 1;
  else
    num_threads = atoi(argv[3]);

  #pragma omp parallel for num_threads(num_threads) schedule(static, 1) 
  for (int i = 0; i < passwords.size(); i++)
  {
    std::string check_pass;
    switch (hash_mode)
    {
    case 0:
      check_pass = md5(passwords[i]);
      break;
    case 1:
      check_pass = sha1(passwords[i]);
      break;
    case 2:
      check_pass = sha256(passwords[i]);
      break;
    }

    // If we have a match then set the flag
    if (actual_pass == check_pass)
    {
      pass_found = true;
      std::cout << "Pass found by process " << omp_get_thread_num() << endl
        << "Password is " << passwords[i] << endl
        << "Hash is " << check_pass << endl;
      break;
    }
  }

  // TODO: Make sure that I can get the timing right, should stop the clock when the process reports
  // that it has found the password.
  // This means that I need to take the wall time of the process that finds the password.
  finish = omp_get_wtime();
  std::cout << "Time to find the password is " << (finish - start) << " seconds" << endl;

  cin.ignore();
  return 0;
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
