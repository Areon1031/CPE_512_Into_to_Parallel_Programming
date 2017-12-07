/*
* Serial_Password_Crack.cpp
*
* Program designed to perform the dictionary password attack by brute force.
*
* Note:
*   The hashing algorithms used in this study are from zedwood.com
*   Refer to the license file attached.
*
* Compile
*   g++ md5.cpp sha1.cpp sha256.cpp Serial_Password_Crack.cpp -o Serial_Password_Crack
*
* Run:
*   ./Serial_Password_Crack HASH_MODE DICTIONARY ACTUAL_PASSWORD NUM_THREADS
*
* Author: Kyle Ray
* CPE 512 Intro to Parallel Programming
* Project
* December 5, 2017
*/

#include <iostream>
#include <vector>
#include <fstream>
#include <ctime>
#include <cstdlib> // atoi
using namespace std;

// Hashing Function Include
#include "md5.h"
#include "sha1.h"
#include "sha256.h"

// Prototypes
int readDictionary(string filename, std::vector<string>& passwords);

int main(int argc, char* argv[])
{
  // Variables
  double start, finish;

  if (argc < 4)
  {
    std::cout << "Usage: programName hashMode dictionaryFile passwordToFind" << endl;
    exit(EXIT_FAILURE);
  }

  int hash_mode = atoi(argv[1]); // 0 - MD5, 1 - SHA-1, 2 - SHA-256

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
  std::cout << "Reading Dictionary" << endl;
  if (!readDictionary(argv[2], passwords))
  {
    // Something went wrong with the file
    cout << "Error reading the dictionary file" << endl;
    exit(EXIT_FAILURE);
  }

  // Set the total number of passwords to make each process aware
  total_num_passes = passwords.size();

  // Display the number of Passwords
  cout << "Finished reading Dictionary with " << total_num_passes << " passwords!" << endl;

  cout << "Starting the Timer " << endl;
  cout << "Processing...\n\n\n\n";
  int start_sec = clock();
  int stop_sec;

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
      stop_sec = clock();
      std::cout << "Pass found!" << endl
        << "Password is " << passwords[i] << endl
        << "Hash is " << check_pass << endl;
      break;
    }
  }

  // TODO: Make sure that I can get the timing right, should stop the clock when the process reports
  // that it has found the password.
  // This means that I need to take the wall time of the process that finds the password.
  std::cout << "Time to find the password is " << (stop_sec - start_sec)/double(CLOCKS_PER_SEC) << " seconds" << endl;

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