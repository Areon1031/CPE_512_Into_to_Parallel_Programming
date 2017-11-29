#include "md5.h"
#include "sha1.h"
#include "sha256.h"
#include <iostream>
#include <vector>
using namespace std;

int main(int argc, char* argv[])
{
  cout << "MD5    of " << argv[1] << " is: " << md5(argv[1]) << endl;
  cout << "SHA1   of " << argv[1] << " is: " << sha1(argv[1]) << endl;
  cout << "SHA256 of " << argv[1] << " is: " << sha256(argv[1]) << endl;

  std::vector<char> test;


  cin.ignore();
  return 0;
}