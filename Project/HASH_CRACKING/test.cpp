#include <iostream>
#include <string>
using namespace std;

int main(int argc, char* argv[])
{
  string test = "testthisstring";

  int num_test = static_cast<int>(test.c_str());
}