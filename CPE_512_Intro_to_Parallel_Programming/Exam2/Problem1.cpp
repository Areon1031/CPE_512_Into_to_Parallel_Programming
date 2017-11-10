#include <iostream>
using namespace std;

int main(int argc, char* argv[])
{
	int i; double sum = 0.0;
#pragma omp parallel for
	for (i = 1; i <= 4; i++)
		sum = sum + i;
	cout << "The sum is " << sum << endl;


	cin.ignore();
	return 0;
}


