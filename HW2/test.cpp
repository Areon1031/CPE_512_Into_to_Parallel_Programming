//test.cpp

#include <iostream>
#include <stdlib.h>
using namespace std;



int main(int argc, char* argv[])
{
	int numtasks = atoi(argv[1]);
	int data_size = atoi(argv[2]);

	int base = data_size / numtasks;
	int extra = data_size % numtasks;

	int rank[4];
	int* ranksize = new int[numtasks];

	int* group = new int[numtasks];

	for (int i = 0; i < numtasks; i++)
	{
		rank[i] = i;

		group[i] = rank[i] < extra ? base + 1 : base;
	}

	cout << "Base: " << base << endl;
	cout << "Extra: " << extra << endl;

	for (int i = 0; i < numtasks; i++)
		cout << "Rank[" << i << "]: " << group[i] << endl;

	delete[] ranksize;
	delete[] group;
	return 0;
}
