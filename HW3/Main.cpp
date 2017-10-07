#include <iostream>
#include <math.h>
using namespace std;


// This is the implementation of the scatter to the processes
// TODO: Must also send the number of multiples and probably the current column to the calculation so that it can perform the correct calculations
// TODO: Implement the dot product for each
// TODO: Implement the gather routine and assemble the final matrix

int a[6]  = { 1, 2, 3, 4, 5, 6 }; // 2 x 3
int b[12] = { 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12 }; // 3 x 4

//int a[9] = { 1, 2, 3, 4, 5, 6, 7, 8, 9 }; // 3 x 3
//int b[6] = { 1, 2, 3, 4, 5, 6 }; // 3 x 2

//int a[32] = { 1, 2, 3, 4, 5, 6, 7, 8, 9 , 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 26, 27, 28, 29, 30, 31, 32 }; // 5 x 6
//int b[48] = { 1, 2, 3, 4, 5, 6, 7, 8, 9 , 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48 }; // 6 x 8

//int a[2] = { 1, 2 }; // 1 x 2
//int b[6] = { 1, 2, 3, 4, 5, 6 }; // 2 x 3

//int a[6] = { 1, 2, 3, 4, 5, 6 }; // 3 x 2
//int b[2] = { 1, 2 }; // 2 x 1

int dim_l = 2;
int dim_m = 3;
int dim_n = 4;


int main(int argc, char* argv[])
{
	//int num_mults = dim_n*dim_l;
  int num_mults;
	int numtasks = 3;

  // Set the number of multiplies required to get to the result
  if (dim_l == 1)
    num_mults = dim_n;
  else if (dim_n == 1)
    num_mults = dim_l;
  else
    num_mults = dim_l * dim_n;

	int row_ind = 0;
	int col_ind = 0;
	int local_a[200]; // going for maximum size, don't care about space it's cheap
	int local_b[200];

  // Variables to keep up with what row and column we are reading from
	int begin_row = 0;
	int begin_column = 0;

  // Loop over the MPI tasks
	for (int j = 0; j < numtasks; j++)
	{
    // Calculate the base number of multiplies for each task
    int base = num_mults / numtasks;

    // Calculate the extra if it is not an even distribution
    int extra = num_mults % numtasks;

    // If there are any extra assign them to the first tasks up to rank of extra
    if (extra != 0)
    {
      // If rank is less than the number of extra items, then this process gets an extra multiply to process
      if (j < extra)
        base = (num_mults / numtasks) + 1;
    }

    // Each task gets at least one row to work with
    int num_rows = 1;

    // Variables to keep up with navigating the matrices
    int count = base;
    int curr_col = begin_column;

    // Calculate the number of rows that we need to send to the process.
    while (count > 0)
    {
      // Get the current distance from the end of dim_n
      int diff = dim_n - curr_col;

      // If we have more multiplies to process then left for this row, we must add another row
      if (count > diff)
      {
        num_rows++;
        count -= diff;
        curr_col = 0;
        continue;
      }

      count -= dim_n;
    }

    // Reset our column check and multiply count
    curr_col = begin_column; 
    count = base;

		// Row assignment
    for(int r = 0; r < num_rows; r++)
		{
      // Fill up the local matrix with values from the main A matrix
			for (int t = 0; t < dim_m; t++)
			{
				local_a[r*dim_m + t] = a[begin_row*dim_m + t];
			}

      // Calculate the distance from the end of dim_n for this set of row multiplications
      int diff = dim_n - curr_col;

      // If we still have more to process we must add the next row to this processes variables
      if (diff <= count)
      {
        begin_row++;
        count -= diff;
        curr_col = 0;
      }
		}

		// Column Assignment
    // This will put the columns in a matrix starting with the last column used for a multiplication
    // This could be cleaned up by broadcasting the entire matrix to each process but I feel that this might be faster with larger 
    // matrices as there won't be as much duplication.
		for (int i = 0; i < base; i++)
		{
      // Fill up the local matrix with values from the B matrix
			for (int k = 0; k < dim_m; k++)
			{
				local_b[i*dim_m + k] = b[begin_column + k*dim_n];
			}

      // Keep up with the current dim_n column that we processed
			if (begin_column != 0 && (begin_column % (dim_n-1) == 0))
				begin_column = 0;
			else if ((begin_column + 1) != dim_n) // can't exceed the dim_n for current column
				begin_column++;
		}
	}

	return 0;
}