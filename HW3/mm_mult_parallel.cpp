/******************************************************************/
/* Matrix Matrix Multiplication Program Example -- serial version */
/* September 2016 -- B. Earl Wells -- University of Alabama       */
/*                                    in Huntsville               */
/******************************************************************/
// mm_mult_serial.cpp
// compilation:
//   gnu compiler
//      g++ mm_mult_serial.cpp -o mm_mult_serial -O3 -lm
// Note: to compile a parallel MPI program version which is named 
//   mm_mult_mpi.cpp
//   then execute the following command
//      gnu compiler
//         mpic++ mm_mult_mpi.cpp -o mm_mult_MPI_gnu -lm  -O3
/*
   This program is designed to perform matrix matrix multiplication
   A x B = C, where A is an lxm matrix, B is a m x n matrix and
   C is a l x n matrix. The program is designed to be a template
   serial program that can be expanded into a parallel multiprocess
   and/or a multi-threaded program.

   The program randomly assigns the elements of the A and B matrix
   with values between 0 and a MAX_VALUE. It then multiples the
   two matrices with the result being placed in the C matrix.
   The program prints out the A, B, and C matrices.

   The program is executed using one or three command line parameters.
   These parameters represent the dimension of the matrices. If only
   one parameter is used then then it is assumed that square matrices are
   to be created and multiplied together that have the specified
   dimension. In cases where three command line parameters are entered
   then the first parameter is the l dimension, the second the m, and
   the third is the n dimension.

   To execute:
   mm_mult_serial [l_parameter] <m_parameter n_parameter>
*/

using namespace std;
#include <iostream>
#include <iomanip>
#include <sstream>
#include <stdlib.h>
#include <string.h>
//#include <sys/time.h>
#include <time.h>
#include <mpi.h>

#define MX_SZ 320
#define SEED 2397           /* random number seed */
#define MAX_VALUE  100.0    /* maximum size of array elements A, and B */

/* copied from mpbench */
#define TIMER_CLEAR     (tv1.tv_sec = tv1.tv_usec = tv2.tv_sec = tv2.tv_usec = 0)
#define TIMER_START     gettimeofday(&tv1, (struct timezone*)0)
#define TIMER_ELAPSED   ((tv2.tv_usec-tv1.tv_usec)+((tv2.tv_sec-tv1.tv_sec)*1000000))
#define TIMER_STOP      gettimeofday(&tv2, (struct timezone*)0)
//struct timeval tv1,tv2;

// Defines so that I can compile the code in visual studio
#define srand48(s) srand(s)
#define drand48() (((double)rand())/((double)RAND_MAX))


/*
This declaration facilitates the creation of a two dimensional
dynamically allocated arrays (i.e. the lxm A array, the mxn B
array, and the lxn C array).  It allows pointer arithmetic to
be applied to a single data stream that can be dynamically allocated.
To address the element at row x, and column y you would use the
following notation:  A(x,y),B(x,y), or C(x,y), respectively.
Note that this differs from the normal C notation if A were a
two dimensional array of A[x][y] but is still very descriptive
of the data structure.
*/
float *a, *b, *c;
#define A(i,j) *(a+i*dim_m+j)
#define B(i,j) *(b+i*dim_n+j)
#define C(i,j) *(c+i*dim_n+j)

/*
   Routine to retrieve the data size of the numbers array from the
   command line or by prompting the user for the information
*/
void get_index_size(int argc, char *argv[], int *dim_l, int *dim_m, int *dim_n, int rank) {
  if (argc != 2 && argc != 4) {
    if (rank == 0)
    {
      cout << "usage:  mm_mult_serial [l_dimension] <m_dimension n_dimmension>"
        << endl;
      MPI_Finalize();
      exit(1);
    }
  }
  else {
    if (argc == 2) {
      *dim_l = *dim_n = *dim_m = atoi(argv[1]);
    }
    else {
      *dim_l = atoi(argv[1]);
      *dim_m = atoi(argv[2]);
      *dim_n = atoi(argv[3]);
    }
  }
  if (rank == 0)
  {
    if (*dim_l <= 0 || *dim_n <= 0 || *dim_m <= 0) {
      cout << "Error: number of rows and/or columns must be greater than 0"
        << endl;
      MPI_Finalize();
      exit(1);
    }
  }
}

/*
   Routine that fills the number matrix with Random Data with values
   between 0 and MAX_VALUE
   This simulates in some way what might happen if there was a
   single sequential data acquisition source such as a single file
*/
void fill_matrix(float *array, int dim_m, int dim_n)
{
  int i, j;
  for (i = 0; i < dim_m; i++) {
    for (j = 0; j < dim_n; j++) {
      array[i*dim_n + j] = drand48()*MAX_VALUE;
    }
  }
}

/*
   Routine that outputs the matrices to the screen
*/
void print_matrix(float *array, int dim_m, int dim_n)
{
  int i, j;
  for (i = 0; i < dim_m; i++) {
    for (j = 0; j < dim_n; j++) {
      cout << array[i*dim_n + j] << " ";
    }
    cout << endl;
  }
}

/*
ONE-TO-ALL BROADCAST COMMUNICATION ROUTINE
Routine to transfer from the root MPI process the value of
the 'int_num' parameter to all other MPI processes in the system.
*/
void broadcast_int(int *int_num, int root, int rank, int numtasks) {
  MPI_Status status;

  int type = 123;

  // root send value of int_num to each of the other processes
  // using a locally blocking point-to-point send
  if (rank == root) {
    for (int mpitask = 0; mpitask < numtasks; mpitask++) {
      if (mpitask != root) {
        MPI_Send(int_num, 1, MPI_INT,
          mpitask, type, MPI_COMM_WORLD);
      }
    }
  }
  // if not root process execute a blocking point-to-point receive
  // with the source being to root process and direct this data to
  // the local copy of 'int_num'
  else {
    MPI_Recv(int_num, 1, MPI_INT,
      root, type, MPI_COMM_WORLD, &status);
  }

}

/* ONE-TO-ALL SCATTER ROUTINE
Routine to divide and scatter the number data array that resides on the
root MPI process to all other MPI processes in the system.
The number data size is given by the'num_size' parameter its source
address is given by the '*numbers' parameter, and the destination
group data associated with the current process is given by the
'*group' parameter.  */
void scatter(float* a, float* b, float *group_a, float* group_b, int num_size, int root, int rank, int numtasks, int dim_l, int dim_m, int dim_n, int* begin_column)
{
  MPI_Status status;
  int type = 234;

  int num_mults = 0;
  if (dim_l == 1)
    num_mults = dim_n;
  else if (dim_n == 1)
    num_mults = dim_l;
  else
    num_mults = dim_l * dim_n;

  int local_a[200];
  int local_b[200];

  int row_ind = 0;
  int col_ind = 0;

  // Variables to keep up with what row and column we are reading from
  int begin_row = 0;
  begin_column = 0;

  if (rank == root)
  {
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
      int curr_col = *begin_column;

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
      curr_col = *begin_column;
      count = base;

      // Row assignment
      for (int r = 0; r < num_rows; r++)
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
          local_b[i*dim_m + k] = b[*begin_column + k*dim_n];
        }

        // Keep up with the current dim_n column that we processed
        if (*begin_column != 0 && (*begin_column % (dim_n - 1) == 0))
          *begin_column = 0;
        else if ((*begin_column + 1) != dim_n) // can't exceed the dim_n for current column
          (*begin_column)++;
      }

      // This should work for every other process, but I need to make sure that the root keeps what it needs
      if (j != root)
      {
        // Send the data to the other processes
        MPI_Send(local_a, num_rows*dim_m, MPI_FLOAT, j, type, MPI_COMM_WORLD);
        MPI_Send(local_b, base*dim_m, MPI_FLOAT, j, type, MPI_COMM_WORLD);
        MPI_Send(begin_column, 1, MPI_INT, j, type, MPI_COMM_WORLD);
      }
    }
  }
  else
  {
    MPI_Recv(a, dim_l*dim_m, MPI_FLOAT, root, type, MPI_COMM_WORLD, &status);
    MPI_Recv(b, dim_l*dim_m, MPI_FLOAT, root, type, MPI_COMM_WORLD, &status);
    MPI_Recv(begin_column, 1, MPI_INT, root, type, MPI_COMM_WORLD, &status);
  }
}
/*
ALL-TO-ONE Reduce ROUTINE
Routine to accumulate the result of the local summation associated
with each MPI process. This routine takes these partial sums and
produces a global sum on the root MPI process (0)
Input arguments to routine include variable name of local partial
sum of each MPI process. The function returns to MPI root process 0,
the global sum (summation of all partial sums).
*/
void reduce(float* sum, float* partial_sum, int root, int rank, int numtasks)
{
  MPI_Status status;
  int type = 123;
  // if MPI root process sum up results from the other p-1 processes
  if (rank == root) {
    *sum = *partial_sum;
    for (int mpitask = 0; mpitask < numtasks; mpitask++) {
      if (mpitask != root) {
        MPI_Recv(partial_sum, 1, MPI_DOUBLE,
          mpitask, type, MPI_COMM_WORLD, &status);
        (*sum) += (*partial_sum);
      }
    }
  }
  // if not root MPI root process then send partial sum to the root  
  else {
    MPI_Send(partial_sum, 1, MPI_DOUBLE,
      root, type, MPI_COMM_WORLD);
  }
}

/*
   MAIN ROUTINE: summation of a number list
*/

int main(int argc, char *argv[])
{
  float dot_prod;
  int dim_l, dim_n, dim_m;
  int i, j, k;

  int num_mults, group_size, num_group, i;
  int numtasks, rank, num;
  int start_column;
  MPI_Status status;

  // Main Routine

  MPI_Init(&argc, &argv); // initalize MPI environment
  MPI_Comm_size(MPI_COMM_WORLD, &numtasks); // get total number of MPI processes
  MPI_Comm_rank(MPI_COMM_WORLD, &rank); // get unique task id number 


  // get matrix sizes
  get_index_size(argc, argv, &dim_l, &dim_m, &dim_n, rank);

  // The root process fills the matrices and then passes them to the othe processes
  if (rank == 0)
  {
    // dynamically allocate from heap the numbers in the memory space
    // for the a,b, and c matrices 
    a = new (nothrow) float[dim_l*dim_m];
    b = new (nothrow) float[dim_m*dim_n];
    c = new (nothrow) float[dim_l*dim_n];
    if (a == 0 || b == 0 || c == 0)
    {
      cout << "ERROR:  Insufficient Memory" << endl;
      exit(1);
    }

    start_column = 0;

    /*
       initialize numbers matrix with random data
    */
    srand48(SEED);
    fill_matrix(a, dim_l, dim_m);
    fill_matrix(b, dim_m, dim_n);

    /*
      output numbers matrix
    */
    cout << "A matrix =" << endl;
    print_matrix(a, dim_l, dim_m);
    cout << endl;

    cout << "B matrix =" << endl;
    print_matrix(b, dim_m, dim_n);
    cout << endl;


    // Broadcast the number of multiplies to each process.
    //broadcast_int(&num_mults, 0, 0, numtasks);
    //MPI_Bcast(&num_mults, 1, MPI_INT, 0, MPI_COMM_WORLD);
  }

  // broad cast the data size, which is really the number of multiplies
  if (dim_l == 1)
    num_mults = dim_n;
  else if (dim_n == 1)
    num_mults = dim_l;
  else
    num_mults = dim_l * dim_n;

  int base = num_mults / numtasks;
  int extra = num_mults % numtasks;

  if (rank > 0)
  {
    // Just need to allocate enough space to hold the rows that the process will receive
    // The below will need to change, but this is fine for now.


    // Each process granted that num_mults >= numtasks should have at least one row and one column to perform a dot product on.
    //a = new (nothrow) float[rank < extra ? base*dim_m : (base - 1)*dim_m];
    //b = new (nothrow) float[rank < extra ? base*dim_m : (base - 1)*dim_m];
    a = new (nothrow) float[dim_l * dim_n]; // don't care about the size right now, memory is cheaper than run time
    b = new (nothrow) float[dim_l * dim_n]; 
    c = new (nothrow) float[dim_l * dim_n];

    start_column = 0;

    if (a == 0 || b == 0)
    {
      cout << "ERROR:  Insufficient Memory" << endl;
      exit(1);
    }
  }

  // Scatter the Data
  // The root process needs to scatter the correct amount of data to each process.
  scatter(a, b, 1, 0, rank, numtasks, dim_l, dim_m, dim_n, start_column);


  // Each process will start working on the data here

  /*
  Start recording the execution time
  */
  /*TIMER_CLEAR;
  TIMER_START;*/

  // multiply local part of matrix
  for (i = 0; i < dim_l; i++) {
    for (j = 0; j < dim_n; j++) {
      dot_prod = 0.0;
      // Need to use this loop to send the arrays to both for the dot product
      for (k = 0; k < dim_m; k++) {
        dot_prod += A(i, k)*B(k, j);
      }
      C(i, j) = dot_prod;
    }
  }

  // Gather, each process will send the results of the dot product back to the main routine
  // The main routine will then put the final lxn matrix back together.

  /*
     stop recording the execution time
  */
  //TIMER_STOP;

  cout << "C matrix =" << endl;
  print_matrix(c, dim_l, dim_n);
  cout << endl;
  //cout << "time=" << setprecision(8) <<  TIMER_ELAPSED/1000000.0 
       //<< " seconds" << endl;


  // KRR 
  cin.ignore();
}


