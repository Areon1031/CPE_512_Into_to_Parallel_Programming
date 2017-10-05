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
float *a,*b,*c;
#define A(i,j) *(a+i*dim_m+j)
#define B(i,j) *(b+i*dim_n+j)
#define C(i,j) *(c+i*dim_n+j)

/*
   Routine to retrieve the data size of the numbers array from the 
   command line or by prompting the user for the information
*/
void get_index_size(int argc,char *argv[],int *dim_l,int *dim_m,int *dim_n) {
   if(argc!=2 && argc!=4) {
      cout<<"usage:  mm_mult_serial [l_dimension] <m_dimension n_dimmension>"
           << endl;
      exit(1);
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
   if (*dim_l<=0 || *dim_n<=0 || *dim_m<=0) {
      cout<<"Error: number of rows and/or columns must be greater than 0"
          << endl;
      exit(1);
   }
}

/*
   Routine that fills the number matrix with Random Data with values
   between 0 and MAX_VALUE
   This simulates in some way what might happen if there was a 
   single sequential data acquisition source such as a single file
*/
void fill_matrix(float *array,int dim_m,int dim_n)
{
   int i,j;
   for(i=0;i<dim_m;i++) {
      for (j=0;j<dim_n;j++) {
         array[i*dim_n+j]=drand48()*MAX_VALUE;
      }
   }
}

/*
   Routine that outputs the matrices to the screen 
*/
void print_matrix(float *array,int dim_m,int dim_n)
{
   int i,j;
   for(i=0;i<dim_m;i++) {
      for (j=0;j<dim_n;j++) {
         cout << array[i*dim_n+j] << " ";
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
void scatter(float *numbers, float* group, int num_size, int root, int rank, int numtasks)
{
  MPI_Status status;
  int type = 234;

  // determine number of elements in subarray groups to be processed by
  // each MPI process assuming a perfectly even distribution of elements 
  int number_elements_per_section = num_size / numtasks;

  // if root MPI process send portion of numbers array to each of the
  // the other MPI processes as well as make a copy of the portion
  // of the numbers array that is slated for the root MPI process
  if (rank == root) {
    int begin_element = 0;

    for (int mpitask = 0; mpitask < numtasks; mpitask++) {

      // in MPI root process case just copy the appropriate subsection
      // locally from the numbers array over to the group array
      if (mpitask == root) {
        for (int i = 0; i < number_elements_per_section; i++)
          group[i] = numbers[i + begin_element];
      }
      // if not the root process send the subsection data to
      // the next MPI process
      else {
        MPI_Send(&numbers[begin_element], number_elements_per_section,
          MPI_DOUBLE, mpitask, type, MPI_COMM_WORLD);
      }
      // point to next unsent or uncopied data in numbers array
      begin_element += number_elements_per_section;
    }
  }
  // if a non root process just receive the data
  else {
    MPI_Recv(group, number_elements_per_section, MPI_DOUBLE,
      root, type, MPI_COMM_WORLD, &status);
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

int main( int argc, char *argv[])
{
   float dot_prod;
   int dim_l,dim_n,dim_m;
   int i,j,k;

   int num_mults, group_size, num_group, i;
   int numtasks, rank, num;
   MPI_Status status;

   /* 
   get matrix sizes
   */
   get_index_size(argc,argv,&dim_l,&dim_m,&dim_n);

   // Main Routine

   MPI_Init(&argc, &argv); // initalize MPI environment
   MPI_Comm_size(MPI_COMM_WORLD, &numtasks); // get total number of MPI processes
   MPI_Comm_rank(MPI_COMM_WORLD, &rank); // get unique task id number 
  
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

     // broad cast the data size, which is really the number of multiplies
     num_mults = dim_l * dim_n;
     broadcast_int(&num_mults, 0, 0, numtasks);
   }
   else
   {
     // Just need to allocate enough space to hold the rows that the process will receive
     // The below will need to change, but this is fine for now.
     a = new (nothrow) float[dim_l*dim_m];
     b = new (nothrow) float[dim_m*dim_n];
     if (a == 0 || b == 0) 
     {
       cout << "ERROR:  Insufficient Memory" << endl;
       exit(1);
     }
   }


   // Scatter the Data
   scatter(a, b, 1, 0, rank, numtasks);
   // Each process will start working on the data here

   /*
   Start recording the execution time
   */
   /*TIMER_CLEAR;
   TIMER_START;*/

   // multiply local part of matrix
   for (i=0;i<dim_l;i++) {
      for (j=0;j<dim_n;j++) {
         dot_prod = 0.0;
         for (k=0;k<dim_m;k++) {
            dot_prod += A(i,k)*B(k,j);
         }
         C(i,j) = dot_prod;
      }
   }
   
   // Gather, each process will send the results of the dot product back to the main routine
   // The main routine will then put the final lxn matrix back together.

   /*
      stop recording the execution time
   */ 
   //TIMER_STOP;

   cout << "C matrix =" << endl;
   print_matrix(c,dim_l,dim_n);
   cout << endl;
   //cout << "time=" << setprecision(8) <<  TIMER_ELAPSED/1000000.0 
        //<< " seconds" << endl;


   // KRR 
   cin.ignore();
}


