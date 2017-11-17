#include <iostream>
#ifdef _OPENMP
#include <omp.h>
#endif
using namespace std;

int test[100][100];

int main(int argc, char* argv[])
{
  int count[4] = { 0, 0, 0, 0 };
#ifdef _OPENMP
  double start = omp_get_wtime();
#endif
#pragma omp parallel for num_threads(4) schedule(dynamic, 20)
  for (int i = 0; i < 99; i++)
  {
    //count[omp_get_thread_num()]++;
    for (int j = i+1; j < 100; j++)
    {
      count[omp_get_thread_num()]++;
      test[i][j] = 42;
    }
  }

#ifdef _OPENMP
  double finish = omp_get_wtime();
  double time = finish - start;
  cout << "Time: " << time << " seconds" << endl;
#endif

  for (int i = 0; i < 4; i++)
  {
    cout << "Count for thread " << i << ": " << count[i] << endl;
  }
  cin.ignore();
  return 0;
}