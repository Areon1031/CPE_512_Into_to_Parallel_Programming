// Travelling Salesman Program -- Serial Version
// B. Earl Wells -- November 2017
//
// Travelling Salesman Problem: "Given a specified number of "cities"
// along with the cost of travel between each pair of them, find the 
// cheapest way of visiting all the cities and returning to the first 
// city visited."
//
// For an n city tour this program examines all (n-1)! tours and 
// returns a tour with the least travel cost.
//
using namespace std;
#include <iostream>
#include <iomanip>
#include <stdlib.h>
//#include <sys/time.h>
#include <time.h>

#define TIMER_CLEAR     (tv1.tv_sec = tv1.tv_usec = tv2.tv_sec = tv2.tv_usec = 0)
#define TIMER_START     gettimeofday(&tv1, (struct timezone*)0)
#define TIMER_ELAPSED   ((long long) (tv3.tv_usec)+((long long) (tv3.tv_sec)*1000000))
#define TIMER_STOP      {gettimeofday(&tv2, (struct timezone*)0);timersub(&tv2,&tv1,&tv3);}
//struct timeval tv1,tv2,tv3;

#define MX_CITIES 30       // maximum number of cities in tour
#define TRUE 1
#define FALSE 0
#define SEED (long) 178937 // random number seed
#define MAX_INT 0x7fffffff // maximum integer possible

// Defines so that I can compile the code in visual studio
#define srand48(s) srand(s)
#define drand48() (((double)rand())/((double)RAND_MAX))

// generate randomly a cost matrix that represents the cost of travel
// between each two cities. This generates a nonsymetric cost matrix
// which means the cost of travel between any two cities may be 
// depend on which city is the source and which is the destination
void fill_cost_matrix(int cost_matrix[][MX_CITIES],int num_cities) {
   int i,j;
   srand48(SEED);
   if (num_cities<=MX_CITIES) {
      for (i=0;i<num_cities;i++) {
         for (j=0;j<num_cities;j++) {
            if (i==j) cost_matrix[i][j]=0;
            else cost_matrix[i][j]= (int) (drand48()*100);
         }
      }
   }
   else {
      cout << "Error: Too many cities -- increase MX_CITIES parameter" << endl;
      cout << " and recompile" << endl;
      exit(1);
   }
}

// This routine outputs to the screen the city cost matrix that has
// been generated 
void print_cost_matrix(int cost_matrix[][MX_CITIES],int num_cities) {
   int i,j;
   for (i=0;i<num_cities;i++) {
      for (j=0;j<num_cities;j++) {
         cout << setw(5) << cost_matrix[i][j]; 
      }
      cout << endl;
   }
   cout << endl;
}

// This routine outputs the City visit order from left to right
// assuming that we are always starting our tour with City 0.
// Note: there is no loss of generality in this since any tour
// order can be rotated to make any of the cities the start and
// end point.
void print_city_visit_order(unsigned int city_order[],unsigned int num_cities_m1) {
   int i;
   cout << "[City 0]->";
   for (i=0;i<num_cities_m1;i++) {
      cout << "[City " << city_order[i]+1 << "]->";
   }
   cout << "[City 0]" << endl;
}

// routine that computes the cost of each tour where the tour is
// described in the city_order array that should have a valid 
// city order premutation of the n-1 remaining cities after the
// first city, City 0, is assumed as the starting point. The
// tour starts at City 0 and goes to the first city in the city_order
// array. The tour is then processed in the sequence dictated by
// the city_order permutation after which the tour proceeds from the
// last city in this permutation back to City 0.
// Note: there is no loss of generality in this since any tour
// order can be rotated to make any of the cities the start and
// end point.
int tour_cost(unsigned int num_cities,unsigned int *city_order,
              int cost_matrix[MX_CITIES][MX_CITIES]) {
   int i,last_city,cost;
   last_city=city_order[0]+1;
   cost = cost_matrix[0][last_city];
   for (i=1;i<num_cities-1;i++) {
      cost += cost_matrix[last_city][city_order[i]+1];
      last_city = city_order[i]+1;
   }
   cost += cost_matrix[last_city][0];
   return cost;
}

// Routine to save the current city_order permutation to another
// data structure so that it can be displayed later by the 
// print_city_visit_order routine
void save_order (unsigned int *city_order_sv,unsigned int num_cities,
                 unsigned int * city_order) {
   int i;
   for (i=0;i<num_cities;i++) {
      city_order_sv[i]=city_order[i];
   }
}

// data structure that is used to store the best city visit order
unsigned int best_order[MX_CITIES],best_cost=MAX_INT; //set to worst value poss.

unsigned long int num_perms(unsigned int num) {
   unsigned long int nm_perms = 1;
   for (unsigned int i = 1; i <= num; i++) {
      nm_perms *= (unsigned long int) i;
   }
   return nm_perms;
}
 
// routine to determine a unique permutation ordering that is
// associated with the specified permutation number assuming
// the specified number of elements that are to be permuted
// i.e. num_symbols. 
// See https://en.wikipedia.org/wiki/Factorial_number_system
// for more information concerning the method that was employed
void permutation(unsigned int *perm, unsigned long int perm_num, 
                 unsigned int num_symbols) {
   unsigned int factoradic[MX_CITIES];
   for (unsigned int i = 1; i <= num_symbols; i++) {
      factoradic[i-1] = perm_num%(long) i;
      perm_num /= (long) i;
   }

   // compute lexical order
   {
      int perm_mask[MX_CITIES];
      for (int i = 0; i<num_symbols; i++) {
         perm_mask[i]=0;
      }
      for (int i = num_symbols-1; i>=0; i--) {
         int lex_cnt = -1; 
         for (int j = 0; j<num_symbols;j++) {
            if (perm_mask[j]==0) lex_cnt++; 
            if (lex_cnt==factoradic[i]) {
               perm_mask[j]=1;
               perm[i] = j;
               break;
            }
         }
      }
   } 
}   

// krr good candidate for parallelizing
void tour_search(unsigned int nm_cities,int city_cost_matrix[][MX_CITIES]) {
   // create a permutation matrix that list the cities in tour visit order
   // This one dimensional matrix is called perm_order and it holds the
   // index number of the nm_cities in the tour minus one because one
   // city is assumed to be held in place with the other city order being
   // permuted around it. Note that this program assumes that the city 
   // held constant is city 0, and this city starts and ends the tour cycle.
   // This means that city 0 will preceed the city that is placed in slot 0 
   // of the perm_order matrix and this city will be the final destination
   // city that is present in slot nm_cities - 2
   unsigned long int total_tours=num_perms(nm_cities-1); 
   for (unsigned long int tour_nm=0;tour_nm<total_tours;tour_nm++)
   {
      unsigned int perm_order[MX_CITIES];
      // get next permutation 
      permutation(perm_order, tour_nm, nm_cities-1);

      int tour_cst=tour_cost(nm_cities,perm_order,city_cost_matrix);
      if (tour_cst<best_cost) {
         save_order (best_order,nm_cities-1,perm_order);
         best_cost = tour_cst;
      }
      // uncomment to view all city tours
      // print_city_visit_order(perm_order,nm_cities-1);
      // cout << "tour cost = " << tour_cst << endl;
   }
}

// routine to obtain the number of cities from the user
// via the command line or by prompting the user for input
unsigned int get_city_number(int argc,char *argv[]) {
   int num_cities;
   if (argc==2) {
      num_cities = atoi(argv[1]);
   }
   else {
      if (argc==1) { 
         // input number of cities 
         cout << "Enter number of cities:" << endl;
         cin >> num_cities;
         cout << endl; 
      }
      else {
         cout << "usage: tsp_serial <number of cities>" << endl;
         exit(1);
      }
   }
   if ((num_cities<2) || (num_cities>MX_CITIES)) {
      cout << "Error: Number of Cities too large or too small" << endl;
      exit(1);
   }
   return (unsigned int) num_cities;
}
 
int main(int argc, char *argv[]) {

   int city_cost_matrix[MX_CITIES][MX_CITIES];
   unsigned int num_cities;

   // get number of cities from the user
   num_cities = get_city_number(argc,argv);

   // fill city cost matrix
   fill_cost_matrix(city_cost_matrix,num_cities);
   
   // print city cost matrix
   print_cost_matrix(city_cost_matrix,num_cities);

   //TIMER_CLEAR;
   //TIMER_START;

   // search through all possible orderings of the
   // cities
   tour_search(num_cities,city_cost_matrix);

   //TIMER_STOP;
   // print city visit order
   cout << "Best City Tour:"  << endl; 
   print_city_visit_order(best_order,num_cities-1);

   cout << "tour cost = "  << 
      tour_cost(num_cities,best_order,city_cost_matrix)
      << endl;
  // cout << num_cities << "  " <<  TIMER_ELAPSED/1000000.0 << endl;
   cin.ignore();
   return 0;
}

