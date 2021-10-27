
///////////////////////////////////////////////////////////////////////
// 1D Advection Equation
// Carolyn Wendeln
// 10/24/2021
///////////////////////////////////////////////////////////////////////


#include <iostream>
using std::cout; using std::endl;
#include <vector>
using std::vector;
#include<cmath>


///////////////////////////////////////////////////////////////////////
// Print Array Function
///////////////////////////////////////////////////////////////////////

void print_array(double arr[], int n) 
{ 
	for (int i = 0; i < n; i++) 
		cout << arr[i] << ", "; 
} 

///////////////////////////////////////////////////////////////////////
// Main Function
///////////////////////////////////////////////////////////////////////

int main () {

	//////////////////////////////////////////////////////////////////
	//Input Data
	//////////////////////////////////////////////////////////////////

	double x_start = 0;
	double x_end = 1;
	int N = 100;
	double delta_x = (x_end-x_start) / N;
  	//double delta = (end - start) / (num - 1);


	double t_start = 0;
	double t_end = 10;
	int n = 100;
	double delta_t = (t_end-t_start) / n;

	double u = 0.1;

	//cout << delta_x << endl;
	//cout << delta_t << endl;

	//////////////////////////////////////////////////////////////////
	//Generate Arrays
	//////////////////////////////////////////////////////////////////
  
	double t_dist[n];
  double x_dist[N];

  double temp[N];
  double a_read[N];
	double a_write[N];

	//////////////////////////////////////////////////////////////////
	//Initialize Arrays
	//////////////////////////////////////////////////////////////////

  for(int i=0; i<(N); ++i)
  {
  	x_dist[i] = x_start + delta_x * i;
    a_read[i] = exp(-1 * pow(x_dist[i] - 0.5,2) * 50.0);
  }

  for(int j=0; j<(n); ++j)
  {
    t_dist[j] = t_start + delta_t * j;
  }

	//print_array(x_dist,N);
	//print_array(t_dist,n);
	//print_array(a_read,N);


	//////////////////////////////////////////////////////////////////
	//Apply Advection Equation
	//////////////////////////////////////////////////////////////////
  	
 	for(int j=0; j<n; ++j)

  {
  	
  	a_write[0] = a_read[0] - (((u * delta_t) / (2 * delta_x)) * (a_read[1] - a_read[99]));
  	
  	for(int i=1; i<(N-2); ++i)
   		{
   			a_write[i] = a_read[i] - (((u * delta_t) / (2 * delta_x)) * (a_read[i+1] - a_read[i-1]));
  		}
  	
  	a_write[99] = a_read[99] - (((u * delta_t) / (2 * delta_x)) * (a_read[0] - a_read[98]));

  	for(int i=1; i<N; ++i)
  	{
  		temp[i] = a_read[i];
  		a_read[i] = a_write[i];
  		a_write[i] = temp[i];
  	}




  }

	print_array(a_write,N);



}

