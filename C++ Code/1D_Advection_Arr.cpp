
///////////////////////////////////////////////////////////////////////
// 1D Advection Equation
// Carolyn Wendeln
// 10/24/2021
///////////////////////////////////////////////////////////////////////


#include <iostream>
#include <fstream>
#include <string>
#include <cmath>

using std::cout; using std::endl;
using std::to_string;
using namespace std;


///////////////////////////////////////////////////////////////////////
// Print Array Function Functions
///////////////////////////////////////////////////////////////////////

void print_array(double arr[], int n) 
{ 
	for (int i = 0; i < n; i++) 
		cout << arr[i] << ", "; 
} 

void output_array(double arr[], int N, int n)
{
  ofstream file;
  string file_name = "time_step_" + to_string(n) + ".txt";
  file.open (file_name);
  if (file.is_open())
  {
  file << "a_write = [";
  for (int i = 0; i < N; i++)
  	file << arr[i] << ", ";
  file << "] \n";
  file.close();
  }
  else cout << "can not open a file";
}


///////////////////////////////////////////////////////////////////////
// Main Function
///////////////////////////////////////////////////////////////////////

int main () {

	/////////////////////////////////////////////////////////////////////
	//Input Data
	/////////////////////////////////////////////////////////////////////


	double x_start = 0;
	double x_end = 1;
	int N = 100;
	double delta_x = (x_end-x_start) / (N-1);

	double t_start = 0;
	double t_end = 10;
	int n = 10;
	double delta_t = (t_end-t_start) / (n-1);

	double u = 0.1;

	/////////////////////////////////////////////////////////////////////
	//Generate Arrays
	/////////////////////////////////////////////////////////////////////
  
  double x_dist[N];
  double temp[N];
  double a_read[N];
	double a_write[N];

	/////////////////////////////////////////////////////////////////////
	//Initialize Arrays
	/////////////////////////////////////////////////////////////////////

  for(int i=0; i<(N); ++i)
  {
  	x_dist[i] = x_start + delta_x * i;
    a_read[i] = exp(-1 * pow(x_dist[i] - 0.5,2) * 50.0);
  }

	/////////////////////////////////////////////////////////////////////
	//Apply Advection Equation
	/////////////////////////////////////////////////////////////////////
  	
 	for(int j=0; j<n; ++j)
  {
  	//first case
  	a_write[0] = a_read[0] - (((u * delta_t) / (2 * delta_x)) * (a_read[1] - a_read[99]));
  	
  	//bulk
  	for(int i=1; i<(N-2); ++i)
   		{
   			a_write[i] = a_read[i] - (((u * delta_t) / (2 * delta_x)) * (a_read[i+1] - a_read[i-1]));
  		}
  	
  	//last case
  	a_write[99] = a_read[99] - (((u * delta_t) / (2 * delta_x)) * (a_read[0] - a_read[98]));

  	for(int i=1; i<N; ++i)
  	{
  		temp[i] = a_read[i];
  		a_read[i] = a_write[i];
  		a_write[i] = temp[i];
  	}
  	
  	output_array(a_write,N,j);
  }

}

