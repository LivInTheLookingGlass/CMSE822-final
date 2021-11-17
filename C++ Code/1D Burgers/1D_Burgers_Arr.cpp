
///////////////////////////////////////////////////////////////////////
// 1D Burgers Equation
// Carolyn Wendeln
// 11/05/2021
///////////////////////////////////////////////////////////////////////


#include <iostream>
#include <fstream>
#include <string>
#include <cmath>

using std::cout; using std::endl;
using std::to_string;
using namespace std;


///////////////////////////////////////////////////////////////////////
//Output Array Functions
///////////////////////////////////////////////////////////////////////

void print_array(double* arr, int n) 
{ 
	for (int i = 0; i < n; i++) 
		cout << arr[i] << ", "; 
} 

void output_array(double* arr, int N, int n, double* x_dist)
{
  ofstream file;
  string file_name = "1D_time_step_" + to_string(n) + ".csv";
  file.open (file_name);
  if (file.is_open())
  {
  	file << "x dist" << " ," << "func value" << endl;

  	for (int i = 0; i < N; i++)
  		file << x_dist[i] << " ," << arr[i] << endl;
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
	int N = 10;
	double delta_x = (x_end-x_start) / (N-1);

	double t_start = 0;
	double t_end = 10;
	int n = 40;
	double delta_t = (t_end-t_start) / (n-1);

	double u = 0.1;

	/////////////////////////////////////////////////////////////////////
	//Generate Arrays
	/////////////////////////////////////////////////////////////////////

	double* x_dist = new double[N];
	double* a_read = new double[N];
	double* a_write = new double[N];
	double* temp; 

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
  	
 	for(int t=0; t<n; ++t)
  {
  	//first case
  	a_write[0] = a_read[0] - (((a_read[0] * delta_t) / (2 * delta_x)) * (a_read[1] - a_read[N-1]));

  	//bulk
  	for(int i=1; i<(N-2); ++i)
   		{
   			a_write[i] = a_read[i] - (((a_read[i] * delta_t) / (2 * delta_x)) * (a_read[i+1] - a_read[i-1]));
  		}
  	
  	//last case
  	a_write[N-1] = a_read[N-1] - (((a_read[N-1] * delta_t) / (2 * delta_x)) * (a_read[0] - a_read[N-2]));

  	temp = a_read;
    a_read = a_write;
    a_write = temp;

  	output_array(a_write,N,t,x_dist);
  }


	delete [] x_dist;
  delete [] a_read;
	delete [] a_write;
}

