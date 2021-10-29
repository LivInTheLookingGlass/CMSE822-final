
///////////////////////////////////////////////////////////////////////
// 2D Advection Equation
// Carolyn Wendeln
// 10/24/2021
///////////////////////////////////////////////////////////////////////


#include <iostream>
using std::cout; using std::endl;
#include<cmath>

//at + (ua)x + (va)y = 0

///////////////////////////////////////////////////////////////////////
// Output Array Functions
///////////////////////////////////////////////////////////////////////

void print_array(double** arr, int n) 
{ 
	cout << " x " << " y " << " a_write " << endl;
	for (int i = 0; i < n; i++)
	{
		for (int j = 0; j < n; j++)
		{
			cout << i << " " << j << " " << arr[i][j] << endl; 
		}
	}
} 
// change the way it prints
// x  y  value 
// x  y  value
	
///////////////////////////////////////////////////////////////////////
// Main Function
///////////////////////////////////////////////////////////////////////

int main () {

	//////////////////////////////////////////////////////////////////
	//Input Data
	//////////////////////////////////////////////////////////////////

	double dist_start = 0;
	double dist_end = 1;
	int N = 100;
	double delta_dist = (dist_end-dist_start) / (N-1);

	double t_start = 0;
	double t_end = 10;
	int n = 10;
	double delta_t = (t_end-t_start) / (n-1);

	double u = 0.1;

	//////////////////////////////////////////////////////////////////
	//Generate Arrays
	//////////////////////////////////////////////////////////////////
  
  	double dist_x[N];
 	double dist_y[N];

	//double** dist_x = new double*[N];
	//double** dist_y = new double*[N];

	double** temp = new double*[N];
	double** a_read = new double*[N];
	double** a_write = new double*[N];
	for(int i = 0; i < N; ++i)
	{
    	temp[i] = new double[N];
    	a_read[i] = new double[N];
    	a_write[i] = new double[N];

	}

	//////////////////////////////////////////////////////////////////
	//Initialize Arrays
	//////////////////////////////////////////////////////////////////

  	for(int i=0; i<(N); ++i)
  	{
  		dist_x[i] = dist_start + delta_dist * i;
  		dist_y[i] = dist_start + delta_dist * i;
  	}


  	// Initial condition
	for(int i=0; i<(N-1); ++i)
  	{
  		for(int j=0; j<(N-1); ++j)
  		{
  			a_read[i][j] = exp(-1 * (pow(dist_x[i] - 0.5,2) + pow(dist_y[j] - 0.5,2)) * 50.0);
  		}
 	 }

	//////////////////////////////////////////////////////////////////
	/////////////////////////////////////////////////////////////////
	//Apply Advection Equation
	//////////////////////////////////////////////////////////////////
  	//////////////////////////////////////////////////////////////////

 	 //loop through time steps
 	 for(int t=0; t<n; ++t)
 	 {

		/////////////////////////////////////////////////////////////////
		//four corners
		//////////////////////////////////////////////////////////////////

	  	//i = 0  ; j = 0
	  	a_write[0][0] = a_read[0][0] - ((((u * delta_t) / (2 * delta_dist))) * (((a_read[1][0] - a_read[99][0])) + ((a_read[0][1] - a_read[0][99]))));

		//i = 0  ; j = 99
	  	a_write[0][99] = a_read[0][99] - ((((u * delta_t) / (2 * delta_dist))) * (((a_read[1][99] - a_read[99][99])) + ((a_read[0][0] - a_read[0][98]))));

		//i = 99 ; j = 0
		 a_write[99][0] = a_read[99][0] - ((((u * delta_t) / (2 * delta_dist))) * (((a_read[0][0] - a_read[98][0])) + ((a_read[99][1] - a_read[99][99]))));

		//i = 99 ; j = 99
	  	a_write[99][99] = a_read[99][99] - ((((u * delta_t) / (2 * delta_dist))) * (((a_read[0][99] - a_read[98][99])) + ((a_read[99][0] - a_read[99][98]))));
		
		/////////////////////////////////////////////////////////////////
		//first few cases
		//////////////////////////////////////////////////////////////////

	 	 //loop along first row
	 	for(int j=1; j<(N-2); ++j)
	 	{
	  		a_write[0][j] = a_read[0][j] - ((((u * delta_t) / (2 * delta_dist))) * (((a_read[1][j] - a_read[99][j])) + ((a_read[0][j+1] - a_read[0][j-1]))));
	 	}

	 	//loop along first column
		for(int i=1; i<(N-2); ++i)
		{
	  		a_write[i][0] = a_read[i][0] - ((((u * delta_t) / (2 * delta_dist))) * (((a_read[i+1][0] - a_read[i-1][0])) + ((a_read[i][1] - a_read[i][99]))));
		}

		/////////////////////////////////////////////////////////////////
		//bulk
		//////////////////////////////////////////////////////////////////

		for(int i=1; i<(N-2); ++i)
	  	{
	  		for(int j=1; j<(N-2); ++j)
	  		{
	  			a_write[i][j] = a_read[i][j] - ((((u * delta_t) / (2 * delta_dist))) * (((a_read[i+1][j] - a_read[i-1][j])) + ((a_read[i][j+1] - a_read[i][j-1]))));
	  		}
	  	}

		/////////////////////////////////////////////////////////////////
		//last few cases
		//////////////////////////////////////////////////////////////////

		//loop along last row
	 	for(int j=1; j<(N-2); ++j)
	 	{
	  		a_write[99][j] = a_read[99][j] - ((((u * delta_t) / (2 * delta_dist))) * (((a_read[0][j] - a_read[98][j])) + ((a_read[99][j+1] - a_read[99][j-1]))));
	 	}

	 	//loop along last column
		for(int i=1; i<(N-2); ++i)
		{
	  		a_write[i][99] = a_read[i][99] - ((((u * delta_t) / (2 * delta_dist))) * (((a_read[i+1][99] - a_read[i-1][99])) + ((a_read[i][0] - a_read[i][98]))));
		}


		/////////////////////////////////////////////////////////////////
		//swap arrays
		//////////////////////////////////////////////////////////////////

        temp = a_read;
        a_read= a_write;
        a_write=temp;

	}


	print_array(a_write,N);




	for(int i = 0; i < N; ++i) 
	{
		delete [] temp[i];
    	delete [] a_read[i];
    	delete [] a_write[i];
	}

	//delete [] dist_x;
	//delete [] dist_y;

	delete [] temp;
    delete [] a_read;
	delete [] a_write;

}


