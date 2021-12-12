
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// 2D Viscous Burgers Equation (Parallel)
// da/dt + a*(da/dx) + a*(da/dy) = nu ((d^2a/dx^2) + (d^2a/dy^2))
// Carolyn Wendeln
// 12/04/2021
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include <omp.h>


using std::cout; using std::endl;
using std::to_string;
using namespace std;

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Output Array Functions
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void print_array(double** arr, int N, double* x_dist, double* y_dist) 
{ 
	//cout << " x " << " y " << " a_write " << endl;
	for (int i = 0; i < N; i++)
	{
		for (int j = 0; j < N; j++)
		{
			cout << x_dist[i] << " " << y_dist[j] << " " << arr[i][j] << endl; 
		}
	}
} 

void output_array(double** arr, int N, int n, double* x_dist, double* y_dist)
{
    ofstream file;
    string file_name = "2D_Viscous_Burgers_Time_Step_" + to_string(n) + ".csv";
    file.open (file_name);
    file << "x dist" << " ," << "y dist" << " ,"<< "func value" << endl;
    if (file.is_open())
    {
  		for (int i = 0; i < N; i++)
  		{
  			for (int j = 0; j < N; j++)
  			{
  				file << x_dist[i] << " ," << y_dist[j] << " ," << arr[i][j] << endl;
  			}
  		}
  	file.close();
    }
    else cout << "can not open a file";
}


///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Main Function
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

int main (int argc, char* argv[]) {
	
	if (argc != 4)
    	cout << "Please give step size, time step, and thread count" << endl;
  	
  	else
  	{
  	    
  		///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		//Input Data
		///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

		double dist_start = 0;
		double dist_end = 1;
		int N = atoi(argv[1]);

		double delta_dist = (dist_end-dist_start) / (N-1);

		double t_start = 0;
		double t_end = 0.3;
		int n = atof(argv[2]);
		double delta_t = (t_end-t_start) / (n-1);

		double nu = delta_dist;
				
	    int thrd_cnt = atof(argv[3]);
	    
	    int chunk_size = 10;

// 		cout << "delta x: "<< delta_dist << endl; 	
// 		cout << "delta t: " << delta_t << endl; 
	    
		///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		//Generate Arrays
		///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		  
		double* x_dist = new double[N];
		double* y_dist = new double[N];

		double** a_read = new double*[N];
		double** a_write = new double*[N];

		double* loop_time = new double[n];

		double** temp;

		for(int i = 0; i < N; ++i)
		{
	    	a_read[i] = new double[N];
	    	a_write[i] = new double[N];
		}

		///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		//Initialize Arrays
		///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	  	for(int i=0; i<(N); ++i)
	  	{
	  		x_dist[i] = dist_start + delta_dist * i;
	  		y_dist[i] = dist_start + delta_dist * i;
	  	}


	  	// Initial condition
		for(int i=0; i<(N-1); ++i)
	  	{
	  		for(int j=0; j<(N-1); ++j)
	  		{
	  			a_read[i][j] = exp(-1 * (pow(x_dist[i] - 0.5,2) + pow(y_dist[j] - 0.5,2)) * 50.0);
	  		}
	 	 }


		///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		//Apply Advection Equation
		///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	  	///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	 	//loop through time steps
	 	for(int t=0; t<n; ++t)
	 	{
	 		
			///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
			//four corners
			///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

		  	//i = 0  ; j = 0
		  	a_write[0][0] = a_read[0][0]
		  	+ ( ((nu * delta_t) / pow(delta_dist,2)) * (a_read[1][0] - (2 * a_read[0][0]) + a_read[N-1][0]) )
		  	+ ( ((nu * delta_t) / pow(delta_dist,2)) * (a_read[0][1] - (2 * a_read[0][0]) + a_read[0][N-1]) )
		  	- ( ((a_read[0][0] * delta_t)/(2 * delta_dist)) * (a_read[1][0] - a_read[N-1][0]) )
		  	- ( ((a_read[0][0] * delta_t)/(2 * delta_dist)) * (a_read[0][1] - a_read[0][N-1]) )
		  	;

			//i = 0  ; j = 99
		  	a_write[0][N-1] = a_read[0][N-1]
		  	+ ( ((nu * delta_t) / pow(delta_dist,2)) * (a_read[1][N-1] - (2 * a_read[0][N-1]) + a_read[N-1][N-1]) )
		  	+ ( ((nu * delta_t) / pow(delta_dist,2)) * (a_read[0][0] - (2 * a_read[0][N-1]) + a_read[0][N-2]) )
		  	- ( ((a_read[0][N-1] * delta_t)/(2 * delta_dist)) * (a_read[1][N-1] - a_read[N-1][N-1]) )
		  	- ( ((a_read[0][N-1] * delta_t)/(2 * delta_dist)) * (a_read[0][0] - a_read[0][N-2]) )
		  	;

			//i = 99 ; j = 0
		  	a_write[N-1][0] = a_read[N-1][0]
		  	+ ( ((nu * delta_t) / pow(delta_dist,2)) * (a_read[0][0] - (2 * a_read[N-1][0]) + a_read[N-2][0]) )
		  	+ ( ((nu * delta_t) / pow(delta_dist,2)) * (a_read[N-1][1] - (2 * a_read[N-1][0]) + a_read[N-1][N-1]) )
		  	- ( ((a_read[N-1][0] * delta_t)/(2 * delta_dist)) * (a_read[0][0] - a_read[N-2][0]) )
		  	- ( ((a_read[N-1][0] * delta_t)/(2 * delta_dist)) * (a_read[N-1][1] - a_read[N-1][N-1]) )
		  	;

			//i = 99 ; j = 99
		  	a_write[N-1][N-1] = a_read[N-1][N-1]
		  	+ ( ((nu * delta_t) / pow(delta_dist,2)) * (a_read[0][N-1] - (2 * a_read[N-1][N-1]) + a_read[N-2][N-1]) )
		  	+ ( ((nu * delta_t) / pow(delta_dist,2)) * (a_read[N-1][0] - (2 * a_read[N-1][N-1]) + a_read[N-1][N-2]) )
		  	- ( ((a_read[N-1][N-1] * delta_t)/(2 * delta_dist)) * (a_read[0][N-1] - a_read[N-2][N-1]) )
		  	- ( ((a_read[N-1][N-1] * delta_t)/(2 * delta_dist)) * (a_read[N-1][0] - a_read[N-1][N-2]) )
		  	;
			
			///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
			//bulk
			///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	        
	        double start = omp_get_wtime();
	        
	        #pragma omp parallel for schedule(static) num_threads(thrd_cnt)
	        
			for(int i=1; i<(N-2); ++i)
		  	{
		  		for(int j=1; j<(N-2); ++j)
		  		{
		  			a_write[i][j] = a_read[i][j]
		  			+ ( ((nu * delta_t) / pow(delta_dist,2)) * (a_read[i+1][j] - (2 * a_read[i][j]) + a_read[i-1][j]) )
		  			+ ( ((nu * delta_t) / pow(delta_dist,2)) * (a_read[i][j+1] - (2 * a_read[i][j]) + a_read[i][j-1]) )
		  			- ( ((a_read[i][j] * delta_t)/(2 * delta_dist)) * (a_read[i+1][j] - a_read[i-1][j]) )
		  			- ( ((a_read[i][j] * delta_t)/(2 * delta_dist)) * (a_read[i][j+1] - a_read[i][j-1]) )
		  			;
		  		}
		  	}
		  	
	        double end = omp_get_wtime();
	        double time_taken = end-start;
	       // cout << "Time taken is " << time_taken << " for time step " << t << endl;

	        loop_time[t] = time_taken;
	        

			///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
			//first & last row 
			///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

		 	for(int j=1; j<(N-2); ++j)
		 	{
		  		//loop along first row
		  		a_write[0][j] = a_read[0][j]
		  		+ ( ((nu * delta_t) / pow(delta_dist,2)) * (a_read[1][j] - (2 * a_read[0][j]) + a_read[N-1][j]) )
		  		+ ( ((nu * delta_t) / pow(delta_dist,2)) * (a_read[0][j+1] - (2 * a_read[0][j]) + a_read[0][j-1]) )
		  		- ( ((a_read[0][j] * delta_t)/(2 * delta_dist)) * (a_read[1][j] - a_read[N-1][j]) )
		  		- ( ((a_read[0][j] * delta_t)/(2 * delta_dist)) * (a_read[0][j+1] - a_read[0][j-1]) )
		  		;

		 	
		 		//loop along last row
		  		a_write[N-1][j] = a_read[N-1][j]
		  		+ ( ((nu * delta_t) / pow(delta_dist,2)) * (a_read[0][j] - (2 * a_read[N-1][j]) + a_read[N-2][j]) )
		  		+ ( ((nu * delta_t) / pow(delta_dist,2)) * (a_read[N-1][j+1] - (2 * a_read[N-1][j]) + a_read[N-1][j-1]) )
		  		- ( ((a_read[N-1][j] * delta_t)/(2 * delta_dist)) * (a_read[0][j] - a_read[N-2][j]) )
		  		- ( ((a_read[N-1][j] * delta_t)/(2 * delta_dist)) * (a_read[N-1][j+1] - a_read[N-1][j-1]) )
		  		;
		 	}


			///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
			//first & last column  
			///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

			for(int i=1; i<(N-2); ++i)
			{
				//loop along first column 
				a_write[i][0] = a_read[i][0]
		  		+ ( ((nu * delta_t) / pow(delta_dist,2)) * (a_read[i+1][0] - (2 * a_read[i][0]) + a_read[i-1][0]) )
		  		+ ( ((nu * delta_t) / pow(delta_dist,2)) * (a_read[i][1] - (2 * a_read[i][0]) + a_read[i][N-1]) )
		  		- ( ((a_read[i][0] * delta_t)/(2 * delta_dist)) * (a_read[i+1][0] - a_read[i-1][0]) )
		  		- ( ((a_read[i][0] * delta_t)/(2 * delta_dist)) * (a_read[i][1] - a_read[i][N-1]) )
		  		;
			
		 		//loop along last column 
		  		a_write[i][N-1] = a_read[i][N-1]
		  		+ ( ((nu * delta_t) / pow(delta_dist,2)) * (a_read[i+1][N-1] - (2 * a_read[i][N-1]) + a_read[i-1][N-1]) )
		  		+ ( ((nu * delta_t) / pow(delta_dist,2)) * (a_read[i][0] - (2 * a_read[i][N-1]) + a_read[i][N-2]) )
		  		- ( ((a_read[i][N-1] * delta_t)/(2 * delta_dist)) * (a_read[i+1][N-1] - a_read[i-1][N-1]) )
		  		- ( ((a_read[i][N-1] * delta_t)/(2 * delta_dist)) * (a_read[i][0] - a_read[i][N-2]) )
		  		;
			}


			///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
			//swap arrays & output
			///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	        temp = a_read;
	        a_read = a_write;
	        a_write = temp;


	// 		output_array(a_write, N, t, x_dist, y_dist);
			
	        
		}
		
		///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		//Time to Run Main For Loop
		///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	    double sum = 0;
	    for(int t=0; t<n; ++t)
	    {
	        sum=sum+loop_time[t];
	    }
	    double average_time = sum / n;
	   // cout << "Average Time: " << average_time << endl;
	   // cout << "Total Time: " << sum << endl;
	    cout << sum << endl;


		///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		//We're Done, Release the Memory 
		///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

		for(int i = 0; i < N; ++i) 
		{
	    	delete [] a_read[i];
	    	delete [] a_write[i];
		}

		delete [] x_dist;
		delete [] y_dist;

	    delete [] a_read;
		delete [] a_write;

		delete [] loop_time;

	}

}

	



