#include <fftw3.h>
#include <iostream>
#include <cmath>

using namespace std;

//Macros for real and imaginary parts
#define REAL 0
#define IMAG 1

int tutorial1()
{
	// Define the length of the complex arrays
	const int n = 5;
	// Input array
	fftw_complex x[n];
	// This is equivalent to: double x[n][2];

	// Output array
	fftw_complex y[n];

	// Fill the first array with some data

	for (size_t i = 0; i < n; i++)
	{
		x[i][REAL] = i + 1;
		x[i][IMAG] = 0;
	}

	// Plant the 1D FFT and execute it 
	fftw_plan plan = fftw_plan_dft_1d(n, x, y, FFTW_FORWARD, FFTW_ESTIMATE);
	fftw_execute(plan);

	// Do sine cleaning
	fftw_destroy_plan(plan);
	fftw_cleanup();

	// Display the results
	cout << "FFT = " << endl;
	for (size_t i = 0; i < n; i++)
	{
		if (y[i][IMAG] < 0)
			cout << y[i][REAL] << "-" << abs(y[i][IMAG]) << "i" << endl;
		else
			cout << y[i][REAL] << "+" << abs(y[i][IMAG]) << "i" << endl;
	}

	fftw_free(x);
	fftw_free(y);

	return 0;
}