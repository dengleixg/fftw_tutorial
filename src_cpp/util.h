#ifndef UTIL_H
#define UTIL_H

#include <iostream>
#include <string>
#include <math.h>

#include <complex>
#include <fftw3.h>

#define M_PI      3.14159265358979323846   // pi
//Macros for real and imaginary parts
#define REAL 0
#define IMAG 1

using namespace std;

void fill_zero_1d_real(int n, double *arr) {
    for (int i = 0; i < n; ++i) {
            arr[i] = 0.0;
        }
}

void fill_zero_1d_cplx(int n, fftw_complex *arr) {
    for (int i = 0; i < n; ++i) {
        arr[i][IMAG] = 0.0;
    }
}

void fill_zero_2d_real(int rows, int cols, double *arr) {
    for (int i = 0; i < rows; ++i) {
        for (int j = 0; j < cols; ++j) {
            arr[i * cols + j] = 0.0;
        }
    }
}

void fill_zero_2d_cplx(int rows, int cols, fftw_complex *arr) {
    for (int i = 0; i < rows; ++i) {
        for (int j = 0; j < cols; ++j) {
            arr[i * cols + j][IMAG] = 0.0;
        }
    }
}

void fill_random_1d_real(int n, double *arr) {
    for (int i = 0; i < n; ++i) {
        arr[i] = rand() / ((double) RAND_MAX);
    }
}

void fill_random_1d_cplx(int n, fftw_complex *arr) {
    double real, imag;
    for (int i = 0; i < n; ++i) {
        real = rand() / ((double) RAND_MAX);
        imag = rand() / ((double) RAND_MAX);
        arr[i][REAL] = real;
        arr[i][IMAG] = imag;
    }
}

void fill_random_2d_real(int rows, int cols, double *arr) {
    for (int i = 0; i < rows; ++i) {
        for (int j = 0; j < cols; ++j) {
            arr[i * cols + j] = rand() / ((double) RAND_MAX);
        }
    }
}

void fill_random_2d_cplx(int rows, int cols, fftw_complex *arr) {
    double real, imag;
    for (int i = 0; i < rows; ++i) {
        for (int j = 0; j < cols; ++j) {
            real = rand() / ((double) RAND_MAX);
            imag = rand() / ((double) RAND_MAX);
            arr[i * cols + j][REAL] = real;
            arr[i * cols + j][IMAG] = imag;
        }
    }
}

int compare_1d_real(int n, double *ref, double *arr, double eps) {
    int status = 0;
    double delta;

    for (int i = 0; i < n; ++i) {
        printf("compare %d (len = %d)\n", i, n);

        delta = arr[i] - ref[i];

        if (fabs(delta) < eps) {
            printf("  real ok (delta = %g)\n", delta);
        } else {
            printf("  real: expected %g, got %g (delta = %g)\n", ref[i], arr[i],
                    delta);
            status = 1;
        }
    }

    if (status == 0) {
        printf("=> all ok\n");
    } else {
        printf("=> errors\n");
    }

    return status;
}

int compare_1d_cplx(int n, complex<double>* ref, complex<double>* arr, double eps) {
    int status = 0;
    double delta_real, delta_imag;

    for (int i = 0; i < n; ++i) {
        printf("compare %d (len = %d)\n", i, n);

        delta_real = real(arr[i]) - real(ref[i]);
        delta_imag = imag(arr[i]) - imag(ref[i]);

        if (fabs(delta_real) < eps) {
            printf("  real ok (delta = %g)\n", delta_real);
        } else {
            printf("  real: expected %g, got %g (delta = %g)\n", real(ref[i]),
                    real(arr[i]), delta_real);
            status = 1;
        }

        if (fabs(delta_imag) < eps) {
            printf("  imag ok (delta = %g)\n", delta_imag);
        } else {
            printf("  imag: expected %g, got %g (delta = %g)\n", imag(ref[i]),
                    imag(arr[i]), delta_imag);
            status = 1;
        }
    }

    if (status == 0) {
        printf("=> all ok\n");
    } else {
        printf("=> errors\n");
    }

    return status;
}

int compare_2d_real(int rows, int cols, double *ref, double *arr, double eps) {
    int status = 0, index;
    double delta;

    for (int i = 0; i < rows; ++i) {
        for (int j = 0; j < cols; ++j) {
            printf("compare (%d,%d) (shape = %dx%d)\n", i, j, rows, cols);

            index = i * cols + j;

            delta = arr[index] - ref[index];

            if (fabs(delta) < eps) {
                printf("  real ok (delta = %g)\n", delta);
            } else {
                printf("  real: expected %g, got %g (delta = %g)\n", ref[index],
                        arr[index], delta);
                status = 1;
            }
        }
    }

    if (status == 0) {
        printf("=> all ok\n");
    } else {
        printf("=> errors\n");
    }

    return status;
}

int compare_2d_cplx(int rows, int cols, complex<double>* ref, complex<double>* arr,
        double eps) {
    int status = 0, index;
    double delta_real, delta_imag;

    for (int i = 0; i < rows; ++i) {
        for (int j = 0; j < cols; ++j) {
            printf("compare (%d,%d) (shape = %dx%d)\n", i, j, rows, cols);

            index = i * cols + j;

            delta_real = real(arr[index]) - real(ref[index]);
            delta_imag = imag(arr[index]) - imag(ref[index]);

            if (fabs(delta_real) < eps) {
                printf("  real ok (delta = %g)\n", delta_real);
            } else {
                printf("  real: expected %g, got %g (delta = %g)\n",
                        real(ref[index]), real(arr[index]), delta_real);
                status = 1;
            }

            if (fabs(delta_imag) < eps) {
                printf("  imag ok (delta = %g)\n", delta_imag);
            } else {
                printf("  imag: expected %g, got %g (delta = %g)\n",
                        imag(ref[index]), imag(arr[index]), delta_imag);
                status = 1;
            }
        }
    }

    if (status == 0) {
        printf("=> all ok\n");
    } else {
        printf("=> errors\n");
    }

    return status;
}

// manual implementation of a general DFT (shift and non-linear phase)
// formula from https://en.wikipedia.org/wiki/Discrete_Fourier_transform#Generalized_DFT_(shifted_and_non-linear_phase) (accessed 2021-03-23)
// a: shift of input
// b: shift of output
void dft_1d_cplx(int n, complex<double>* in, complex<double>* out, double a, double b) {
    complex<double> phi;
    for (int k = 0; k < n; ++k) {
        out[k] = 0.0;
        for (int j = 0; j < n; ++j) {
            phi = -2.0i * M_PI * (k + b) * (j + a) / ((double)n);
            out[k] += in[j] * exp(phi);
        }
    }
}



#endif // UTIL_H
