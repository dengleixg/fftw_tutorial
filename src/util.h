#ifndef UTIL_H
#define UTIL_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include <fftw3.h>
#include <complex.h>

#define M_PI      3.14159265358979323846   // pi
//Macros for creal and imaginary parts
#define REAL 0
#define IMAG 1

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
    double creal, cimag;
    for (int i = 0; i < n; ++i) {
        creal = rand() / ((double) RAND_MAX);
        cimag = rand() / ((double) RAND_MAX);
        arr[i][REAL] = creal;
        arr[i][IMAG] = cimag;
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
    double creal, cimag;
    for (int i = 0; i < rows; ++i) {
        for (int j = 0; j < cols; ++j) {
            creal = rand() / ((double) RAND_MAX);
            cimag = rand() / ((double) RAND_MAX);
            arr[i * cols + j][REAL] = creal;
            arr[i * cols + j][IMAG] = cimag;
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
            printf("  creal ok (delta = %g)\n", delta);
        } else {
            printf("  creal: expected %g, got %g (delta = %g)\n", ref[i], arr[i],
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

int compare_1d_cplx(int n, _Dcomplex* ref, _Dcomplex* arr, double eps) {
    int status = 0;
    double delta_real, delta_imag;

    for (int i = 0; i < n; ++i) {
        printf("compare %d (len = %d)\n", i, n);

        delta_real = creal(arr[i]) - creal(ref[i]);
        delta_imag = cimag(arr[i]) - cimag(ref[i]);

        if (fabs(delta_real) < eps) {
            printf("  creal ok (delta = %g)\n", delta_real);
        } else {
            printf("  creal: expected %g, got %g (delta = %g)\n", creal(ref[i]),
                    creal(arr[i]), delta_real);
            status = 1;
        }

        if (fabs(delta_imag) < eps) {
            printf("  cimag ok (delta = %g)\n", delta_imag);
        } else {
            printf("  cimag: expected %g, got %g (delta = %g)\n", cimag(ref[i]),
                    cimag(arr[i]), delta_imag);
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
                printf("  creal ok (delta = %g)\n", delta);
            } else {
                printf("  creal: expected %g, got %g (delta = %g)\n", ref[index],
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

int compare_2d_cplx(int rows, int cols, _Dcomplex* ref, _Dcomplex* arr,
        double eps) {
    int status = 0, index;
    double delta_real, delta_imag;

    for (int i = 0; i < rows; ++i) {
        for (int j = 0; j < cols; ++j) {
            printf("compare (%d,%d) (shape = %dx%d)\n", i, j, rows, cols);

            index = i * cols + j;

            delta_real = creal(arr[index]) - creal(ref[index]);
            delta_imag = cimag(arr[index]) - cimag(ref[index]);

            if (fabs(delta_real) < eps) {
                printf("  creal ok (delta = %g)\n", delta_real);
            } else {
                printf("  creal: expected %g, got %g (delta = %g)\n",
                        creal(ref[index]), creal(arr[index]), delta_real);
                status = 1;
            }

            if (fabs(delta_imag) < eps) {
                printf("  cimag ok (delta = %g)\n", delta_imag);
            } else {
                printf("  cimag: expected %g, got %g (delta = %g)\n",
                        cimag(ref[index]), cimag(arr[index]), delta_imag);
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
void dft_1d_cplx(int n, _Dcomplex* in, _Dcomplex* out, double a, double b) {
    double phi;
    _Dcomplex ii = _Cbuild(0.0, 1.0);
    for (int k = 0; k < n; ++k) {
        out[k] = _Cbuild(0.0, 0.0);
        for (int j = 0; j < n; ++j) {
            phi = -2.0 * M_PI * (k + b) * (j + a) / ((double)n);
            _Dcomplex tmp = _Cmulcc(in[j], cexp(_Cmulcr(ii, phi)));
            double _real = creal(out[k]) + creal(tmp);
            double _imag = cimag(out[k]) + cimag(tmp);
            out[k] = _Cbuild(_real, _imag);
        }
    }
}



#endif // UTIL_H
