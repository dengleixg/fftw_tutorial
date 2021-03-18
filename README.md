# FFTW Tutorial

This is a basic C project (Makefile, but also for the Eclipse IDE) I use for exploring FFTW.

One- and two-dimensional discrete Fourier transforms (DFTs) of random data are computed using both FFTW and straight-forward naive algorithms
in order to illustrate explicitly what kind of symmetries and scaling properties FFTW implies in its inputs and outputs.

## One-Dimensional Examples
This tutorial starts by computing one-dimensional (1D) DFTs of random input data.

### 1D complex-to-complex
The first example is basically a self-contained version of the [corresponding example in the FFTW manual](http://fftw.org/fftw3_doc/Complex-One_002dDimensional-DFTs.html#Complex-One_002dDimensional-DFTs). 

We want to compute the complex-valued one-dimensional DFT here, which is specified in [section 4.8.1 of the FFTW reference manual](http://fftw.org/fftw3_doc/The-1d-Discrete-Fourier-Transform-_0028DFT_0029.html#The-1d-Discrete-Fourier-Transform-_0028DFT_0029).

![complex-valued DFT](eqn/complex_dft.png)

The sign in the exponent of the basis function specifies the direction in which the Fourier transform is to be computed:
`-1` indicates a "forward" transform and `+1` indicates a backward transform. These values are available via the `FFTW_FORWARD` and `FFTW_BACKWARD` preprocessor macros.

In order to compute the DFT, complex-valued products of the following form need to be evaluated:

![complex product](eqn/c2c_product.png)

Eulers formula comes in handy now (where *i* is the imaginary unit with *i*^2=-1):

![Eulers formula](eqn/euler.png)

The angle argument `phi` can be identified in above formulas for the DFT:

![phi in DFT](eqn/phi.png)

Now the complex-valued product can be computed using only real-valued variables:

![complex product using real numbers](eqn/c2c_product_real.png)

FFTW implements all this math internally and the explicit formulation was only used to build a basis for the computations to follow.
Below is the example code showing how to compute the 1d `c2c` DFT using both FFTW and a manual implementation.
The size of the DFT is specified via the variable `n` and the direction (forward or backward) is specified via the variable `dir`.
Complex-valued arrays `in`, `ref_out` and `fftw_out` are allocated to hold the input (*X_k*) and outputs (*Y_k*) of the DFT.
A plan for the corresponding DFT is created using `fftw_plan_dft_1d`.
Only after this, the input data is written into the `in` array.
The reference output is computed now (before calling `fftw_execute`), since in some later examples (most notably multi-dimensional `c2r` transforms),
FFTW overwrites the input data and for good practise, we keep this in mind already now.
Note that this is an out-of-place transform, since `in` and `fftw_out` are allocated to be separate arrays.
Next the FFTW transform can be executed via `fftw_execute`.
This fills the corresponding output array `fftw_out`, which is subsequently compared against the reference output in `ref_out`.
A conservative tolerance of `1e-12` is specified to make the example work also for weird input data (as generated by the PRNG).
The actual deviations are usually much smaller and can be observed in the screen output from `compare_1d_cplx`.
Finally, the `fftw_plan` is destroyed and the memory is released using `fftw_free` (which has to be used if the array was allocated using `fftw_alloc_*`).

```C
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include <complex.h>
#include <fftw3.h>

#include "util.h"

void test_1d_c2c() {
    int n = 32;
    int dir = FFTW_BACKWARD;
    double real, imag, phi;

    fftw_complex *in = fftw_alloc_complex(n);
    fftw_complex *ref_out = fftw_alloc_complex(n);
    fftw_complex *fftw_out = fftw_alloc_complex(n);

    fftw_plan p = fftw_plan_dft_1d(n, in, fftw_out, dir, FFTW_ESTIMATE);

    // fill the input array with random data
    fill_random_1d_cplx(n, in);

    // compute the reference output
    for (int k = 0; k < n; ++k) {
        ref_out[k] = 0.0;
        for (int j = 0; j < n; ++j) {
            phi = dir * 2.0 * M_PI * j * k / ((double) n);

            real = creal(in[j]) * cos(phi) - cimag(in[j]) * sin(phi);
            imag = creal(in[j]) * sin(phi) + cimag(in[j]) * cos(phi);
            ref_out[k] += real + I * imag;
        }
    }

    // compute the DFT of in using FFTW
    fftw_execute(p);

    // compare reference output with FFTW output
    double eps = 1e-12;
    compare_1d_cplx(n, ref_out, fftw_out, eps);

    fftw_destroy_plan(p);
    fftw_free(in);
    fftw_free(ref_out);
    fftw_free(fftw_out);
}

int main(int argc, char** argv) {
    test_1d_c2c();
    return 0;
}
```

The code is available in the file [`src/test_1d_c2c.c`](src/test_1d_c2c.c).

### 1D real-to-complex and complex-to-real
The next two examples deal with DFTs of purely real data (`r2c`) and DFTs which *produce* purely real data (`c2r`).
These are covered in the [official FFTW tutorial](http://fftw.org/fftw3_doc/One_002dDimensional-DFTs-of-Real-Data.html#One_002dDimensional-DFTs-of-Real-Data)
as well as in the [FFTW reference manual](http://fftw.org/fftw3_doc/The-1d-Real_002ddata-DFT.html#The-1d-Real_002ddata-DFT).

In case either the input array or the output array are constrained to be purely real, the corresponding complex-valued output or input array
features Hermitian symmetry (where the `n`-periodicity has been included as well):

![Hermitian symmetry](eqn/hermitian.png)

For the case of a backward DFT of real-valued *X_j* and complex-valued *Y_k* with Hermitian symmetry,
the Fourier sum is written out excplicitly as follows:

![Hermitian symmetry in Fourier sum](eqn/hermitian_sum.png)

The figure below illustrates the structure of the complex-valued Fourier space arrays
occuring in the DFT for both even-valued (`n=6`) and odd-valued (`n=7`) sizes of the DFT.

![array structure](img/array_structures.png)

The size required to contain all information required for the transform from or to a real-valued array
is contained in the first `n/2+1` (division by 2 rounded down) entries of the complex array, indicated by the red bars in above figure.
For both even and odd values of `n`, Hermitian symmetry implies *Y_0* = *Y\*_0* and thus *Y_0* is always real.
In the case of even `n`, we can intuitively observe that furthermore *Y_m* = *Y\*_-m* where *m*=`n/2+1` (the element at the Nyquist frequency),
and thus also the last element of the complex-valued Fourier-space array is purely real.

The DFT is formulated to include all elements of the Fourier-space array.
For odd `n`, all components of the Fourier-space array except the DC element at *k*=0 have to be weighted with a factor of 2.
For even `n`, all components of the Fourier-space array except the DC element at *k*=0 and the Nyquist element at *k*=`n/2` have to be weighted with a factor of 2. The elements that need to be weighted by a factor of 2 are highlighted by solid blue lines in above illustration.
The redundant elements that are not explicitly needed are indicated by dashed blue lines.

The transformation from complex Fourier space to real-valued real space is demonstrated in [`src/test_1d_c2r.c`](src/test_1d_c2r.c).
The relevant portion of the source code is here:

```C
int nCplx = n / 2 + 1;
for (int k = 0; k < n; ++k) {

    // start with DC component, which is purely real due to Hermitian symmetry
    ref_out[k] = creal(in[0]);

    int loopEnd = nCplx;

    // special case for even n
    if (n % 2 == 0) {
        // Nyquist element is purely real as well
        phi = 2.0 * M_PI * (nCplx - 1) * k / ((double) n);
        ref_out[k] += creal(in[nCplx - 1]) * cos(phi);

        loopEnd = nCplx-1;
    }

    // middle elements are handled the same for even and odd n
    for (int j = 1; j < loopEnd; ++j) {
        phi = 2.0 * M_PI * j * k / ((double) n);

        real = creal(in[j]) * cos(phi) - cimag(in[j]) * sin(phi);
        ref_out[k] += 2.0 * real;
    }
}
```

Note that the integer division used to compute `nCplx` (correctly) rounds down.
The rest is a relatively straight-forward implementation of above verbose algorithm.
The DC component is always taken to be real.
Depending on whether `n` is even or odd, the number of elements to take into account with both real and imaginary component (`loopEnd`) is adjusted.
The (purely real) Nyquist element at `n/2` is added separately if `n` is even.

## Allocation of arrays
Throughout this example collection, the proposed convenience wrapper functions provided by FFTW for allocating real- and complex-valued arrays are used:
```C
int n = 32;
int nOut = n/2+1;
double *in = fftw_alloc_real(n);
fftw_complex *out = fftw_alloc_complex(nOut);
```
where `N` is the real-space size of the DFT and `outN` is the number of Fourier coefficients resulting from a `r2c` DFT.
The corresponding "raw" allocation code would look like this:
```C
double *in = (double*) fftw_malloc(sizeof(double) * n);
fftw_complex *out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * nOut);
```
Note that above code is equivalent to the standard C way of allocating memory using `malloc`:
```C
double *in = (double*) malloc(sizeof(double) * n);
fftw_complex *out = (fftw_complex*) malloc(sizeof(fftw_complex) * nOut);
```
except that the FFTW routines ensure proper memory alignment for exploiting SIMD instructions of modern CPUs.

## Utility functions
In order to keep the examples short, a separate header file [`src/util.h`](src/util.h) is provided.
It contains methods to operate on one- and two-dimensional arrays (the latter stored in row-major order)
of real (`double`) and complex (`fftw_complex`) numbers.
The following operations are supported:
 * fill with random numbers between 0 and 1: e.g. `fill_random_1d_cplx`
 * element-wise check for approximate equality: e.g. `compare_1d_real`
 * write into a text file.

