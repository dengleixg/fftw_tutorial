#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include <fftw3.h>
#include <complex.h>

#include "util.h"

int test_1d_c2c(int n) {

    int dir = FFTW_BACKWARD;
    double _real, _imag, phi;

    _Dcomplex*in = fftw_alloc_complex(n);
    _Dcomplex* ref_out = fftw_alloc_complex(n);
    fftw_complex *fftw_out = fftw_alloc_complex(n);

    fftw_plan p = fftw_plan_dft_1d(n, in, fftw_out, dir, FFTW_ESTIMATE);

    // fill the input array with random data
    fill_random_1d_cplx(n, in);

    // compute the reference output
    for (int k = 0; k < n; ++k) {
        ref_out[k] = _Cbuild(0.0, 0.0);
        for (int j = 0; j < n; ++j) {
            phi = dir * 2.0 * M_PI * j * k / ((double) n);

            _real = creal(in[j]) * cos(phi) - cimag(in[j]) * sin(phi);
            _real += creal(ref_out[k]);
            _imag = creal(in[j]) * sin(phi) + cimag(in[j]) * cos(phi);
            _imag += cimag(ref_out[k]);
            ref_out[k] =  _Cbuild(_real , _imag);
        }
    }

    // compute the DFT of in using FFTW
    fftw_execute(p);

    // compare reference output with FFTW output
    double eps = 1e-12;
    int status = compare_1d_cplx(n, ref_out, fftw_out, eps);

    fftw_destroy_plan(p);
    fftw_free(in);
    fftw_free(ref_out);
    fftw_free(fftw_out);

    return status;
}

int main(int argc, char** argv) {
    int status = 0;
    status += test_1d_c2c(32);
    status += test_1d_c2c(33);
    return status;
}