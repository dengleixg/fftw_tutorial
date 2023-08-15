
#include <complex>
#include <fftw3.h>
#include "util.h"

int test_1d_r2c(int n) {
    double _real, _imag, phi;

    int nCplx = n / 2 + 1;

    printf("n=%d nCplx=%d\n", n, nCplx);

    double *in = fftw_alloc_real(n);
    complex<double>* ref_out = (complex<double>*)fftw_alloc_complex(nCplx);
    fftw_complex *fftw_out = fftw_alloc_complex(nCplx);

    fftw_plan p = fftw_plan_dft_r2c_1d(n, in, fftw_out, FFTW_ESTIMATE);

    // fill the input array with random data
    fill_random_1d_real(n, in);

    // compute the reference output
    for (int k = 0; k < nCplx; ++k) {

        // DC component is always real
        ref_out[k] = in[0];

        for (int j = 1; j < n; ++j) {
            phi = -2.0 * M_PI * j * k / ((double) n);

            _real = in[j] * cos(phi);
            _imag = in[j] * sin(phi);
            ref_out[k] += _real + 1i * _imag;
        }
    }

    // compute the DFT of in using FFTW
    fftw_execute(p);

    // compare reference output with FFTW output
    double eps = 1e-12;
    int status = compare_1d_cplx(nCplx, ref_out, (complex<double>*)fftw_out, eps);

    fftw_destroy_plan(p);
    fftw_free(in);
    fftw_free(ref_out);
    fftw_free(fftw_out);

    return status;
}

int main(int argc, char** argv) {
    int status = 0;
    status += test_1d_r2c(32);
    status += test_1d_r2c(33);
    return status;
}
