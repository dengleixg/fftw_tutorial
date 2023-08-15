
#include <fftw3.h>
#include <complex>
#include "util.h"


int test_1d_c2c(int n) {

    int dir = FFTW_BACKWARD;
    double _real, _imag, _phi;

    complex<double>* in = (complex<double>*)fftw_alloc_complex(sizeof(fftw_complex) * n);
    complex<double>* ref_out = (complex<double>*)fftw_alloc_complex(sizeof(fftw_complex) * n);
    fftw_complex *fftw_out = (fftw_complex*)fftw_alloc_complex(sizeof(fftw_complex) * n);

    fftw_plan p = fftw_plan_dft_1d(n, 
        reinterpret_cast<fftw_complex*>(in), fftw_out, dir, FFTW_ESTIMATE);

    // fill the input array with random data
    fill_random_1d_cplx(n, 
        reinterpret_cast<fftw_complex*>(in));

    // compute the reference output
    for (int k = 0; k < n; ++k) {
        ref_out[k] = 0.0;
        ref_out[k] = 0.0;
        for (int j = 0; j < n; ++j) {
            _phi = dir * 2.0 * M_PI * j * k / ((double) n);

            _real = real(in[j]) * cos(_phi) - imag(in[j]) * sin(_phi);
            _imag = real(in[j]) * sin(_phi) + imag(in[j]) * cos(_phi);
            ref_out[k] += _real + 1i * _imag;
        }
    }

    // compute the DFT of in using FFTW
    fftw_execute(p);

    // compare reference output with FFTW output
    double eps = 1e-12;
    int status = compare_1d_cplx(n, ref_out, (complex<double>*)fftw_out, eps);

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
