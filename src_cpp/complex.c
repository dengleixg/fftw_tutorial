#include <complex.h>
#include <stdio.h>
#include <math.h>

int main(void)
{
    // initilize float complex
    _Fcomplex zf = _FCbuild(1.3f, 4.9f);
    printf("z = %.1f% + .1fi\n", crealf(zf), cimagf(zf));

    // initilize double complex
    double real = 1.3, imag = 4.9;
    _Dcomplex z = { real, imag };
    printf("z = %.1f% + .1fi\n", creal(z), cimag(z));
 
    // complex conjugate
    _Dcomplex conjz = conj(z);
    _Fcomplex conjzf = conjf(zf);
    printf("conj(z) = %.1f% + .1fi\n",creal(conjz),cimag(conjz));
    printf("conj(z) = %.1f% + .1fi\n", crealf(conjzf), cimagf(conjzf));

    // absolute value
    printf("Absolute value z= %.1f\n", cabs(z));
    printf("Absolute value zf= %.1f\n", cabsf(zf));

    // phase angle of a complex number
    printf("Phase Angle z = %.1f radians\n", carg(z));
    printf("Phase Angle zf = %.1f radians\n", cargf(zf));

    // natural exponential of a complex number
    printf("exp z = %.1f% + .1fi\n", creal(cexp(z)), cimag(cexp(z)));
    printf("exp zf = %.1f% + .1fi \n", crealf(cexpf(zf)), cimagf(cexpf(zf)));

    // complex unit
    printf("I = %.1f% + .1fi\n", crealf(I), cimagf(I));

	return 0;
}