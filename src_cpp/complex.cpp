#include <iostream>
#include <complex>

using namespace std;

int main()
{
	std::complex<double> cplx(3., 4.);

	// real()
	cout << "real part:" << real(cplx) << endl;
	// imag()
	cout << "imaginary part:" << imag(cplx) << endl;

	// abs()
	cout << "absolute value of " << cplx << " is: ";
	cout << abs(cplx) << endl;

	// arg()
	cout << "argument of " << cplx << " is: ";
	cout << arg(cplx) << endl;

	// use of polar()
	cout << "The complex whose magnitude is " << 2.0;
	cout << " and phase angle is " << 0.5;
	cout << " is " << polar(2.0, 0.5) << endl;

	// use of norm()
	cout << "The norm of " << cplx << " is "
		<< norm(cplx) << endl;

	// use of conj()
	std::complex<double> mycomplex(10.0, 2.0);
	cout << "The conjugate of " << mycomplex << " is: ";
	cout << conj(mycomplex) << endl;

	// proj()
	std::complex<double> c1(1, 2);
	cout << "proj" << c1 << " = " << proj(c1) << endl;
	std::complex<double> c2(INFINITY, -1);
	cout << "proj" << c2 << " = " << proj(c2) << endl;
	std::complex<double> c3(0, -INFINITY);
	cout << "proj" << c3 << " = " << proj(c3) << endl;

	// use of sqrt()
	cout << "Square root of -4 is "
		<< sqrt(std::complex<double>(-4, 0)) << endl
		<< "Square root of (-4,-0), the other side of the cut, is "
		<< sqrt(std::complex<double>(-4, -0.0)) << endl;

	// use of log()
	cout << "The log of " << mycomplex << " is "
		<< log(mycomplex) << endl;

	// initializing the complex: (-1.0+0.0i)
	mycomplex = complex<double>(0.0, 1.0);

	// use of cos() cos z = (e^(iz) + e^(-iz))/2
	cout << "The cos of " << mycomplex << " is "
		<< cos(mycomplex) << endl;

	// use of sin() sin z = (e^(iz) - e^(-iz))/2i
	cout << "The sin of " << mycomplex << " is "
		<< sin(mycomplex) << endl;

	// use of tan() tan z = i(e^(-iz) - e^(iz)) / (e^(-iz) + e^iz)
	cout << "The tan of " << mycomplex << " is "
		<< tan(mycomplex) << endl;

	// behaves like real cosh, sinh, tanh along the real line;
	// z = a + 0i
	//  cosh(z) = (e ^ z + e ^ (-z)) / 2
	//	sinh(z) = (e ^ z - e ^ (-z)) / 2.
	//	tanh(z) = (e ^ (2z) - 1) / (e ^ (2z) + 1)
	complex<double> z(1, 0);
	cout << "cosh" << z << " = " << cosh(z)
		<< " (cosh(1) = " << cosh(1) << ")" << endl;
	cout << "sinh" << z << " = " << sinh(z)
		<< " (sinh(1) = " << sinh(1) << ")" << endl;
	cout << "tanh" << z << " = " << tanh(z)
		<< " (tanh(1) = " << tanh(1) << ")" << endl;

	// behaves like real cosine,sine,tangent along the imaginary line; z2=0+1i
	complex<double> z2(0, 1);
	cout << "cosh" << z2 << " = " << cosh(z2)
		<< " ( cos(1) = " << cos(1) << ")" << endl;
	cout << "sinh" << z2 << " = " << sinh(z2)
		<< " ( sin(1) = " << sin(1) << ")" << endl;
	cout << "tanh" << z2 << " = " << tanh(z2)
		<< " ( tan(1) = " << tan(1) << ")" << endl;

	std::complex<double> z1 = 1i * 1i; // imaginary unit squared
	std::cout << "i * i = " << z1 << '\n';

	const double PI = std::acos(-1); // or std::numbers::pi in C++20
	std::complex<double> z3 = std::exp(1i * PI); // Euler's formula
	std::cout << "exp(i * pi) = " << z3 << '\n';

	std::complex<double> z4 = 1.0 + 2i, z5 = 1.0 - 2i; // conjugates
	std::cout << "(1 + 2i) * (1 - 2i) = " << z4 * z5 << '\n';

	return 0;
}