#include <valarray>
#include <cmath>
#include <iostream>
#include <complex>
#include <fstream>
#include <functional>
#include <unistd.h>

typedef std::valarray<std::complex<double>> carray;

carray fft(carray &x) {
  int N = x.size();
  if(N <= 1) return x;

  carray X(N);
  carray even = x[std::slice(0, N/2, 2)];
  carray odd  = x[std::slice(1, N/2, 2)];

  even = fft(even);
  odd = fft(odd);

  for(int k = 0; k < N/2; k++) {
    std::complex<double> q = std::polar(1., -2*M_PI*k/N)*odd[k];
    X[k] = even[k] + q;
    X[k + N/2] = even[k] - q;
  }
  return X;
}

carray ifft(carray &X) {
  carray XC = X.apply(std::conj);
  carray x = fft(XC);
  return x.apply(std::conj)/x.size();
}

double f(double t) {
  double cf = 2*M_PI;
  return cos(cf*t) + 2*cos(cf*2*t) + 16*sin(cf*3*t) + 11;
}

carray sample_func(double x0, double x, int N, std::function<double(double)> func) {
  carray c(N);
  double h = (x - x0)/N;
  for(int i = 0; i < N; i++) {
    c[i] = std::complex<double>(func(i*h), 0.);
  }
  return c;
}

carray desample_func(double x0, double x, carray &data) {
  int N = data.size();
  carray c(N);
  double h = (x - x0)/N;
  for(int i = 0; i < N; i++) {
    c[i] = std::complex<double>(h*i, data[i].real());
  }
  return c;
}

carray freq_representation(double x0, double x, carray &ffdata) {
  int N = ffdata.size();
  carray c(N);
  double fres = 1./(x - x0);
  for(int i = 0; i < N/2; i++) {
    c[i] = std::complex<double>((i)*fres, std::abs(ffdata[i]));
  }
  for(int i = N/2; i < N; i++) {
    c[i] = std::complex<double>((i - N)*fres, std::abs(ffdata[i]));
  }
  return c;
}

carray gamp_filt(carray &ffdata, double amp) {
  for(int i = 0; i < ffdata.size(); i++) {
    if(std::abs(ffdata[i]) < amp) {
      //ffdata[i] = std::complex<double>(ffdata[i].real(), 0.);
      ffdata[i] = 0;
    }
  }
  return ffdata;
}


int main(int argc, char *argv[]) {
  double xl = atof(argv[1]), \
         xr = atof(argv[2]);
  int N = (int)pow(2, atoi(argv[3]));

  carray data = sample_func(xl, xr, N, f);
  carray ft = fft(data);
  carray ft_hz = freq_representation(xl, xr, ft);

  std::ofstream ofi("fft.dat");
  for(int i = 0; i < ft.size(); i++) {
    ofi << ft_hz[i].real() << " " << ft_hz[i].imag() << std::endl;
  }

  //ft = gamp_filt(ft, 100000.);
  //ft_hz = freq_representation(xl, xr, ft);
  carray ret = ifft(ft);
  ret = desample_func(xl, xr, ret);

  std::ofstream of2("rfft.dat");
  std::ofstream of3("fftt.dat");
  for(int i = 0; i < ft.size(); i++) {
    of2 << ret[i].real() << " " << ret[i].imag() << " " << f(ret[i].real()) << std::endl;
    of3 << ft_hz[i].real() << " " << ft_hz[i].imag() << std::endl;
  }

  (void)execl("/usr/bin/gnuplot", "/usr/bin/gnuplot", "plot.gp", (char *)0);

  return 0;
}