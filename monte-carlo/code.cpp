#include <iostream>
#include <random>
#include <fstream>
#include <cmath>
#include "int.h"

double f(double x) {
  return  14*sin(x)*x;
}

int main(int argc, char *argv[]) {
  if(argc != 6) {
    std::cout << "usage: ./a.out xmin xmax ymin ymax N" << std::endl;
    return -1;
  } 

  double xmin = atof(argv[1]), \
         xmax = atof(argv[2]), \
         ymin = atof(argv[3]), \
         ymax = atof(argv[4]);

  int N = atoi(argv[5]), Kp = 0, Kn = 0;

  std::random_device rd;
  std::mt19937 gen(rd());
  std::uniform_real_distribution<> xdis(xmin, xmax);
  std::uniform_real_distribution<> ydis(ymin, ymax);

  std::ofstream valid("valid.dat");
  std::ofstream rest("rest.dat");
  std::ofstream graph("graph.dat");

  for(int i = 0; i < N; i++) {
    double xi = xdis(gen), yi = ydis(gen);
    if(f(xi)*yi > 0) {
      if(fabs(f(xi)) - fabs(yi) > 0 && f(xi) > 0) {
        Kp++;
        valid << xi << " " << yi << std::endl;
        graph << xi << " " << f(xi) << std::endl;
        continue;
      } else if (fabs(f(xi)) - fabs(yi) > 0 && f(xi) < 0) {
        Kn++;
        valid << xi << " " << yi << std::endl;
        graph << xi << " " << f(xi) << std::endl;
        continue;
      }
      
    }
    rest << xi << " " << yi << std::endl;
    graph << xi << " " << f(xi) << std::endl;
  }

  double I = (double)(Kp - Kn) / N * (xmax - xmin)*(ymax - ymin);

  double eps = 10e-3;
  double res = 0, err = 0;

  solve(res, err, xmin, xmax, eps, f);


  std::cout << "Integral value using Monte-Carlo method:" << std::endl;
  std::cout << I << std::endl;
  std::cout << "Integral value using 3/8 rule:" << std::endl;
  std::cout << res << " +- " << err << std::endl;
  std::cout << "difference: " << fabs(I - res) << std::endl;

  FILE *gp = popen("gnuplot -persist", "w");
  fprintf(gp, "set grid x y\n");
  fprintf(gp, "show grid\n");
  fprintf(gp, "plot 'valid.dat', 'rest.dat', 'graph.dat'\n");
  fprintf(gp, "pause mouse close\n");
  pclose(gp);

  return 0;
}