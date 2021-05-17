#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <random>
#include <algorithm>
#include <vector>
#include <cmath>
#include <cfloat>
#include <unistd.h>

struct point {
  double x, y;
  int cluster;
  double mdist = DBL_MAX;

  point() {
    x = 0.; y = 0.;
    cluster = -1;
  }

  point(double _x, double _y) {
    x = _x; y = _y;
    cluster = -1;
  }

  double dist(point b) {
    return sqrt((b.x - x)*(b.x - x) + (b.y - y)*(b.y - y));
  }
};

std::vector<point> parse_csv(int xcol, int ycol, char *name) {
  std::vector<point> pts;
  std::string line;
  std::ifstream file(name);

  int chflag = 0;
  if(ycol < xcol) {
    std::swap(xcol, ycol);
    chflag++;
  }

  while(getline(file, line)) {
    std::stringstream lstr(line);
    std::string tmp;
    double x, y;

    for(int i = 0; i <= xcol; i++) {
      getline(lstr, tmp, ',');
    }
    x = std::stod(tmp);
    for(int i = xcol + 1; i <= ycol; i++) {
      getline(lstr, tmp, ',');
    }
    y = std::stod(tmp);

    if(chflag) std::swap(x, y);
    pts.push_back(point(x, y));
  }

  return pts;  
}

void kmeans(std::vector<point> &pts, int iters, int k) {
  //assign centroids randomly
  std::vector<point> centroids(k);
  std::random_device rd;
  std::mt19937 gen(rd());
  std::uniform_int_distribution<> dis(0, pts.size() - 1);
  for(int i = 0; i < k; i++) {
    centroids[i] = pts[dis(gen)];
  }


  //initial clustering by distance
  for(auto c = centroids.begin(); c != centroids.end(); c++) {
    for(auto p = pts.begin(); p != pts.end(); p++) {
      double d = c->dist(*p);
      if(d < p->mdist) {
        p->mdist = d;
        p->cluster = c - centroids.begin();
      }
    }
  }

  //main loop
  for(int t = 0; t < iters; t++) {
    //recalculate centroids
    std::vector<int> npts(k, 0);
    std::vector<double> sumx(k, 0.), sumy(k, 0.);

    for(auto p = pts.begin(); p != pts.end(); p++) {
      npts[p->cluster]++;
      sumx[p->cluster] += p->x;
      sumy[p->cluster] += p->y;

      p->mdist = DBL_MAX;
    }

    for(auto c = centroids.begin(); c != centroids.end(); c++) {
      c->x = sumx[c - centroids.begin()] / npts[c - centroids.begin()];
      c->y = sumy[c - centroids.begin()] / npts[c - centroids.begin()];
    }

    //rearrange clusters
    for(auto c = centroids.begin(); c != centroids.end(); c++) {
      for(auto p = pts.begin(); p != pts.end(); p++) {
        double d = c->dist(*p);
        if(d < p->mdist) {
          p->mdist = d;
          p->cluster = c - centroids.begin();
        }
      }
    }
  }

  //write resulting clusters to a file
  std::sort(pts.begin(), pts.end(), [](const point &a, const point &b) {
    return a.cluster < b.cluster;
  });
  std::ofstream out("clusters.dat");

  auto p = pts.begin();
  for(int i = 0; i < k; i++) {
    out << "\"cluster " << i << "\""<< std::endl;
    while(p != pts.end() && p->cluster == i) {
      out << p->x << " " << p->y << std::endl;
      p++;
    }
    out << std::endl << std::endl;
  }

}

int main(int argc, char *argv[]) {
  char *fname = argv[1];
  int xcol = atoi(argv[2]) - 1;
  int ycol = atoi(argv[3]) - 1;
  int k = atoi(argv[4]);
  int iters = atoi(argv[4]);

  std::vector<point> pts = parse_csv(xcol, ycol, fname);

  kmeans(pts, iters, k);

  (void)execl("/usr/bin/gnuplot", "/usr/bin/gnuplot", "plot.gp", (char *)0);

  return 0;
}