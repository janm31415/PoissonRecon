#include "poisson_recon.h"

#include "libpoisson/estimate_normals.h"
#include "libpoisson/poisson_reconstruction_screened.h"

#include "test_assert.h"

#include <cmath>
#include <iostream>

namespace {

  void test_poisson_1() {
    std::vector<std::array<double, 3>> points;
    
    const double pi = 3.1415926535897;
    
    double r = 1.0;
    int max_i = 10;
    int max_j = 20;
    for (int i = 0; i < max_i; ++i) {
      double theta = (double)i/(double)max_i * pi;
      for (int j = 0; j < max_j; ++j) {
        double phi = (double)j/(double)max_j * 2.0 * pi;
        double x = r*std::cos(phi)*std::sin(theta);
        double y = r*std::sin(phi)*std::sin(theta);
        double z = r*std::cos(theta);
        points.push_back({{x, y, z}});
      }
    }
    
    
    auto normals = libpoisson::estimate_normals((const double*)points.data(), (uint32_t)points.size());
    
    std::vector<std::array<double, 3>> vertices;
    std::vector<std::array<uint32_t, 3>> triangles;
    
    libpoisson::poisson_reconstruction_screened_parameters pars;
    pars.output_stream = &std::cout;
    
    libpoisson::poisson_reconstruction_screened(vertices, triangles, (const double*)points.data(), (const double*)normals.data(), (uint32_t)points.size(), pars);
    
    TEST_ASSERT(vertices.size() > 100);
    TEST_ASSERT(triangles.size() > 100);
    for (const auto& v : vertices) {
      double length = sqrt(v[0]*v[0]+v[1]*v[1]+v[2]*v[2]);
      TEST_EQ_CLOSE(length, r, 1e-1);
    }
  }

}

void run_all_poisson_reconstruction_tests() {
  test_poisson_1();
}
