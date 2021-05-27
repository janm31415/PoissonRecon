#include "estimate_normals.h"

#include "libpoisson/estimate_normals.h"

#include "test_assert.h"

namespace {
void test_estimate_normals_1() {
  std::vector<std::array<float, 3>> pts;
  
  for (int i = 0; i < 20; ++i)
    pts.push_back({{(float)(i%5), (float)(i/5), 0.f}});
    
  auto normals = libpoisson::estimate_normals((const float*)pts.data(), (uint32_t)pts.size());
  
  for (const auto& n : normals) {
    TEST_EQ(0.f, n[0]);
    TEST_EQ(0.f, n[1]);
    TEST_EQ(1.f, n[2]);
  }
  
  }
  
}


void run_all_estimate_normals_tests() {
  test_estimate_normals_1();
}
