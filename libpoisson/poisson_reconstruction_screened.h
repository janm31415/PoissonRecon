#pragma once

#include "libpoisson_api.h"

#include <stdint.h>
#include <array>
#include <vector>

namespace libpoisson
{

struct poisson_reconstruction_screened_parameters
  {
  poisson_reconstruction_screened_parameters()
    {
    depth = 8;
    octree_depth = 5;
    samples_per_node = 1.5;
    solver_divide = 0;
    output_stream = nullptr;
    }

  int depth;
  int octree_depth;
  double samples_per_node;
  int solver_divide;
  void* output_stream;
  };

LIBPOISSON_API void poisson_reconstruction_screened(
  std::vector<std::array<float, 3>>& vertices,
  std::vector<std::array<uint32_t, 3>>& triangles,
  const float* pts3d,
  const float* normals3d,
  uint32_t number_of_points,
  const poisson_reconstruction_screened_parameters& par);
  
LIBPOISSON_API void poisson_reconstruction_screened(
  std::vector<std::array<float, 3>>& vertices,
  std::vector<std::array<uint32_t, 3>>& triangles,
  std::vector<uint32_t>& vertex_colors,
  const float* pts3d,
  const float* normals3d,
  const uint32_t* colors,
  uint32_t number_of_points,
  const poisson_reconstruction_screened_parameters& par);

LIBPOISSON_API void poisson_reconstruction_screened(
  std::vector<std::array<double, 3>>& vertices,
  std::vector<std::array<uint32_t, 3>>& triangles,
  const double* pts3d,
  const double* normals3d,
  uint32_t number_of_points,
  const poisson_reconstruction_screened_parameters& par);
  
LIBPOISSON_API void poisson_reconstruction_screened(
  std::vector<std::array<double, 3>>& vertices,
  std::vector<std::array<uint32_t, 3>>& triangles,
  std::vector<uint32_t>& vertex_colors,
  const double* pts3d,
  const double* normals3d,
  const uint32_t* colors,
  uint32_t number_of_points,
  const poisson_reconstruction_screened_parameters& par);

} // namespace libpoisson
