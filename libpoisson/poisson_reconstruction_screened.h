#pragma once

#include "jtk/vec.h"
#include "libpoisson_api.h"
#include <stdint.h>
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

LIBPOISSON_API void poisson_reconstruction_screened(std::vector<jtk::vec3<float>>& vertices, std::vector<jtk::vec3<uint32_t>>& triangles, const std::vector<jtk::vec3<float>>& pts, const std::vector<jtk::vec3<float>>& normals, const poisson_reconstruction_screened_parameters& par);

LIBPOISSON_API void poisson_reconstruction_screened(std::vector<jtk::vec3<float>>& vertices, std::vector<jtk::vec3<uint32_t>>& triangles, std::vector<jtk::vec3<float>>& vertex_colors, const std::vector<jtk::vec3<float>>& pts, const std::vector<jtk::vec3<float>>& normals, const std::vector<uint32_t>& colors, const poisson_reconstruction_screened_parameters& par);

LIBPOISSON_API void poisson_reconstruction_screened(std::vector<jtk::vec3<double>>& vertices, std::vector<jtk::vec3<uint32_t>>& triangles, const std::vector<jtk::vec3<double>>& pts, const std::vector<jtk::vec3<double>>& normals, const poisson_reconstruction_screened_parameters& par);

LIBPOISSON_API void poisson_reconstruction_screened(std::vector<jtk::vec3<double>>& vertices, std::vector<jtk::vec3<uint32_t>>& triangles, std::vector<jtk::vec3<double>>& vertex_colors, const std::vector<jtk::vec3<double>>& pts, const std::vector<jtk::vec3<double>>& normals, const std::vector<uint32_t>& colors, const poisson_reconstruction_screened_parameters& par);

} // namespace libpoisson