#pragma once

#include <array>
#include <vector>

namespace libpoisson
{

std::vector<std::array<float, 3>> estimate_normals(const float* pts3d, uint32_t number_of_points, uint32_t k=10);

std::vector<std::array<double, 3>> estimate_normals(const double* pts3d, uint32_t number_of_points, uint32_t k=10);

}
