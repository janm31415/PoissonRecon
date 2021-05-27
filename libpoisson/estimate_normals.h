#pragma once

#include "jtk/vec.h"
#include <vector>

namespace libpoisson
{


std::vector<jtk::vec3<float>> estimate_normals(const std::vector<jtk::vec3<float>>& pts, uint32_t k);

}
