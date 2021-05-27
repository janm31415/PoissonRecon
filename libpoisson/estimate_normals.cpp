#include "estimate_normals.h"

#include "jtk/point_tree.h"
#include "jtk/containers.h"
#include "jtk/fitting.h"

namespace libpoisson
{

namespace {
  struct tree_point
    {
    jtk::vec3<float> pt;
    uint32_t idx;
    float operator [](int i) const
      {
      return pt[i];
      }
    float& operator [](int i)
      {
      return pt[i];
      }
    };
    
  struct point_tree_traits
    {
    typedef float value_type;
    enum { dimension = 3 };
    typedef tree_point point_type;
    };
    
  uint64_t make_edge(uint32_t v0, uint32_t v1)
    {
    uint64_t e = v0;
    e <<= 32;
    e |= (uint64_t)v1;
    return e;
    }

  uint32_t edge_to_v0(uint64_t e)
    {
    e >>= 32;
    return (uint32_t)(e & 0xffffffff);
    }

  uint32_t edge_to_v1(uint64_t e)
    {
    return (uint32_t)(e & 0xffffffff);
    }
}

std::vector<jtk::vec3<float>> estimate_normals(const std::vector<jtk::vec3<float>>& pts, uint32_t k) {
  std::vector<jtk::vec3<float>> normals;
  normals.reserve(pts.size());

  jtk::point_tree<point_tree_traits> tree;
  std::vector<tree_point> vert;
  vert.reserve(pts.size());
  for (uint32_t v = 0; v < pts.size(); ++v)
    {
    tree_point pt;
    pt.pt = pts[v];
    pt.idx = v;
    vert.push_back(pt);
    }
  tree.efficient_build_tree(vert.begin(), vert.end());

  for (uint32_t v = 0; v < pts.size(); ++v)
    {
    tree_point tp;
    tp.pt = pts[v];;
    std::vector < tree_point > pts = tree.find_k_nearest((int)k, tp);
    jtk::vec3<float> origin, normal;
    float eig;
    std::vector<jtk::vec3<float>> raw_pts;
    raw_pts.reserve(pts.size());
    for (auto p : pts)
      raw_pts.push_back(p.pt);
    jtk::fit_plane(origin, normal, eig, raw_pts);
    normals.push_back(normal);
    }

  std::vector<bool> treated(pts.size(), false);

  jtk::hashed_heap<uint64_t, float> heap;

  uint32_t v = 0;
  while (true)
    {
    while (v < pts.size() && treated[v])
      ++v;
    if (v == pts.size())
      break;

    treated[v] = true;

    tree_point tp;
    tp.pt = pts[v];
    std::vector < tree_point > kpts = tree.find_k_nearest((int)k, tp);
    for (const auto& pt : kpts)
      {
      if (pt.idx != v && !treated[pt.idx])
        {
        float score = fabs(dot(normals[v], normals[pt.idx]));
        heap.push(std::pair<uint64_t, float>(make_edge(v, pt.idx), score));
        }
      }

    while (!heap.empty())
      {
      auto top_element = heap.top();
      heap.pop();
      uint32_t v0 = edge_to_v0(top_element.first);
      uint32_t v1 = edge_to_v1(top_element.first);
      if (!treated[v1])
        {
        treated[v1] = true;
        if (dot(normals[v0], normals[v1]) < 0.f)
          normals[v1] = -normals[v1];

        tp.pt = pts[v1];
        std::vector < tree_point > pts = tree.find_k_nearest((int)k, tp);
        for (const auto& pt : pts)
          {
          if (pt.idx != v1 && !treated[pt.idx])
            {
            float score = fabs(dot(normals[v1], normals[pt.idx]));
            heap.push(std::pair<uint64_t, float>(make_edge(v1, pt.idx), score));
            }
          }

        }
      }
    }

  return normals;
}

}
