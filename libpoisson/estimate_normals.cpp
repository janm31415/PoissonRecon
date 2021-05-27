#include "estimate_normals.h"

#include "jtk/vec.h"
#include "jtk/point_tree.h"
#include "jtk/containers.h"
#include "jtk/fitting.h"

namespace libpoisson
{

namespace {
template <class T>
struct tree_point
{
  jtk::vec3<T> pt;
  uint32_t idx;
  T operator [](int i) const
  {
    return pt[i];
  }
  T& operator [](int i)
  {
    return pt[i];
  }
};

template <class T>
struct point_tree_traits
{
  typedef T value_type;
  enum { dimension = 3 };
  typedef tree_point<T> point_type;
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

template <class T>
std::array<T, 3> vec_to_array(const jtk::vec3<T>& p)
{
  std::array<T, 3> arr = {{p[0], p[1], p[2]}};
  return arr;
}

template <class T>
T dot(const std::array<T, 3>& left, const std::array<T, 3>& right) {
  return left[0]*right[0]+left[1]*right[1]+left[2]*right[2];
}

template <class T>
std::vector<std::array<T, 3>> _estimate_normals(const T* pts3d, uint32_t number_of_points, uint32_t k) {
  std::vector<std::array<T, 3>> normals;
  normals.reserve(number_of_points);
  
  jtk::point_tree<point_tree_traits<T>> tree;
  std::vector<tree_point<T>> vert;
  vert.reserve(number_of_points);
  for (uint32_t v = 0; v < number_of_points; ++v)
  {
    tree_point<T> pt;
    pt.pt[0] = *(pts3d+v*3);
    pt.pt[1] = *(pts3d+v*3 + 1);
    pt.pt[2] = *(pts3d+v*3 + 2);
    pt.idx = v;
    vert.push_back(pt);
  }
  tree.efficient_build_tree(vert.begin(), vert.end());
  
  for (uint32_t v = 0; v < number_of_points; ++v)
  {
    tree_point<T> tp;
    tp.pt[0] = *(pts3d+v*3);
    tp.pt[1] = *(pts3d+v*3 + 1);
    tp.pt[2] = *(pts3d+v*3 + 2);
    std::vector < tree_point<T> > pts = tree.find_k_nearest((int)k, tp);
    jtk::vec3<T> origin, normal;
    T eig;
    std::vector<jtk::vec3<T>> raw_pts;
    raw_pts.reserve(pts.size());
    for (auto p : pts)
      raw_pts.push_back(p.pt);
    jtk::fit_plane(origin, normal, eig, raw_pts);
    normals.push_back(vec_to_array(normal));
  }
  
  std::vector<bool> treated(number_of_points, false);
  
  jtk::hashed_heap<uint64_t, float> heap;
  
  uint32_t v = 0;
  while (true)
  {
    while (v < number_of_points && treated[v])
      ++v;
    if (v == number_of_points)
      break;
    
    treated[v] = true;
    
    tree_point<T> tp;
    tp.pt[0] = *(pts3d+v*3);
    tp.pt[1] = *(pts3d+v*3 + 1);
    tp.pt[2] = *(pts3d+v*3 + 2);
    std::vector < tree_point<T> > kpts = tree.find_k_nearest((int)k, tp);
    for (const auto& pt : kpts)
    {
      if (pt.idx != v && !treated[pt.idx])
      {
        T score = fabs(dot(normals[v], normals[pt.idx]));
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
        if (dot(normals[v0], normals[v1]) < 0)
        {
          normals[v1][0] = -normals[v1][0];
          normals[v1][1] = -normals[v1][1];
          normals[v1][2] = -normals[v1][2];
        }
        
        tp.pt[0] = *(pts3d+v1*3);
        tp.pt[1] = *(pts3d+v1*3 + 1);
        tp.pt[2] = *(pts3d+v1*3 + 2);
        std::vector < tree_point<T> > pts = tree.find_k_nearest((int)k, tp);
        for (const auto& pt : pts)
        {
          if (pt.idx != v1 && !treated[pt.idx])
          {
            T score = fabs(dot(normals[v1], normals[pt.idx]));
            heap.push(std::pair<uint64_t, float>(make_edge(v1, pt.idx), score));
          }
        }
        
      }
    }
  }
  
  return normals;
}

}

std::vector<std::array<float, 3>> estimate_normals(const float* pts3d, uint32_t number_of_points, uint32_t k) {
  return _estimate_normals<float>(pts3d, number_of_points, k);
}

std::vector<std::array<double, 3>> estimate_normals(const double* pts3d, uint32_t number_of_points, uint32_t k) {
  return _estimate_normals<double>(pts3d, number_of_points, k);
}
}
