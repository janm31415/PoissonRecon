#include "poisson_reconstruction_screened.h"

#pragma warning(push)
#pragma warning(disable:4100)
#pragma warning(disable:4127)
#pragma warning(disable:4189)
#pragma warning(disable:4244)
#pragma warning(disable:4245)
#pragma warning(disable:4456)
#pragma warning(disable:4458)
#pragma warning(disable:4701)
#pragma warning(disable:4702)
#pragma warning(disable:4706)

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>

#ifdef _WIN32
#include <Windows.h>
#include <Psapi.h>
#endif

#include "poisson_screened/MyTime.h"
#include "poisson_screened/MarchingCubes.h"
#include "poisson_screened/Octree.h"
#include "poisson_screened/SparseMatrix.h"
#include "poisson_screened/CmdLineParser.h"
#include "poisson_screened/PPolynomial.h"
#include "poisson_screened/PlyVertexMini.h"
#include "poisson_screened/MemoryUsage.h"
#ifdef _LIBPOISSON_OPENMP
#include "omp.h"
#endif // _LIBPOISSON_OPENMP
#include "poisson_screened/MultiGridOctreeData.h"

#include <iostream>
#include <cassert>

#include "jtk/vec.h"

namespace libpoisson
{

namespace
{
///////////////////////////////////////////////////////////////////////////////

template <class Real>
Point3D<Real> Convert(const Real* pt)
{
  return Point3D<Real>(pt[0], pt[1], pt[2]);
}

///////////////////////////////////////////////////////////////////////////////

template <class Real>
class MyPointStream : public OrientedPointStream < Real >
{
public:
  MyPointStream(const Real* pts, const Real* normals, uint32_t number_of_points) : p_pts(pts), p_normals(normals), _number_of_points(number_of_points){}
  void reset()
  {
    _current = 0;
  }
  virtual bool nextPoint(OrientedPoint3D< Real >& p)
  {
    if (_current >= _number_of_points)
      return false;
    p.p = Convert<Real>(p_pts + _current*3);
    p.n = Convert<Real>(p_normals + _current*3);
    ++_current;
    return true;
  }
private:
  const Real* p_pts;
  const Real* p_normals;
  uint32_t _number_of_points;
  int _current;
};

///////////////////////////////////////////////////////////////////////////////

template <class Real>
class MyColoredPointStream : public OrientedPointStreamWithData < Real, Point3D<Real> >
{
public:
  MyColoredPointStream(const Real* pts, const Real* normals, const uint32_t* colors, uint32_t number_of_points)
  : p_pts(pts), p_normals(normals), p_colors(colors), _number_of_points(number_of_points) {}
  void reset()
  {
    _current = 0;
  }
  virtual bool nextPoint(OrientedPoint3D< Real >& p, Point3D<Real>& d)
  {
    if (_current >= _number_of_points)
      return false;
    p.p = Convert(p_pts+_current*3);
    p.n = Convert(p_normals+_current*3);
    uint32_t clr = *(p_colors + _current);
    uint32_t red = clr & 255;
    uint32_t green = (clr >> 8) & 255;
    uint32_t blue = (clr >> 16) & 255;
    d[0] = red;
    d[1] = green;
    d[2] = blue;
    ++_current;
    return true;
  }
private:
  const Real* p_pts;
  const Real* p_normals;
  const uint32_t* p_colors;
  uint32_t _number_of_points;
  int _current;
};

///////////////////////////////////////////////////////////////////////////////

void DumpOutput(void* output, const char* format, ...)
{
  if (output)
  {
    char buf[4096];
    va_list marker;
    va_start(marker, format);
    
    vsprintf(buf, format, marker);
    va_end(marker);
    
    *((std::ostream*)output) << buf;
  }
}

void DumpOutput2(void* output, std::vector< char* >& comments, const char* format, ...)
{
  if (output)
  {
    char buf[4096];
    va_list marker;
    va_start(marker, format);
    
    vsprintf(buf, format, marker);
    va_end(marker);
    *((std::ostream*)output) << buf;
  }
}

  template <class T>
  jtk::boundingbox3d<T> bounding_volume_3d(const T* pts3d, uint32_t number_of_points)
    {
    jtk::boundingbox3d<T> result;
    if (number_of_points == 0)
      {
      result.min[0] = (T)1;
      result.min[1] = (T)1;
      result.min[2] = (T)1;
      result.max[0] = (T)0;
      result.max[1] = (T)0;
      result.max[2] = (T)0;
      return result;
      }
    result.min[0] = *(pts3d+0);
    result.min[1] = *(pts3d+1);
    result.min[2] = *(pts3d+2);
    result.max[0] = *(pts3d+0);
    result.max[1] = *(pts3d+1);
    result.max[2] = *(pts3d+2);
    for (uint32_t j = 1; j < number_of_points; ++j)
      {
      const T* pt = pts3d + j*3;
      for (int i = 0; i < 3; ++i)
        {
        if (result.min[i] > (T)(pt[i]))
          result.min[i] = (T)(pt[i]);
        if (result.max[i] < (T)(pt[i]))
          result.max[i] = (T)(pt[i]);
        }
      }
    return result;
    }


template< class Real>
XForm4x4<Real> GetPointStreamScale(const Real* pts3d, uint32_t number_of_points, Real expFact)
{
  
  auto bb = bounding_volume_3d<Real>(pts3d, number_of_points);
  
  Real scale = std::max<Real>(bb.max[0] - bb.min[0], std::max<Real>(bb.max[1] - bb.min[1], bb.max[2] - bb.min[2]));
  scale *= expFact;
  //Real scale = bb.Dim()[bb.MaxDim()] * expFact;
  auto cen = center(bb);
  for (int i = 0; i < 3; i++) cen[i] -= scale / 2;
  XForm4x4< Real > tXForm = XForm4x4< Real >::Identity(), sXForm = XForm4x4< Real >::Identity();
  for (int i = 0; i < 3; i++) sXForm(i, i) = (Real)(1. / scale), tXForm(3, i) = -cen[i];
  return sXForm * tXForm;
}

template <class Real>
void _poisson_reconstruction_screened(
  std::vector<std::array<Real, 3>>& vertices,
  std::vector<std::array<uint32_t, 3>>& triangles,
  const Real* pts3d,
  const Real* normals,
  uint32_t number_of_points,
  const poisson_reconstruction_screened_parameters& par)
{
  vertices.clear();
  triangles.clear();
  int MaxDepthVal = 8;
  int MaxSolveDepthVal = -1;
  int KernelDepthVal = -1;
  int MinDepthVal = 0;
  int FullDepthVal = 5;
  Real SamplesPerNodeVal = 1.5f;
  Real ScaleVal = 1.1f;
  bool ConfidenceFlag = false;
  bool CleanFlag = false;
  bool DensityFlag = false;
  Real PointWeightVal = 4.f;
  int AdaptiveExponentVal = 1;
  int BoundaryTypeVal = 1;
  bool CompleteFlag = false;
  bool NonManifoldFlag = false;
  bool ShowResidualFlag = false;
  int CGDepthVal = 0;
  int ItersVal = 8;
  Real CSSolverAccuracyVal = 1e-3f;
  
  bool VerboseFlag = false;
  int ThreadsVal = omp_get_num_procs();
  bool LinearFitFlag = false;
  Real LowResIterMultiplierVal = 1.f;
  Real ColorVal = 16.f;
  
#define Degree 2
  
  typedef typename Octree< Real >::template DensityEstimator< WEIGHT_DEGREE > DensityEstimator;
  typedef typename Octree< Real >::template InterpolationInfo< false > InterpolationInfo;
  //typedef OrientedPointStreamWithData< Real, Point3D< Real > > PointStreamWithData;
  typedef TransformedOrientedPointStream< Real > XPointStream;
  Reset< Real >();
  std::vector< char* > comments;
  
  XForm4x4< Real > xForm = GetPointStreamScale<Real>(pts3d, number_of_points, ScaleVal);
  XForm4x4< Real > iXForm = xForm.inverse();
  DumpOutput2(par.output_stream, comments, "Running Screened Poisson Reconstruction (Version 9.0)\n");
  double startTime = Time();
  
  OctNode< TreeNodeData >::SetAllocator(MEMORY_ALLOCATOR_BLOCK_SIZE);
  Octree< Real > tree;
  //OctreeProfiler< Real > profiler(tree);
  tree.threads = ThreadsVal;
  if (MaxSolveDepthVal < 0) MaxSolveDepthVal = MaxDepthVal;
  
  
  
  
  //	int kernelDepth = KernelDepth.set ? KernelDepth.value : Depth.value-2;
  if (KernelDepthVal < 0) KernelDepthVal = MaxDepthVal - 2;
  if (KernelDepthVal > MaxDepthVal)
  {
    printf("kernelDepth cannot be greateer Depth.value\n");
    return;
  }
  
  int pointCount;
  
  Real pointWeightSum;
  std::vector< typename Octree< Real >::PointSample >* samples = new std::vector< typename Octree< Real >::PointSample >();
  std::vector< ProjectiveData< Point3D< Real >, Real > >* sampleData = NULL;
  DensityEstimator* density = NULL;
  SparseNodeData< Point3D< Real >, NORMAL_DEGREE >* normalInfo = NULL;
  Real targetValue = (Real)0.5;
  
  // Read in the samples (and color data)
  {
    
    sampleData = nullptr;// new std::vector< ProjectiveData< Point3D< Real >, Real > >();
    MyPointStream my_pointStream(pts3d, normals, number_of_points);
    my_pointStream.reset();
    XPointStream _pointStream(xForm, my_pointStream);
    pointCount = tree.template init< Point3D< Real > >(_pointStream, MaxDepthVal, ConfidenceFlag, *samples, sampleData);
    
#pragma omp parallel for num_threads( ThreadsVal )
    for (int i = 0; i < (int)samples->size(); i++) (*samples)[i].sample.data.n *= (Real)-1;
    
    DumpOutput(par.output_stream, "Input Points / Samples: %d / %d\n", pointCount, samples->size());
    //profiler.dumpOutput2(comments, "# Read input into tree:");
  }
  
  DenseNodeData< Real, Degree > solution;
  
  {
    DenseNodeData< Real, Degree > constraints;
    InterpolationInfo* iInfo = NULL;
    int solveDepth = MaxSolveDepthVal;
    
    tree.resetNodeIndices();
    
    // Get the kernel density estimator [If discarding, compute anew. Otherwise, compute once.]
    {
      //profiler.start();
      density = tree.template setDensityEstimator< WEIGHT_DEGREE >(*samples, KernelDepthVal, SamplesPerNodeVal);
      //profiler.dumpOutput2(comments, "#   Got kernel density:");
    }
    
    // Transform the Hermite samples into a vector field [If discarding, compute anew. Otherwise, compute once.]
    {
      //profiler.start();
      normalInfo = new SparseNodeData< Point3D< Real >, NORMAL_DEGREE >();
      *normalInfo = tree.template setNormalField< NORMAL_DEGREE >(*samples, *density, pointWeightSum, BOUNDARY_NEUMANN == BOUNDARY_NEUMANN);
      //profiler.dumpOutput2(comments, "#     Got normal field:");
    }
    
    if (!DensityFlag) delete density, density = NULL;
    
    // Trim the tree and prepare for multigrid
    {
      //profiler.start();
      std::vector< int > indexMap;
      
      constexpr int MAX_DEGREE = NORMAL_DEGREE > Degree ? NORMAL_DEGREE : Degree;
      tree.template inalizeForBroodedMultigrid< MAX_DEGREE, Degree, BOUNDARY_NEUMANN >(FullDepthVal, typename Octree< Real >::template HasNormalDataFunctor< NORMAL_DEGREE >(*normalInfo), &indexMap);
      
      if (normalInfo) normalInfo->remapIndices(indexMap);
      if (density) density->remapIndices(indexMap);
      //profiler.dumpOutput2(comments, "#       Finalized tree:");
    }
    
    // Add the FEM constraints
    {
      //profiler.start();
      constraints = tree.template initDenseNodeData< Degree >();
      tree.template addFEMConstraints< Degree, BOUNDARY_NEUMANN, NORMAL_DEGREE, BOUNDARY_NEUMANN >(FEMVFConstraintFunctor< NORMAL_DEGREE, BOUNDARY_NEUMANN, Degree, BOUNDARY_NEUMANN >(1., 0.), *normalInfo, constraints, solveDepth);
      //profiler.dumpOutput2(comments, "#  Set FEM constraints:");
    }
    
    // Free up the normal info [If we don't need it for subseequent iterations.]
    delete normalInfo, normalInfo = NULL;
    
    // Add the interpolation constraints
    if (PointWeightVal > 0)
    {
      //profiler.start();
      iInfo = new InterpolationInfo(tree, *samples, targetValue, AdaptiveExponentVal, (Real)PointWeightVal * pointWeightSum, (Real)0);
      tree.template addInterpolationConstraints< Degree, BOUNDARY_NEUMANN >(*iInfo, constraints, solveDepth);
      //profiler.dumpOutput2(comments, "#Set point constraints:");
    }
    
    DumpOutput(par.output_stream, "Leaf Nodes / Active Nodes / Ghost Nodes: %d / %d / %d\n", (int)tree.leaves(), (int)tree.nodes(), (int)tree.ghostNodes());
    DumpOutput(par.output_stream, "Memory Usage: %.3f MB\n", Real(MemoryInfo::Usage()) / (1 << 20));
    
    // Solve the linear system
    {
      //profiler.start();
      typename Octree< Real >::SolverInfo solverInfo;
      solverInfo.cgDepth = CGDepthVal, solverInfo.iters = ItersVal, solverInfo.cgAccuracy = CSSolverAccuracyVal, solverInfo.verbose = VerboseFlag, solverInfo.showResidual = ShowResidualFlag, solverInfo.lowResIterMultiplier = std::max< double >(1., LowResIterMultiplierVal);
      solution = tree.template solveSystem< Degree, BOUNDARY_NEUMANN >(FEMSystemFunctor< Degree, BOUNDARY_NEUMANN >(0, 1., 0), iInfo, constraints, solveDepth, solverInfo);
      //profiler.dumpOutput2(comments, "# Linear system solved:");
      if (iInfo) delete iInfo, iInfo = NULL;
    }
  }
  
  typedef PlyVertex< Real > Vertex;
  
  CoredVectorMeshData< Vertex > mesh;
  
  {
    //profiler.start();
    double valueSum = 0, weightSum = 0;
    typename Octree< Real >::template MultiThreadedEvaluator< Degree, BOUNDARY_NEUMANN > evaluator(&tree, solution, ThreadsVal);
#pragma omp parallel for num_threads( ThreadsVal ) reduction( + : valueSum , weightSum )
    for (int j = 0; j < samples->size(); j++)
    {
      ProjectiveData< OrientedPoint3D< Real >, Real >& sample = (*samples)[j].sample;
      Real w = sample.weight;
      if (w > 0) weightSum += w, valueSum += evaluator.value(sample.data.p / sample.weight, omp_get_thread_num(), (*samples)[j].node) * w;
    }
    Real isoValue = (Real)(valueSum / weightSum);
    //		if( samples ) delete samples , samples = NULL;
    //profiler.dumpOutput("Got average:");
    DumpOutput(par.output_stream, "Iso-Value: %e\n", isoValue);
    
    //profiler.start();
    SparseNodeData< ProjectiveData< Point3D< Real >, Real >, DATA_DEGREE >* colorData = NULL;
    if (sampleData)
    {
      colorData = new SparseNodeData< ProjectiveData< Point3D< Real >, Real >, DATA_DEGREE >();
      *colorData = tree.template setDataField< DATA_DEGREE, false >(*samples, *sampleData, (DensityEstimator*)NULL);
      delete sampleData, sampleData = NULL;
      for (const OctNode< TreeNodeData >* n = tree.tree().nextNode(); n; n = tree.tree().nextNode(n))
      {
        ProjectiveData< Point3D< Real >, Real >* clr = (*colorData)(n);
        if (clr)
          (*clr) *= (Real)pow(ColorVal, tree.depth(n));
      }
    }
    tree.template getMCIsoSurface< Degree, BOUNDARY_NEUMANN, WEIGHT_DEGREE, DATA_DEGREE >(density, colorData, solution, isoValue, mesh, !LinearFitFlag, !NonManifoldFlag, false /*PolygonMesh.set*/);
    DumpOutput(par.output_stream, "Vertices / Polygons: %d / %d\n", mesh.outOfCorePointCount() + mesh.inCorePoints.size(), mesh.polygonCount());
    //profiler.dumpOutput2(comments, "#        Got triangles:");
  }
  
  //        FreePointer( solution );
  
  
  mesh.resetIterator();
  int nr_vertices = int(mesh.outOfCorePointCount() + mesh.inCorePoints.size());
  int nr_faces = mesh.polygonCount();
  
  vertices.clear();
  triangles.clear();
  
  for (int i = 0; i<int(mesh.inCorePoints.size()); ++i)
  {
    PlyVertex< Real > vertex = mesh.inCorePoints[i];
    
    Point3D<Real> pp = iXForm * vertex.point;
    
    vertices.push_back({{pp[0], pp[1], pp[2]}});
  }
  for (int i = 0; i < mesh.outOfCorePointCount(); ++i)
  {
    PlyVertex< Real > vertex;
    mesh.nextOutOfCorePoint(vertex);
    Point3D<Real> pp = iXForm * vertex.point;
    
    vertices.push_back({{pp[0], pp[1], pp[2]}});
  }
  std::vector< CoredVertexIndex > polygon;
  while (mesh.nextPolygon(polygon))
  {
    assert(polygon.size() == 3);
    int indV[3];
    for (int i = 0; i<int(polygon.size()); i++)
    {
      if (polygon[i].inCore) indV[i] = polygon[i].idx;
      else                    indV[i] = polygon[i].idx + int(mesh.inCorePoints.size());
    }
    std::array<uint32_t, 3> tria = {{(uint32_t)indV[0], (uint32_t)indV[1], (uint32_t)indV[2]}};
    triangles.push_back(tria);
  }
  
  if (density) delete density, density = NULL;
  DumpOutput2(par.output_stream, comments, "#          Total Solve: %9.1f (s), %9.1f (MB)\n", Time() - startTime, tree.maxMemoryUsage());
  
  
}


template <class Real>
void _poisson_reconstruction_screened(
  std::vector<std::array<Real, 3>>& vertices,
  std::vector<std::array<uint32_t, 3>>& triangles,
  std::vector<uint32_t>& vertex_colors,
  const Real* pts3d,
  const Real* normals,
  const uint32_t* colors,
  uint32_t number_of_points,
  const poisson_reconstruction_screened_parameters& par)
{
  vertices.clear();
  triangles.clear();
  vertex_colors.clear();
  int MaxDepthVal = 8;
  int MaxSolveDepthVal = -1;
  int KernelDepthVal = -1;
  int MinDepthVal = 0;
  int FullDepthVal = 5;
  Real SamplesPerNodeVal = 1.5f;
  Real ScaleVal = 1.1f;
  bool ConfidenceFlag = false;
  bool CleanFlag = false;
  bool DensityFlag = false;
  Real PointWeightVal = 4.f;
  int AdaptiveExponentVal = 1;
  int BoundaryTypeVal = 1;
  bool CompleteFlag = false;
  bool NonManifoldFlag = false;
  bool ShowResidualFlag = false;
  int CGDepthVal = 0;
  int ItersVal = 8;
  Real CSSolverAccuracyVal = 1e-3f;
  
  bool VerboseFlag = false;
  int ThreadsVal = omp_get_num_procs();
  bool LinearFitFlag = false;
  Real LowResIterMultiplierVal = 1.f;
  Real ColorVal = 16.f;
  
#define Degree 2
  
  typedef typename Octree< Real >::template DensityEstimator< WEIGHT_DEGREE > DensityEstimator;
  typedef typename Octree< Real >::template InterpolationInfo< false > InterpolationInfo;
  //typedef OrientedPointStreamWithData< Real, Point3D< Real > > PointStreamWithData;
  typedef TransformedOrientedPointStreamWithData< Real, Point3D<Real> > XPointStreamWithData;
  Reset< Real >();
  std::vector< char* > comments;
  
  XForm4x4< Real > xForm = GetPointStreamScale<Real>(pts3d, number_of_points, ScaleVal);
  XForm4x4< Real > iXForm = xForm.inverse();
  DumpOutput2(par.output_stream, comments, "Running Screened Poisson Reconstruction (Version 9.0)\n");
  double startTime = Time();
  
  OctNode< TreeNodeData >::SetAllocator(MEMORY_ALLOCATOR_BLOCK_SIZE);
  Octree< Real > tree;
  //OctreeProfiler< Real > profiler(tree);
  tree.threads = ThreadsVal;
  if (MaxSolveDepthVal < 0) MaxSolveDepthVal = MaxDepthVal;
  
  
  
  
  //	int kernelDepth = KernelDepth.set ? KernelDepth.value : Depth.value-2;
  if (KernelDepthVal < 0) KernelDepthVal = MaxDepthVal - 2;
  if (KernelDepthVal > MaxDepthVal)
  {
    printf("kernelDepth cannot be greateer Depth.value\n");
    return;
  }
  
  int pointCount;
  
  Real pointWeightSum;
  std::vector< typename Octree< Real >::PointSample >* samples = new std::vector< typename Octree< Real >::PointSample >();
  std::vector< ProjectiveData< Point3D< Real >, Real > >* sampleData = NULL;
  DensityEstimator* density = NULL;
  SparseNodeData< Point3D< Real >, NORMAL_DEGREE >* normalInfo = NULL;
  Real targetValue = (Real)0.5;
  
  // Read in the samples (and color data)
  {
    
    sampleData = new std::vector< ProjectiveData< Point3D< Real >, Real > >();
    MyColoredPointStream my_pointStream(pts3d, normals, colors, number_of_points);
    my_pointStream.reset();
    XPointStreamWithData _pointStream(xForm, my_pointStream);
    pointCount = tree.template init< Point3D< Real > >(_pointStream, MaxDepthVal, ConfidenceFlag, *samples, sampleData);
    
#pragma omp parallel for num_threads( ThreadsVal )
    for (int i = 0; i < (int)samples->size(); i++) (*samples)[i].sample.data.n *= (Real)-1;
    
    DumpOutput(par.output_stream, "Input Points / Samples: %d / %d\n", pointCount, samples->size());
    //profiler.dumpOutput2(comments, "# Read input into tree:");
  }
  
  DenseNodeData< Real, Degree > solution;
  
  {
    DenseNodeData< Real, Degree > constraints;
    InterpolationInfo* iInfo = NULL;
    int solveDepth = MaxSolveDepthVal;
    
    tree.resetNodeIndices();
    
    // Get the kernel density estimator [If discarding, compute anew. Otherwise, compute once.]
    {
      //profiler.start();
      density = tree.template setDensityEstimator< WEIGHT_DEGREE >(*samples, KernelDepthVal, SamplesPerNodeVal);
      //profiler.dumpOutput2(comments, "#   Got kernel density:");
    }
    
    // Transform the Hermite samples into a vector field [If discarding, compute anew. Otherwise, compute once.]
    {
      //profiler.start();
      normalInfo = new SparseNodeData< Point3D< Real >, NORMAL_DEGREE >();
      *normalInfo = tree.template setNormalField< NORMAL_DEGREE >(*samples, *density, pointWeightSum, BOUNDARY_NEUMANN == BOUNDARY_NEUMANN);
      //profiler.dumpOutput2(comments, "#     Got normal field:");
    }
    
    if (!DensityFlag) delete density, density = NULL;
    
    // Trim the tree and prepare for multigrid
    {
      //profiler.start();
      std::vector< int > indexMap;
      
      constexpr int MAX_DEGREE = NORMAL_DEGREE > Degree ? NORMAL_DEGREE : Degree;
      tree.template inalizeForBroodedMultigrid< MAX_DEGREE, Degree, BOUNDARY_NEUMANN >(FullDepthVal, typename Octree< Real >::template HasNormalDataFunctor< NORMAL_DEGREE >(*normalInfo), &indexMap);
      
      if (normalInfo) normalInfo->remapIndices(indexMap);
      if (density) density->remapIndices(indexMap);
      //profiler.dumpOutput2(comments, "#       Finalized tree:");
    }
    
    // Add the FEM constraints
    {
      //profiler.start();
      constraints = tree.template initDenseNodeData< Degree >();
      tree.template addFEMConstraints< Degree, BOUNDARY_NEUMANN, NORMAL_DEGREE, BOUNDARY_NEUMANN >(FEMVFConstraintFunctor< NORMAL_DEGREE, BOUNDARY_NEUMANN, Degree, BOUNDARY_NEUMANN >(1., 0.), *normalInfo, constraints, solveDepth);
      //profiler.dumpOutput2(comments, "#  Set FEM constraints:");
    }
    
    // Free up the normal info [If we don't need it for subseequent iterations.]
    delete normalInfo, normalInfo = NULL;
    
    // Add the interpolation constraints
    if (PointWeightVal > 0)
    {
      //profiler.start();
      iInfo = new InterpolationInfo(tree, *samples, targetValue, AdaptiveExponentVal, (Real)PointWeightVal * pointWeightSum, (Real)0);
      tree.template addInterpolationConstraints< Degree, BOUNDARY_NEUMANN >(*iInfo, constraints, solveDepth);
      //profiler.dumpOutput2(comments, "#Set point constraints:");
    }
    
    DumpOutput(par.output_stream, "Leaf Nodes / Active Nodes / Ghost Nodes: %d / %d / %d\n", (int)tree.leaves(), (int)tree.nodes(), (int)tree.ghostNodes());
    DumpOutput(par.output_stream, "Memory Usage: %.3f MB\n", Real(MemoryInfo::Usage()) / (1 << 20));
    
    // Solve the linear system
    {
      //profiler.start();
      typename Octree< Real >::SolverInfo solverInfo;
      solverInfo.cgDepth = CGDepthVal, solverInfo.iters = ItersVal, solverInfo.cgAccuracy = CSSolverAccuracyVal, solverInfo.verbose = VerboseFlag, solverInfo.showResidual = ShowResidualFlag, solverInfo.lowResIterMultiplier = std::max< double >(1., LowResIterMultiplierVal);
      solution = tree.template solveSystem< Degree, BOUNDARY_NEUMANN >(FEMSystemFunctor< Degree, BOUNDARY_NEUMANN >(0, 1., 0), iInfo, constraints, solveDepth, solverInfo);
      //profiler.dumpOutput2(comments, "# Linear system solved:");
      if (iInfo) delete iInfo, iInfo = NULL;
    }
  }
  
  typedef PlyColorVertex< Real > Vertex;
  
  CoredVectorMeshData< Vertex > mesh;
  
  {
    //profiler.start();
    double valueSum = 0, weightSum = 0;
    typename Octree< Real >::template MultiThreadedEvaluator< Degree, BOUNDARY_NEUMANN > evaluator(&tree, solution, ThreadsVal);
#pragma omp parallel for num_threads( ThreadsVal ) reduction( + : valueSum , weightSum )
    for (int j = 0; j < samples->size(); j++)
    {
      ProjectiveData< OrientedPoint3D< Real >, Real >& sample = (*samples)[j].sample;
      Real w = sample.weight;
      if (w > 0) weightSum += w, valueSum += evaluator.value(sample.data.p / sample.weight, omp_get_thread_num(), (*samples)[j].node) * w;
    }
    Real isoValue = (Real)(valueSum / weightSum);
    //		if( samples ) delete samples , samples = NULL;
    //profiler.dumpOutput("Got average:");
    DumpOutput(par.output_stream, "Iso-Value: %e\n", isoValue);
    
    //profiler.start();
    SparseNodeData< ProjectiveData< Point3D< Real >, Real >, DATA_DEGREE >* colorData = NULL;
    if (sampleData)
    {
      colorData = new SparseNodeData< ProjectiveData< Point3D< Real >, Real >, DATA_DEGREE >();
      *colorData = tree.template setDataField< DATA_DEGREE, false >(*samples, *sampleData, (DensityEstimator*)NULL);
      delete sampleData, sampleData = NULL;
      for (const OctNode< TreeNodeData >* n = tree.tree().nextNode(); n; n = tree.tree().nextNode(n))
      {
        ProjectiveData< Point3D< Real >, Real >* clr = (*colorData)(n);
        if (clr)
          (*clr) *= (Real)pow(ColorVal, tree.depth(n));
      }
    }
    tree.template getMCIsoSurface< Degree, BOUNDARY_NEUMANN, WEIGHT_DEGREE, DATA_DEGREE >(density, colorData, solution, isoValue, mesh, !LinearFitFlag, !NonManifoldFlag, false /*PolygonMesh.set*/);
    DumpOutput(par.output_stream, "Vertices / Polygons: %d / %d\n", mesh.outOfCorePointCount() + mesh.inCorePoints.size(), mesh.polygonCount());
    //profiler.dumpOutput2(comments, "#        Got triangles:");
    
    if (colorData) delete colorData, colorData = NULL;
  }
  
  //        FreePointer( solution );
  
  
  mesh.resetIterator();
  int nr_vertices = int(mesh.outOfCorePointCount() + mesh.inCorePoints.size());
  int nr_faces = mesh.polygonCount();
  
  vertices.clear();
  triangles.clear();
  vertex_colors.clear();
  
  for (int i = 0; i<int(mesh.inCorePoints.size()); ++i)
  {
    Vertex vertex = mesh.inCorePoints[i];
    
    Point3D<Real> pp = iXForm * vertex.point;
    
    vertices.push_back({{pp[0], pp[1], pp[2]}});
    
    uint32_t clr = 0xff000000 | ((uint32_t)vertex.color[2] << 16) | ((uint32_t)vertex.color[1] << 8) | ((uint32_t)vertex.color[0]);
    
    vertex_colors.push_back(clr);
  }
  for (int i = 0; i < mesh.outOfCorePointCount(); ++i)
  {
    Vertex vertex;
    mesh.nextOutOfCorePoint(vertex);
    Point3D<Real> pp = iXForm * vertex.point;
    
    vertices.push_back({{pp[0], pp[1], pp[2]}});
    
    uint32_t clr = 0xff000000 | ((uint32_t)vertex.color[2] << 16) | ((uint32_t)vertex.color[1] << 8) | ((uint32_t)vertex.color[0]);
    vertex_colors.push_back(clr);
  }
  std::vector< CoredVertexIndex > polygon;
  while (mesh.nextPolygon(polygon))
  {
    assert(polygon.size() == 3);
    int indV[3];
    for (int i = 0; i<int(polygon.size()); i++)
    {
      if (polygon[i].inCore) indV[i] = polygon[i].idx;
      else                    indV[i] = polygon[i].idx + int(mesh.inCorePoints.size());
    }
    std::array<uint32_t, 3> tria = {{(uint32_t)indV[0], (uint32_t)indV[1], (uint32_t)indV[2]}};
    triangles.push_back(tria);
  }
  
  if (density) delete density, density = NULL;
  DumpOutput2(par.output_stream, comments, "#          Total Solve: %9.1f (s), %9.1f (MB)\n", Time() - startTime, tree.maxMemoryUsage());
  
  
  
  
}

}


void poisson_reconstruction_screened(
  std::vector<std::array<float, 3>>& vertices,
  std::vector<std::array<uint32_t, 3>>& triangles,
  const float* pts3d,
  const float* normals3d,
  uint32_t number_of_points,
  const poisson_reconstruction_screened_parameters& par) {
  _poisson_reconstruction_screened(vertices, triangles, pts3d, normals3d, number_of_points, par);
}

void poisson_reconstruction_screened(
  std::vector<std::array<float, 3>>& vertices,
  std::vector<std::array<uint32_t, 3>>& triangles,
  std::vector<uint32_t>& vertex_colors,
  const float* pts3d,
  const float* normals3d,
  const uint32_t* colors,
  uint32_t number_of_points,
  const poisson_reconstruction_screened_parameters& par) {
  _poisson_reconstruction_screened(vertices, triangles, vertex_colors, pts3d, normals3d, colors, number_of_points, par);
}

void poisson_reconstruction_screened(
  std::vector<std::array<double, 3>>& vertices,
  std::vector<std::array<uint32_t, 3>>& triangles,
  const double* pts3d,
  const double* normals3d,
  uint32_t number_of_points,
  const poisson_reconstruction_screened_parameters& par) {
  _poisson_reconstruction_screened(vertices, triangles, pts3d, normals3d, number_of_points, par);
}
  
void poisson_reconstruction_screened(
  std::vector<std::array<double, 3>>& vertices,
  std::vector<std::array<uint32_t, 3>>& triangles,
  std::vector<uint32_t>& vertex_colors,
  const double* pts3d,
  const double* normals3d,
  const uint32_t* colors,
  uint32_t number_of_points,
  const poisson_reconstruction_screened_parameters& par) {
  _poisson_reconstruction_screened(vertices, triangles, vertex_colors, pts3d, normals3d, colors, number_of_points, par);
}
}
#pragma warning(pop)
