#ifndef DC_3D_H
#define DC_3D_H

#include "util.h"

#include <map>
#include <vector>
#include <set>
#include <list>

class DepthComplexity2D;

// Usage:
//   TriMesh mesh = ...;
//   int discretSteps = 10;
//   int fboWidth = 400;
//   int fboHeight = 400;
//   DepthComplexity3D cd3d(mesh, fboWidth, fboHeight, discretSteps);  // Computation is done.
//   int maximum = dc3d.maximum();
//   int minimum = dc3d.minimum();
//   c3d.writeHistogram(std::cout);
//   c3d.writeRays(std::cout);
class DepthComplexity3D {
public:
  DepthComplexity3D(int fboWidth, int fboHeight, int discretSteps);
  virtual ~DepthComplexity3D();

  void process(const TriMesh &mesh);

  void setComputeHistogram(bool computeHistogram);
  void setComputeMaximumRays(bool computeMaximumRays);
  void setComputeGoodRays(bool computeGoodRays);
  void setThreshold(unsigned threshold);

  unsigned int maximum() const { return _maximum; }
  const std::set<Segment, classcomp> &maximumRays() const { return _maximumRays; }
  const std::set<Segment, classcomp> &goodRays(int intersect) const { return _goodRays[intersect]; }
  const std::vector<Segment> &usedPlanes() const { return _usedPlanes; }
  const std::vector<Point> &intersectionPoints() const { return _intersectionPoints; }
  unsigned int getThreshold() { return _threshold; }
  const std::list<unsigned int> &intersectionTris(int rayIndex){return *_intersectionTriList.find(rayIndex)->second;}

  void writeHistogram(std::ostream& out);
  void writeRays(std::ostream& out);
  void writeRays(std::ostream& out, const std::set<Segment,classcomp> & _rays);

private:
  enum PlaneAlign{
      AlignZ, AlignY, AlignX
  };

  void processMeshAlign(const PlaneAlign &palign, const PlaneAlign &salign);
  void processMeshPlane(const vec4d& plane, std::vector<Segment> *segments, std::vector<Triangle> *meshTris);
  void processMeshSegment(const Segment& segment, int index);
  bool intersectTriangleSegment(const Segment& segment, const Triangle& tri, Point *pnt);
  bool intersectPlaneTriangle(const vec4d& plane, const Triangle& tri, Segment *seg);
  bool intersectPlaneSegment(const vec4d& plane, const vec3d& p0, const vec3d& p1, vec3d *pt);
  vec4d makePlane(const vec3d& a, const vec3d& b, const vec3d& c);

private:
  DepthComplexity2D *_dc2d;

  // Input
  const TriMesh *_mesh;
  int _fboWidth;
  int _fboHeight;
  int _discretSteps;
  unsigned int _maximum;
  unsigned int _threshold;

  // State
  bool _computeHistogram;
  bool _computeMaximumRays;
  bool _computeGoodRays;

  // Output
  std::set<Segment,classcomp> _maximumRays;
  std::vector< std::set<Segment, classcomp> > _goodRays;
  std::vector<Segment> _usedPlanes;
  std::vector<unsigned long long> _histogram;
  std::vector<Point> _intersectionPoints;
  std::map<int, std::list<unsigned int>* > _intersectionTriList;
  //  std::vector<Segment> _intersectionSegments;

  friend int doInteractive(TriMesh& mesh);
  friend void drawRays();
};


#endif // DC_3D_H
