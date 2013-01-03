/* 
 * File:   RFDepthComplexity3D.h
 * Author: wagner
 *
 * Created on 4 de Dezembro de 2012, 13:51
 */

#ifndef RFDEPTHCOMPLEXITY3D_H
#define	RFDEPTHCOMPLEXITY3D_H

#include "util.h"
#include "ShaderMgr.h"
#include "CImg.h"

#include <map>
#include <set>
#include <vector>
#include <fstream>
#include <sstream>
#include <cstring>
#include <time.h>

#include <GL/glew.h>

using namespace cimg_library;

class RFDepthComplexity3D {
public:
    RFDepthComplexity3D(int fboWidth, int fboHeight, int discretSteps);
        //RDepthComplexity3D(int fboWidith, int fboHeight, int discretSteps, const char* filenameRays);
        
    void process(const TriMesh &mesh);
    
    RFDepthComplexity3D();
    RFDepthComplexity3D(const RFDepthComplexity3D& orig);
    virtual ~RFDepthComplexity3D();
    
  // sets
  void setComputeHistogram(bool computeHistogram);
  void setComputeMaximumRays(bool computeMaximumRays);
  void setComputeGoodRays(bool computeGoodRays);
  void setThreshold(unsigned threshold);
  
  // gets
  unsigned maximum() const { return _maximum; }
  const std::set<Segment,classcomp> &maximumRays() const { return _maximumRays; }
  const std::set<Segment,classcomp> &goodRays(int intersect) const { return _goodRays[intersect]; }
  const std::vector<Segment> &usedPlanes() const { return _usedPlanes; }
  const std::vector<Point> &intersectionPoints() const { return _intersectionPoints; }
  unsigned getThreshold() { return _threshold; }
  
  // 
  void writeHistogram(std::ostream& out);
  void writeRays(std::ostream& out);
  void writeRays(std::ostream& out, const std::set<Segment,classcomp> & _rays);
  void writeRaysSpherical(std::ostream& out, int k);
  void writeMaximumRays(std::ostream& out);
  void writeGoodRays(std::ostream& out);
  void writeRays(std::ostream& out, const std::set<Segment,classcomp> & _rays, int dc);
  
private:
  enum PlaneAlign{
      AlignZ, AlignY, AlignX
  };

  void processMeshSegment(const Segment& segment, std::vector<Point> *points);
  bool intersectTriangleSegment(const Segment& segment, const Triangle& tri, Point *pnt);
  bool intersectPlaneSegment(const vec4d& plane, const vec3d& p0, const vec3d& p1, vec3d *pt);
  vec4d makePlane(const vec3d& a, const vec3d& b, const vec3d& c);

  void initTextureCounter();
  unsigned int renderScene(vec3d point);
  void setShaderClearCounterBuffer();
  void setShaderCountDC();
  unsigned int findMaxValueInCounterBuffer();
  void erodeTriangle(vec3d &v1, vec3d &v2, vec3d &v3);
  
  CImgDisplay dualDisplay;
  
//bool readRaysFromFile(std::istream& in);
  
private:
  //buffers 
  GLuint						_counterBuffId;
  GLuint                                		_rboId;
  
  // Shaders
  GLuint 						_shaderclearBuffer;
  GLuint                                                _shaderCountDC;

  // Input
  const TriMesh *_mesh;
  std::vector<Triangle> _sorted_faces;
  int _fboWidth;
  int _fboHeight;
  int _discretSteps;
  unsigned _maximum;
  unsigned _threshold;
//std::vector<Segment> _raysFromFile;

  // State
  bool _computeHistogram;
  bool _computeMaximumRays;
  bool _computeGoodRays;
  bool _computeRaysFromFile;

  // Output
  std::set<Segment,classcomp> _maximumRays;
  std::vector< std::set<Segment,classcomp> > _goodRays;
  std::vector<Segment> _usedPlanes;
  std::vector<unsigned long long> _histogram;
  std::vector<Point> _intersectionPoints;
  //  std::vector<Segment> _intersectionSegments;

  friend int doInteractive(TriMesh& mesh);
  friend void drawRays();

};

#endif	/* RFDEPTHCOMPLEXITY3D_H */

