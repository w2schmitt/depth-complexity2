#ifndef DC_2D2_H
#define DC_2D2_H

#define GL_GLEXT_PROTOTYPES


#include <GL/glew.h>
#ifdef __APPLE__
#include <GLUT/glut.h>
#else
#include <GL/glut.h>
#include <algorithm>
#endif

#include "util.h"
#include "ShaderMgr.h"
//#include "clip_polygon.h"
#include <map>
#include <set>
#include <list>
#include <cstdlib>
#include <stdio.h>
#include <iomanip>

#include "CImg.h"
using namespace cimg_library;

// Compute the Depth complexity in 2D
//
// Usage:
//   int fboWidth = ..., fboHeight = ...;
//   DepthComplexity2D dc2d(fboWidth, fboHeight);
//
//   Segment from = ...;
//   Segment to = ...;
//   vector<Segment> segments = ...;
//   dc2d.process(from, to, segments);
//
class DepthComplexity2D {
public:
  DepthComplexity2D(const int fboWidth, const int fboHeight);
  ~DepthComplexity2D();

  void setComputeHistogram(bool computeHistogram) { this->_computeHistogram = computeHistogram; }
  void setComputeMaximumRays(bool computeMaximumRays) {this->_computeMaximumRays = computeMaximumRays; }
  void setComputeGoodRays(bool computeGoodRays) {this->_computeGoodRays = computeGoodRays; }
  void setThreshold(unsigned threshold) {this->_threshold = threshold; }
  //void setTex3dSize(vec3d size) { this->_texSize = size;}
  void setMeshBoundingbox(BoundingBox aabb) { this->_aabb = aabb; }
  //void setDCType(bool type) { _dc_type = type; }
  
  void cimg2Tex(unsigned int maxDC);
  // Compute the maximum depth complexity
  void process(
    const Segment &from, const Segment &to, const std::vector<Segment> &segments);

  // Copy stencil buffer to color buffer.
  // Colors are defined in a table.
  void copyStencilToColor();

  // Methods to obtain outputs
  unsigned int                     maximum() const           { return _maximum; }
  std::vector<unsigned long long>  histogram()               { return _histogram; }
  std::set<Segment,classcomp>      maximumRays()             { return _maximumRays; }
  std::set<Segment,classcomp>      goodRays(unsigned i)      { return _goodRays[i]; }
  std::vector<Point>               points()                  { return _points; }
  GLuint                           textureId() const         { return _cboTexId; }
  GLuint                           texture3DId() const       { return _texID; }
private:
  // Create framebuffer objects
  bool initFBO();

  // Primal: point -> Dual: Segment
  Segment computeDualSegmentFromPoint(const Point &p);
		
  void clipPolygon(const Point &p1, const Point &p2, const Point &p3, std::vector<Point> &polygon);
  void clipPolygon(const Point &p1, const Point &p2, const Point &p3, const Point &p4, std::vector<Point> &polygon);
  
  void bresenham_line_3D(vec3d &v1, vec3d &v2, unsigned int color);
  void drawPixel(vec3d pos, unsigned int color);
		
  // Go through the counter buffer to find the maximum value.
  unsigned int findMaxValueInCounterBuffer();
  void setShaderClearCounterBuffer();
  void setShaderCountDC();

  // Compute depth complexity using rays from one segment to the other.
  void findDepthComplexity2D();

  void findMaximumRaysAndHistogram();
  
  // functions for 3d texturing
  void updateTexture3D(Segment line, unsigned int dc);
  //void computeIntersectionPoints(std::list<Point> &pts);
  
private:
  //buffers 
  GLuint                                		_cboTexId;
  GLuint						_counterBuffId;
  GLuint                                		_fboId;
  GLuint                                		_rboId;
  
  // Shaders
  GLuint 						_shaderclearBuffer;
  GLuint                                                _shaderCountDC;

  // State
  bool                                  		_status;
  bool                                  		_computeHistogram;
  bool                                 			_computeMaximumRays;
  bool                                  		_computeGoodRays;
  //bool                                                  _dc_type;

  // Inputs
  int                                   		_fboWidth;
  int                                   		_fboHeight;
  Segment                               		_from;
  Segment                               		_to;
  const std::vector<Segment>*           		_segments;
  const std::vector<Triangle>*                          _meshTris;
  unsigned                              		_threshold;

  // Outputs
  unsigned int                                          _maximum;
  std::vector<unsigned long long>       		_histogram;
  std::set<Segment, classcomp>                          _maximumRays;
  std::vector< std::set<Segment,classcomp> >            _goodRays;
  std::vector<Point>                                    _points;
  
  CImg<float>                                           _tex3D;
  vec3d                                                 _texSize;
  BoundingBox                                           _aabb;
  GLuint                                                _texID;
  
};

#endif // DC_2D2_H
