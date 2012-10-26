#include "dc_3d.h"
#include "dc_2d.h"
#include "util.h"
#include "timer.h"

#include <cassert>
#include <iostream>
#include <iomanip>
#include <math.h>
#include <fstream>

//#define DEBUG_SAVE_HIST_EACH_FRAME

const double EPS = 1.0E-5;

template<class T>
T mix(const T& a, const T& b, double x) {
  const T& p = (1-x)*a;
  const T& q = x*b;
  return p + q;
}

DepthComplexity3D::DepthComplexity3D(int fboWidth, int fboHeight, int discretSteps):
  _fboWidth(fboWidth),
  _fboHeight(fboHeight),
  _discretSteps(discretSteps),
  _maximum(0),
  _computeHistogram(false),
  _computeMaximumRays(false),
  _computeGoodRays(false) {

  _goodRays.resize(1);
  _dc2d = new DepthComplexity2D(_fboWidth, _fboHeight);
}

DepthComplexity3D::~DepthComplexity3D() {
	std::map<int, std::list<unsigned int>* >::iterator it = _intersectionTriList.begin();
	for (; it != _intersectionTriList.end(); ++it)
		delete it->second;
}

void DepthComplexity3D::setComputeHistogram(bool computeHistogram) {
  this->_computeHistogram = computeHistogram;
  assert(_dc2d != 0);
  _dc2d->setComputeHistogram(computeHistogram);
}

void DepthComplexity3D::setComputeMaximumRays(bool computeMaximumRays) {
  this->_computeMaximumRays = computeMaximumRays;
  assert(_dc2d != 0);
  _dc2d->setComputeMaximumRays(computeMaximumRays);
}

void DepthComplexity3D::setComputeGoodRays(bool computeGoodRays) {
  this->_computeGoodRays = computeGoodRays;
  assert(_dc2d != 0);
  _dc2d->setComputeGoodRays(computeGoodRays);
}

void DepthComplexity3D::setThreshold(unsigned threshold) {
  this->_threshold = threshold;
  assert(_dc2d != 0);
  _dc2d->setThreshold(threshold);
}

void DepthComplexity3D::writeHistogram(std::ostream& out) {
  out << "Intersections Frequency\n";
  for (unsigned i=0; i< _histogram.size(); ++i)
    out << std::left << std::setw(14) << i << _histogram[i] << "\n";
}

void DepthComplexity3D::writeRays(std::ostream& out) {
  out << "OFF" << "\n";
  out << (_maximumRays.size()*2) << " "
      << _maximumRays.size() << " 0\n";

  std::set<Segment, classcomp>::const_iterator ite = _maximumRays.begin();
  std::set<Segment, classcomp>::const_iterator end = _maximumRays.end();
  for (; ite != end; ++ite) {
    out << ite->a.x << " " << ite->a.y << " " << ite->a.z << "\n"
        << ite->b.x << " " << ite->b.y << " " << ite->b.z << "\n";
  }

  for(unsigned i=0; i<2*_maximumRays.size(); i+=2)
    out << "2 " << i << " " << (i+1) << "\n";
}

void DepthComplexity3D::writeRays(std::ostream& out, const std::set<Segment,classcomp> & _rays) {
  out << "OFF" << "\n";
  out << (_rays.size()*2) << " "
      << _rays.size() << " 0\n";

  std::set<Segment,classcomp>::const_iterator ite = _rays.begin();
  std::set<Segment,classcomp>::const_iterator end = _rays.end();
  for (; ite != end; ++ite) {
    out << ite->a.x << " " << ite->a.y << " " << ite->a.z << "\n"
        << ite->b.x << " " << ite->b.y << " " << ite->b.z << "\n";
  }

  for(unsigned i=0; i<2*_rays.size(); i+=2)
    out << "2 " << i << " " << (i+1) << "\n";
}

void DepthComplexity3D::writeRaysSpherical(std::ostream& out, int k) {
  assert(_computeGoodRays);
  long long unsigned int total = 0;
  for(int i = 0 ; i <= k ; ++i) {
    int ind = maximum()-i;
    if (ind<0) break;
    const std::set<Segment, classcomp> & _rays = goodRays(ind);
    //std::cout << "GOOD_RAYS_SIZE["<<maximum()-i<<"] = "<<_rays.size() << std::endl;
    total += _rays.size();
  }
  out << total << std::endl;

  /* Simple spherical coordinates
  for(int i = 0 ; i <= k ; ++i) {
    int ind = maximum()-i;
    if (ind<0) break;
    const std::set<Segment,classcomp> & _rays = goodRays(ind);
    std::set<Segment,classcomp>::const_iterator ite = _rays.begin();
    std::set<Segment,classcomp>::const_iterator end = _rays.end();
    for (; ite != end; ++ite) {
      double a, b, c, d;
      double x1 = ite->a.x, x2 = ite->b.x, y1 = ite->a.y, y2 = ite->b.y, z1 = ite->a.z, z2 = ite->b.z;
      a = atan2(y1, x1);
      b = atan2(z1, sqrt(x1*x1 + y1*y1));
      c = atan2(y2, x2);
      d = atan2(z2, sqrt(x2*x2 + y2*y2));
      out << a << " " << b << " " << c << " " << d << " " << maximum()-i << "\n";
    }
  }
   */
  BoundingBox aabb = _mesh->aabb;
  Sphere sph;
  sph.center = aabb.center();
  sph.radius = aabb.extents().length()/2;
  
  // Coordinates of intersection with bounding sphere
  for (int i=0; i <= k; ++i){
      int ind = maximum()-i;
      if (ind<0) break;
      const std::set<Segment,classcomp>& _rays = goodRays(ind);
      std::set<Segment,classcomp>::const_iterator ite = _rays.begin();
      std::set<Segment,classcomp>::const_iterator end = _rays.end();
      for (; ite != end; ++ite) {
          vec3d f, s;
          if(!segmentSphereIntersection3D(*ite,sph,f,s)){
              continue;
          }
          f = cartesianToSpherical(f);
          s = cartesianToSpherical(s);
          out << f.y << " " << f.z << " " << s.y << " " << s.z << " " << maximum()-i << std::endl;
      }
  }  
}

void DepthComplexity3D::process(const TriMesh &mesh) {
  this->_mesh = &mesh;

  _usedPlanes.clear();
  _goodRays.clear();
  _goodRays.resize(1);
  _dc2d->setThreshold(_threshold);
  _maximum = 0;
  
  //std::cout << _fboWidth << " " << _fboHeight << " " << _discretSteps << " " << _maximum << " " << _threshold << std::endl;
  
  processMeshAlign(AlignZ, AlignX);
  processMeshAlign(AlignZ, AlignY);

  processMeshAlign(AlignY, AlignX);
  processMeshAlign(AlignY, AlignZ);

  processMeshAlign(AlignX, AlignY);
  processMeshAlign(AlignX, AlignZ);

}

// Call this varying palign and salign.
void DepthComplexity3D::processMeshAlign(const PlaneAlign &palign, const PlaneAlign &salign) {

  assert(palign != salign);
  BoundingBox aabb = _mesh->aabb;
  aabb.merge(aabb.min - aabb.extents()/10.0);
  aabb.merge(aabb.max + aabb.extents()/10.0);

  vec3d c0 = vec3d(aabb.min.x, aabb.min.y, aabb.min.z);
  vec3d c1 = vec3d(aabb.max.x, aabb.min.y, aabb.min.z);
  vec3d c2 = vec3d(aabb.min.x, aabb.max.y, aabb.min.z);
  vec3d c3 = vec3d(aabb.max.x, aabb.max.y, aabb.min.z);
  vec3d c4 = vec3d(aabb.min.x, aabb.min.y, aabb.max.z);
  vec3d c5 = vec3d(aabb.max.x, aabb.min.y, aabb.max.z);
  vec3d c6 = vec3d(aabb.min.x, aabb.max.y, aabb.max.z);
  vec3d c7 = vec3d(aabb.max.x, aabb.max.y, aabb.max.z);

  // generate all planes varying on z  
  const unsigned steps = _discretSteps;
  //unsigned az = 0;
  //unsigned bz = 0;
  for (unsigned az = 0; az<steps; ++az) {
    //double t = 3*az / (steps - 1.0) - 1; // [-1, 2]
    double t = az / (steps - 1.0);
    for (unsigned bz = 0; bz<steps; ++bz) {
      // double u = 3*bz / (steps - 1.0) - 1; // [-1, 2]
      double u = bz / (steps - 1.0);
      Segment sa, sb;
      switch (palign) {
      case AlignZ: {
        Point a = mix(c0, c4, t); // along Z
        Point b = mix(c7, c3, u); // along Z
        sa.a = sa.b = a;
        sb.a = sb.b = b;
		
        if (salign == AlignX) {
          // extend along X
          sa.b.x = aabb.max.x;
          sb.a.x = aabb.min.x;
        } else {
          // extend along Y
          sa.b.y = aabb.max.y;
          sb.a.y = aabb.min.y;
        }
        break;
      }
      case AlignY: {
        Point a = mix(c0, c2, t); // along Y
        Point b = mix(c7, c5, u); // along Y
        sa.a = sa.b = a;
        sb.a = sb.b = b;

        if (salign == AlignX) {
          // extend along X
          sa.b.x = aabb.max.x;
          sb.a.x = aabb.min.x;
        } else {
          // extend along Z
          sa.b.z = aabb.max.z;
          sb.a.z = aabb.min.z;
        }
        break;
      }
      case AlignX:  {
        Point a = mix(c0, c1, t);
        Point b = mix(c7, c6, u);
        sa.a = sa.b = a;
        sb.a = sb.b = b;

        if (salign == AlignY) {
          // extend along Y
          sa.b.y = aabb.max.y;
          sb.a.y = aabb.min.y;
        } else {
          // extend along Z
          sa.b.z = aabb.max.z;
          sb.a.z = aabb.min.z;
        }
        break;
      }
      }

      Segment saa, sbb;
      saa.a = sa.a;
      saa.b = sb.b;
      sbb.a = sa.b;
      sbb.b = sb.a;
      _usedPlanes.push_back(sa);
      _usedPlanes.push_back(saa);
      _usedPlanes.push_back(sb);
      _usedPlanes.push_back(sbb);

      vec4d plane = makePlane(sa.a, sa.b, sb.a);
      std::vector<Segment> segments;
      processMeshPlane(plane, &segments);

      // make the segments extra-long
      vec3d dsa = sa.b - sa.a; sa.a -= dsa; sa.b += dsa;
      vec3d dsb = sb.b - sb.a; sb.a -= dsb; sb.b += dsb;

	  
      _dc2d->process(sa, sb, segments);

      unsigned int tempMaximum = _dc2d->maximum();
						//if (tempMaximum==3)
								//std::cout << "values: " << az << " " << bz << std::endl;
      if (tempMaximum >= _maximum) {
        if (tempMaximum > _maximum) {
          _maximumRays.clear();
          _goodRays.resize(tempMaximum+1);
          _histogram.resize(tempMaximum+1);
          _intersectionPoints.clear();
          //          _intersectionSegments.clear();
        }
        _maximum = tempMaximum;
        std::set<Segment, classcomp> tempRays = _dc2d->maximumRays();
       // std::set<Segment, classcomp>::iterator it = tempRays.begin();

        // Testing rays and saving intersectin points.
        /*
        for (; it!=tempRays.end(); ++it) {
          for (unsigned s=0; s<segments.size(); ++s) {
            double t1, t2;
            if(segmentIntersection3D(*it, segments[s], &t1, &t2)) {
              _intersectionPoints.push_back(it->a + t1*(it->b-it->a));
            }
          }
        }
        */
        

//        _intersectionSegments.insert(_intersectionSegments.end(), segments.begin(), segments.end());
//        _intersectionPoints.insert(_intersectionPoints.end(), points.begin(), points.end());

        _maximumRays.insert(tempRays.begin(), tempRays.end());
        // Shouldn't the histogram be used without regard to the current tempMaximum? (changed it)
      }
      
      std::vector<unsigned long long> tempHist = _dc2d->histogram();
      for(unsigned i=0; i< tempHist.size(); ++i)
        _histogram[i] += tempHist[i];
      
#ifdef DEBUG_SAVE_HIST_EACH_FRAME
      char filename[100]; sprintf(filename,"hist/hist%d%d.txt",az,bz);
      std::ofstream fhist ( filename );
      writeHistogram(fhist);
      fhist.close();
#endif
      if(_computeGoodRays) {
        //std::cout << "size of goodRays: " << _goodRays.size() << " and _threshold = " << _threshold << std::endl;
        for(unsigned int i = _threshold ; i <= tempMaximum ; ++i) {
          //std::cout << "i = " << i << " and size(i) = " << _dc2d->goodRays(i).size() << std::endl;
          std::set<Segment,classcomp> tempRays = _dc2d->goodRays(i);
          _goodRays[i].insert(tempRays.begin(), tempRays.end());
        }
      }
    }
  }
    //int index=0;
    // check which mesh triangles intersect with each rays
    //std::set<Segment,classcomp>::iterator it = _maximumRays.begin();
    //for (; it!=_maximumRays.end(); ++it)
    //        processMeshSegment(*it,index++);
}

//INPUT: plane -> The normal vector of the plane which we will check overlaps;
//OUTPUT: segments -> A vector containing the segments of the mesh that intersect the plane.
void DepthComplexity3D::processMeshPlane(const vec4d& plane, std::vector<Segment> *segments) {
  assert(segments);

  for (unsigned i=0; i<_mesh->faces.size(); ++i) {
    Segment s;
    if (intersectPlaneTriangle(plane, _mesh->faces[i], &s))
      segments->push_back(s);
  }
}

void DepthComplexity3D::processMeshSegment(const Segment& segment, int rayIndex) {
	//assert(points);
  std::list<unsigned int> *triIndex = new std::list<unsigned int>;
	for (unsigned i=0; i<_mesh->faces.size(); ++i) {
		Point p;
		if (intersectTriangleSegment(segment, _mesh->faces[i], &p)){
			triIndex->push_back(i);
    }
	}
	_intersectionTriList.insert ( std::pair<int,std::list<unsigned int>* >(rayIndex,triIndex) );;
}

bool DepthComplexity3D::intersectTriangleSegment(const Segment& segment, const Triangle& tri, Point *pnt) {
	assert(pnt);

	if(!intersectPlaneSegment(makePlane(tri.a, tri.b, tri.c),segment.a,segment.b,pnt))
		return false;

	vec3d u = tri.b - tri.a;
	vec3d v = tri.c - tri.a;
	vec3d w =	*pnt - tri.a;

	double uu = dot(u,u);
	double uv = dot(u,v);
	double vv = dot(v,v);
	double wu = dot(w,u);
	double wv = dot(w,v);
	
	double den = uv*uv - uu*vv;

	double s = (uv*wv - vv*wu)/den;
	//std::cout << "s=" << s << std::endl;
	if(s<0. || s>1.)
		return false;
	double t = (uv*wu - uu*wv)/den;
	//std::cout << "t=" << t << std::endl;
	if(t<0. || s+t>1.)
		return false;
	
	return true;
}


bool DepthComplexity3D::intersectPlaneTriangle(const vec4d& plane, const Triangle& tri, Segment *seg) {
  assert(seg);

  double da, db, dc;
  da = dot(plane, vec4d(tri.a, 1));
  db = dot(plane, vec4d(tri.b, 1));
  dc = dot(plane, vec4d(tri.c, 1));

  if (fabs(da)< EPS) da=0.f;
  if (fabs(db)< EPS) db=0.f;
  if (fabs(dc)< EPS) dc=0.f;

  if (da == 0 && db == 0 && dc == 0) // if the triangle is inside the plane, there's no intersection
    return false;

  vec3d *p = &seg->a;

  if (std::signbit(da) != std::signbit(db)) {
    intersectPlaneSegment(plane, tri.a, tri.b, p);
    p = &seg->b;
  }

  if (std::signbit(db) != std::signbit(dc)) {
    intersectPlaneSegment(plane, tri.b, tri.c, p);
    p = &seg->b;
  }

  if (std::signbit(dc) != std::signbit(da)) {
    intersectPlaneSegment(plane, tri.c, tri.a, p);
    p = &seg->b;
  }
  //if (p == &seg->b)
    //printf("da = %f, db = %f dc= %f\npx = %f, py = %f, pz = %f ------- (%s)\n\n", da,db,dc,p->x, p->y, p->z, (p== &seg->b)? "true":"false");
  return p == &(seg->b);
}

bool DepthComplexity3D::intersectPlaneSegment(const vec4d& plane, const vec3d& p0, const vec3d& p1, vec3d *pt) {
  double num = -plane.w - dot(plane.xyz(), p0);
  double den = dot(plane.xyz(), p1 - p0);
  double r = num / den;
  *pt = mix(p0, p1, r);
  if (0 <= r && r <= 1)
      return true;
  return false;
}
// GIVEN 3 POINTS --> RETURN A 4D VECTOR: 
vec4d DepthComplexity3D::makePlane(const vec3d& a, const vec3d& b, const vec3d& c) {
    vec3d normal = cross(b-a, c-a);
    normal.normalize();
    double d = dot(a, normal);
    return vec4d(normal, -d);
}
