#include "dc_3d_random.h"
#include "util.h"
#include "timer.h"

#include <cassert>
#include <iostream>
#include <iomanip>
#include <cstdlib>
#include <ctime>
#include <cmath>

const double EPS = 1.0E-5;
double const PI = 4*atan(1); 

template<class T>
T mix(const T& a, const T& b, double x) {
	const T& p = (1-x)*a;
	const T& q = x*b;
	return p + q;
}

void RDepthComplexity3D::writeRays(std::ostream& out, const std::set<Segment,classcomp> & _rays, int dc) {
  //out << "OFF" << "\n";
  //out << (_maximumRays.size()*2) << " "
  //    << _maximumRays.size() << " 0\n";

  std::set<Segment, classcomp>::const_iterator ite = _rays.begin();
  std::set<Segment, classcomp>::const_iterator end = _rays.end();
  for (; ite != end; ++ite) {
    out << ite->a.x << "," << ite->a.y << "," << ite->a.z << "," << ite->b.x << "," << ite->b.y << "," << ite->b.z << "," << dc << "\n";
  }

  //for(unsigned i=0; i<2*_maximumRays.size(); i+=2)
  //  out << "2 " << i << " " << (i+1) << "\n";
}

void RDepthComplexity3D::writeMaximumRays(std::ostream& out){
    writeRays(out, _maximumRays, _maximum);
}

void RDepthComplexity3D::writeGoodRays(std::ostream& out){
    for(unsigned i = _threshold; i <= _maximum; ++i) {
        writeRays(out, _goodRays[i], i);
    }
}

RDepthComplexity3D::RDepthComplexity3D(int fboWidth, int fboHeight, int discretSteps):
	_fboWidth(fboWidth),
	_fboHeight(fboHeight),
	_discretSteps(discretSteps),
	_maximum(0),
	_computeHistogram(false),
	_computeMaximumRays(false),
	_computeGoodRays(false),
	_computeRaysFromFile(false){
	
	_goodRays.resize(1);
}

RDepthComplexity3D::RDepthComplexity3D(int fboWidth, int fboHeight, int discretSteps, const char* filenameRays):
	_fboWidth(fboWidth),
	_fboHeight(fboHeight),
	_discretSteps(discretSteps),
	_maximum(0),
	_computeHistogram(false),
	_computeMaximumRays(false),
	_computeGoodRays(false) {
	
	_goodRays.resize(1);

	if (strcmp(filenameRays, "")!=0){
		try{
			std::ifstream file(filenameRays);
			_computeRaysFromFile = readRaysFromFile(file);
			
		}
		catch (const char* msg)	{
				std::cerr << "Failed: " << msg << std::endl;
		}
		catch (std::string msg) {
				std::cerr << "Failed: " << msg << std::endl;
		}
	}
}

bool RDepthComplexity3D::readRaysFromFile(std::istream& in){
	std::string line;
	unsigned i=0;
	unsigned limit=50;
	while (1){
		Segment s;
		if (i>limit) break;
		if (getline(in, line)){
			//std::cout << line.c_str() << std::endl;
			if (sscanf(line.c_str(), "%lf %lf %lf %lf %lf %lf\n", &s.a.x, &s.a.y, &s.a.z, &s.b.x, &s.b.y, &s.b.z)!=6){
				std::clog << "[ERROR] Wrong file rays format.\n"; 
				return false;
			} 
			else {				
				i++;
				_raysFromFile.push_back(s);
			}				
		} 
		else{
			std::clog << "Loaded " << _raysFromFile.size() << " rays from file!" << std::endl;
			return true;
		}
	}
	std::clog << "Loaded " << _raysFromFile.size() << " rays from file!" << std::endl;
	return true;
}

RDepthComplexity3D::~RDepthComplexity3D() {}

void RDepthComplexity3D::setComputeHistogram(bool computeHistogram) {
	this->_computeHistogram = computeHistogram;
}

void RDepthComplexity3D::setComputeMaximumRays(bool computeMaximumRays) {
	this->_computeMaximumRays = computeMaximumRays;
}

void RDepthComplexity3D::setComputeGoodRays(bool computeGoodRays) {
	this->_computeGoodRays = computeGoodRays;
}

void RDepthComplexity3D::writeRaysSpherical(std::ostream& out, int k) {
  assert(_computeGoodRays);
  long long unsigned int total = 0;
  for(int i = 0 ; i <= k ; ++i) {
    int ind = maximum()-i;
    if (ind<0) break;
    const std::set<Segment,classcomp> & _rays = goodRays(maximum()-i);
    total += _rays.size();
  }
  out << total << "\n";

  /*
  for(int i = 0 ; i <= k ; ++i) {
    int ind = maximum()-i;
    if (ind<0) break;
    const std::set<Segment,classcomp> & _rays = goodRays(maximum()-i);
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

void RDepthComplexity3D::setThreshold(unsigned threshold) {
	this->_threshold = threshold;
}

void RDepthComplexity3D::writeHistogram(std::ostream& out) {
	out << "Intersections Frequency\n";
	for (unsigned i=0; i< _histogram.size(); ++i)
		out << std::left << std::setw(14) << i << _histogram[i] << "\n";
}

void RDepthComplexity3D::writeRays(std::ostream& out) {
	out << "OFF" << "\n";
	out << (_maximumRays.size()*2) << " "
			<< _maximumRays.size() << " 0\n";

	std::set<Segment,classcomp>::const_iterator ite = _maximumRays.begin();
	std::set<Segment,classcomp>::const_iterator end = _maximumRays.end();
	for (; ite != end; ++ite) {
		out << ite->a.x << " " << ite->a.y << " " << ite->a.z << "\n"
				<< ite->b.x << " " << ite->b.y << " " << ite->b.z << "\n";
	}

	for(unsigned i=0; i<2*_maximumRays.size(); i+=2)
		out << "2 " << i << " " << (i+1) << "\n";
}

void RDepthComplexity3D::writeRays(std::ostream& out, const std::set<Segment,classcomp> & _rays) {
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

double uniformRandom() { return rand()/double(RAND_MAX); }

#define CONSTANT_FACTOR 500
void RDepthComplexity3D::process(const TriMesh &mesh) {
	 this->_mesh = &mesh;

	_usedPlanes.clear();
	_goodRays.clear();
	_goodRays.resize(1);
	_histogram.clear();
	_histogram.resize(1);
	_maximum = 0;
	
	/* test */
	Point p;
	Segment seg;
	seg.a = Point(0,0,0);
	seg.b = Point(0,5,0);
	Triangle tri;
	tri.a = Point(-2,3,-2);
	tri.b = Point(2,3,-2);
	tri.c = Point(-1,3,5);
	assert(intersectTriangleSegment(seg,tri,&p));
	
	//std::cout << _fboWidth << " " << _fboHeight << " " << _discretSteps << " " << _maximum << " " << _threshold << std::endl;

	BoundingBox aabb = _mesh->aabb;
	aabb.merge(aabb.min - aabb.extents()/10.0);
	aabb.merge(aabb.max + aabb.extents()/10.0);
	
	double r = aabb.extents().length()/2.0;
	vec3d center = aabb.center();
	
	srand(time(0));

	const unsigned steps = _discretSteps;
	const unsigned Nrays = (_computeRaysFromFile? _raysFromFile.size() : steps*steps*CONSTANT_FACTOR);

	for (unsigned tn = 0; tn < Nrays ; ++tn) {
		Segment randomSegment;

		if (!_computeRaysFromFile){		
			double theta =	2.*uniformRandom()*PI;
			double cosPhi = 2.*uniformRandom() - 1.;
			double sinPhi = sqrt(1. - cosPhi*cosPhi);
			vec3d randomVector1(cos(theta)*sinPhi*r,sin(theta)*sinPhi*r,cosPhi*r);
			
			theta =	2.*uniformRandom()*PI;
			cosPhi = 2.*uniformRandom() - 1.;
			sinPhi = sqrt(1. - cosPhi*cosPhi);
			vec3d randomVector2(cos(theta)*sinPhi*r,sin(theta)*sinPhi*r,cosPhi*r);
			
			randomSegment = Segment(center+randomVector1,center+randomVector2);
		}
		else {
			randomSegment = _raysFromFile[tn];
		}

		std::vector<Point> points;
		processMeshSegment(randomSegment, &points);

		unsigned tempMaximum = points.size();
                if (tempMaximum >= _maximum){
                    _maximum = tempMaximum;
                    _goodRays.resize(tempMaximum+1);
                    _histogram.resize(tempMaximum+1);
                }
                
		if (tempMaximum >= _maximum) {
			if (tempMaximum > _maximum) {
				_maximumRays.clear();
				_goodRays.resize(tempMaximum+1);
				_histogram.resize(tempMaximum+1);
				_intersectionPoints.clear();
			}
			_maximum = tempMaximum;

			//_intersectionPoints.insert(_intersectionPoints.end(), points.begin(), points.end());

                        //if (_maximumRays.size() < 50)
                                _maximumRays.insert(randomSegment);
			// Shouldn't the histogram be used without regard to the current tempMaximum? (changed it)
		}
		
		if(_computeHistogram)
			++_histogram[points.size()];
		
		if(_computeGoodRays && points.size() >= _threshold)// && _goodRays[points.size()].size() < 50)
			_goodRays[points.size()].insert(randomSegment);
	}
}

void RDepthComplexity3D::processMeshSegment(const Segment& segment, std::vector<Point> *points) {
	assert(points);

	for (unsigned i=0; i<_mesh->faces.size(); ++i) {
		Point p;
		if (intersectTriangleSegment(segment, _mesh->faces[i], &p))
			points->push_back(p);
	}
}

bool RDepthComplexity3D::intersectTriangleSegment(const Segment& segment, const Triangle& tri, Point *pnt) {
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

bool RDepthComplexity3D::intersectPlaneSegment(const vec4d& plane, const vec3d& p0, const vec3d& p1, vec3d *pt) {
	double num = -plane.w - dot(plane.xyz(), p0);
	double den = dot(plane.xyz(), p1 - p0);
	if (fabs(den) < EPS)
		return false;
	double r = num / den;
	//std::cout << "r=" << r << "\n";
	*pt = mix(p0, p1, r);
	if (0 <= r && r <= 1)
			return true;
	return false;
}

vec4d RDepthComplexity3D::makePlane(const vec3d& a, const vec3d& b, const vec3d& c) {
		vec3d normal = cross(b-a, c-a);
		normal.normalize();
		double d = dot(a, normal);
		return vec4d(normal, -d);
}
