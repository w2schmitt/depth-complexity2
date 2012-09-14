#ifndef UTIL_H_
#define UTIL_H_

#include <vector>
#include <algorithm>

#include "vector.hpp"



struct Triangle {
    vec3d tca, a, na; 
    vec4d ca;
    vec3d tcb, b, nb;
    vec4d cb;
    vec3d tcc, c, nc;
    vec4d cc;
};


typedef vec3d Point;

struct Plane {
    Point a,b,c,d;
};


struct Segment {
  Segment(const Point &a_, const Point &b_):a(a_), b(b_), active(true) {}
  Segment() : active(true){}
  void sortPoints()
  {
	   double sega[3]={this->a.x, this->a.y, this->a.z};
       double segb[3]={this->b.x, this->b.y, this->b.z};
       if (!std::lexicographical_compare(sega,sega+3,segb,segb+3)){
		   vec3d aux = this->a;
		   this->a = this->b;
		   this->b = aux;
	   }
  }
  vec3d a, b;
  bool active;
};
bool operator<(const Segment &s1, const Segment &s2);

/*struct CuttingSegment: Segment {
  CuttingSegment() : intersect(0) {}
  CuttingSegment(const Segment &s, const int intersect_ = 0):intersect(intersect_) {}
  int intersect;
};*/

struct Ray{
    Segment s;
    unsigned int dc;
};


struct BoundingBox {
    vec3d min, max;
    bool empty;
    
    BoundingBox() : empty(true)
    {}
    
    BoundingBox& merge(const vec3d& p)
    {
        if (empty) {
            empty = false;
            min = max = p;
        } else {
            min.x = std::min(min.x, p.x);
            min.y = std::min(min.y, p.y);
            min.z = std::min(min.z, p.z);
            max.x = std::max(max.x, p.x);
            max.y = std::max(max.y, p.y);
            max.z = std::max(max.z, p.z);
        }
        return *this;
    }

    vec3d extents() const
    {
        if (empty)
            return vec3d(0,0,0);
        else
            return max - min;
    }

    vec3d center() const
    {
        return 0.5 * (min + max);
    }
};

struct TriMesh {
    std::vector<Triangle> faces;
    BoundingBox aabb;
};

struct classcomp {
	bool operator() (const Segment& lhs, const Segment& rhs) const
	{
		float seg1[6] = { lhs.a.x, lhs.a.y, lhs.a.z, 
						  lhs.b.x, lhs.b.y, lhs.b.z};
						  
		float seg2[6] = { rhs.a.x, rhs.a.y, rhs.a.z, 
						  rhs.b.x, rhs.b.y, rhs.b.z}; 
		
		return std::lexicographical_compare(seg1,seg1+6,seg2,seg2+6);	
	}
};

struct classcomp2 {
	bool operator() (const Ray& lhs, const Ray& rhs) const
	{
		float seg1[6] = { lhs.s.a.x, lhs.s.a.y, lhs.s.a.z, 
						  lhs.s.b.x, lhs.s.b.y, lhs.s.b.z};
						  
		float seg2[6] = { rhs.s.a.x, rhs.s.a.y, rhs.s.a.z, 
						  rhs.s.b.x, rhs.s.b.y, rhs.s.b.z}; 
		
		return std::lexicographical_compare(seg1,seg1+6,seg2,seg2+6);	
	}
};


TriMesh loadOFFMesh(std::istream& in);
TriMesh loadOBJMesh(std::istream& in);

// help functions
bool lineIntersection3D(const Segment &line1, const Segment &line2, double *t1, double *t2);
bool lineIntersection2D(const Segment &line1, const Segment &line2, double *t1, double *t2);
bool segmentIntersection2D(const Segment &seg1, const Segment &seg2, double *t1, double *t2);
bool segmentIntersection3D(const Segment &seg1, const Segment &seg2, double *t1, double *t2);


bool intersectPlaneSegment(const vec4d& plane, const vec3d& p0, const vec3d& p1, vec3d *pt);
vec4d makePlane(const vec3d& a, const vec3d& b, const vec3d& c);
bool intersectTriangleSegment(const Segment& segment, const Triangle& tri, Point *pnt);

#endif
