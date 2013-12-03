#ifndef UTIL_H_
#define UTIL_H_

#include <boost/config.hpp>

#include <vector>
#include <utility>
#include <algorithm>
#include <iomanip>

#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/connected_components.hpp>
#include <boost/graph/dijkstra_shortest_paths.hpp>

using namespace boost;

#include "vector.hpp"

// do not change this structure parameters
struct Triangle {
    vec3d tca, a, na;
    vec4d ca;
    double ia;
    vec3d tcb, b, nb;
    vec4d cb;
    double ib;
    vec3d tcc, c, nc;
    vec4d cc;
    double ic;
};




struct Sphere {
  vec3d center;
  double radius;
};

typedef vec3d Point;

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

struct GeoTest{
  Segment segment;
  vec3d src;
  vec3d dst;
  vec3d closeSrc;
  vec3d closeDst;
  std::vector<vec3d> path;
  double cost;
  double dist;
};
/*struct CuttingSegment: Segment {
  CuttingSegment() : intersect(0) {}
  CuttingSegment(const Segment &s, const int intersect_ = 0):intersect(intersect_) {}
  int intersect;
};*/

struct MyGraph {
    std::vector<std::vector< std::pair<int,float> > > edges;
};

// information stored in vertices
struct VertexInformation {
  unsigned component;
};

// information stored in edges
struct EdgeInformation {
  double weight;
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
    std::vector<vec3d> vertices;
    std::vector<Triangle> faces;
    BoundingBox aabb;
};

struct classcomp {
	bool operator() (const Segment& lhs, const Segment& rhs) const
	{
		double seg1[6] = { lhs.a.x, lhs.a.y, lhs.a.z, 
						  lhs.b.x, lhs.b.y, lhs.b.z};
						  
		double seg2[6] = { rhs.a.x, rhs.a.y, rhs.a.z, 
						  rhs.b.x, rhs.b.y, rhs.b.z}; 
		
		return std::lexicographical_compare(seg1,seg1+6,seg2,seg2+6);	
	}
};

TriMesh loadOFFMesh(std::istream& in);
TriMesh loadOBJMesh(std::istream& in);
TriMesh loadOBJMesh(std::istream& in, vec3d rotation);

// help functions
bool lineIntersection3D(const Segment &line1, const Segment &line2, double *t1, double *t2);
bool lineIntersection2D(const Segment &line1, const Segment &line2, double *t1, double *t2);
bool segmentIntersection2D(const Segment &seg1, const Segment &seg2, double *t1, double *t2);
bool segmentIntersection3D(const Segment &seg1, const Segment &seg2, double *t1, double *t2);
bool segmentSphereIntersection3D(const Segment &seg, const Sphere sph, vec3d & out0, vec3d & out1);
vec3d cartesianToSpherical(vec3d point);
void createBoostGraph(TriMesh& mesh);
double boost_dijkstra(int src, int dst, const TriMesh* mesh, GeoTest &debug);
MyGraph *createGraphFromModel(TriMesh& mesh);

#endif
