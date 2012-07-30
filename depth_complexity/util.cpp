#include "util.h"

#include <iostream>
#include <cstdio>
#include <boost/tokenizer.hpp>

const double EPS = 1.0E-5;

template<class T>
T mix(const T& a, const T& b, double x) {
  const T& p = (1-x)*a;
  const T& q = x*b;
  return p + q;
}


TriMesh
loadOFFMesh(std::istream& in){
  std::clog << "Loading OFF file" << std::endl;

  std::string head;
  in >> head;
  if (head != "OFF")
    throw "Does not start with OFF!";

  int nverts, nfaces, nedges;

  if (!(in >> nverts >> nfaces >> nedges))
    throw "Could not read number of vertices, faces, edges";

  TriMesh mesh;
  mesh.faces.reserve(nfaces);

  std::vector<vec3d> vertices(nverts);
  for (int i=0; i<nverts; ++i) {
    in >> vertices[i].x >> vertices[i].y >> vertices[i].z;
    mesh.aabb.merge(vertices[i]); // build AABB
  }

  for (int i=0; i<nfaces; ++i) {
    int sz, a, b, c;
    Triangle t;
    in >> sz >> a >> b >> c;
    t.a = vertices.at(a);
    t.b = vertices.at(b);
    t.c = vertices.at(c);
    mesh.faces.push_back(t);

    for (int j=3; j<sz; ++j) {
      //std::clog << "triangulating face " << i << std::endl;

      b = c;
      in >> c;
      t.b = vertices.at(b);
      t.c = vertices.at(c);
      mesh.faces.push_back(t);
    }
  }

  std::clog << "Loaded " << mesh.faces.size() << " faces" << std::endl;

  // now build the normals
  for (unsigned i=0; i<mesh.faces.size(); ++i) {
    Triangle& t = mesh.faces[i];
    const vec3d ab = t.b - t.a;
    const vec3d ac = t.c - t.a;
    vec3d n = cross(ab, ac);
    n.normalize();
    t.na = t.nb = t.nc = n;
  }

  return mesh;
}


TriMesh
loadOBJMesh(std::istream& in) {
  std::clog << "Loading OBJ file" << std::endl;

  std::string line;

  std::vector<vec3d> vertices;
  std::vector<vec3d> normals;

  TriMesh mesh;

  for(;;) {
    if (!getline(in, line))
      break;
    if (line == "")
      continue;

    switch (line[0]) {
    case 'v':
      if (line[1] == ' ') {
        vec3d p;
        if (sscanf(line.c_str()+2, "%lf %lf %lf", &p.x, &p.y, &p.z) != 3)
          throw std::string("Error reading vertex at line " + line);
	
        vertices.push_back(p);
        mesh.aabb.merge(p);
      } else if (line[1] == 'n') {
        vec3d n;
        if (sscanf(line.c_str()+3, "%lf %lf %lf", &n.x, &n.y, &n.z) != 3)
          throw std::string("Error reading normal at line " + line);
        normals.push_back(n);
      }
      break;

    case 'f':
    {
      std::string content = line.substr(2);
      boost::char_separator<char> sep(" ");
      boost::tokenizer<boost::char_separator<char> > tok(content, sep);
      int a=-2, b=-2, c=-2, na=-2, nb=-2, nc=-2;
      for (boost::tokenizer<boost::char_separator<char> >::iterator i=tok.begin(); i!=tok.end(); ++i) {
        int v=-1, n=-1;
        if (sscanf(i->c_str(), "%d//%d", &v, &n) != 2
            &&
            sscanf(i->c_str(), "%d/%*d/%d", &v, &n) != 2
            &&
            sscanf(i->c_str(), "%d/%*d", &v) != 1 )
          throw std::string("Error reading face at line " + line);

        if (a == -2) {
          // first vertex
          a = v;
          na = n;
        } else if (b == -2) {
          // second vertex
          b = v;
          nb = n;
        } else {
          // third or later vertices, formed a complete face
          c = v;
          nc = n;

          Triangle t;
          t.a = vertices.at(a-1);
          t.b = vertices.at(b-1);
          t.c = vertices.at(c-1);

          // if we don't have a normal, make one
          if (na <= 0) {
            vec3d normal = cross(t.b - t.a, t.c - t.a);
            normal.normalize();
            std::clog << "created normal for face " << mesh.faces.size() << std::endl;
            normals.push_back(normal);
            na = nb = nc = normals.size();
          }
          t.na = normals.at(na-1);
          t.nb = normals.at(nb-1);
          t.nc = normals.at(nc-1);

          t.ca = vec4d(0.057f, 0.25f, 0.42f, 0.35f);
          t.cb = vec4d(0.057f, 0.25f, 0.42f, 0.35f);
          t.cc = vec4d(0.057f, 0.25f, 0.42f, 0.35f);

          mesh.faces.push_back(t);			

          b = c;
          nb = nc;
        }
      }
      break;
    }
    default: {
      std::clog << "ignored line: " << line << std::endl;
    }
    }
  }
  
        // update Bounding Box 
    mesh.aabb.merge(mesh.aabb.min - mesh.aabb.extents()/10.0);
    mesh.aabb.merge(mesh.aabb.max + mesh.aabb.extents()/10.0);    
  
    vec3d desloc = mesh.aabb.min;
    vec3d bblen = mesh.aabb.extents();
    vec3d pos;
  


  // set up 3d texture coordinates based in the position of the triangle
  for (unsigned int i=0; i < mesh.faces.size(); ++i){
      Triangle &t = mesh.faces[i];
      
      //vertex 1
      pos = t.a - desloc;
      pos.x = pos.x/bblen.x; 
      pos.y = pos.y/bblen.y;
      pos.z = pos.z/bblen.z;
      t.tca = pos;
      
      //vertex 2
      pos = t.b - desloc;
      pos.x = pos.x/bblen.x; 
      pos.y = pos.y/bblen.y;
      pos.z = pos.z/bblen.z;
      t.tcb = pos;
      
      //vertex 3
      pos = t.c - desloc;
      pos.x = pos.x/bblen.x; 
      pos.y = pos.y/bblen.y;
      pos.z = pos.z/bblen.z;
      t.tcc = pos;    
      

      //std::cout << t.tca << std::endl;
      //std::cout << t.tcb << std::endl;
      //std::cout << t.tcc << std::endl;
  }

    //std::cout << mesh.aabb.max << std::endl;
    //std::cout << mesh.aabb.min << std::endl;
  
  std::clog << "Loaded " << mesh.faces.size() << " faces" << std::endl;

  return mesh;
}

#define perp(u,v)  ((u).x * (v).y - (u).y * (v).x)  // perp product (2D)
// lineIntersect in 2D
bool lineIntersection2D(const Segment &line1, const Segment &line2, double *t1, double *t2) {
  vec3d u = line1.b - line1.a;
  vec3d v = line2.b - line2.a;
  vec3d w = line1.a - line2.a;
  // just X and Y values are used
  double D = perp(u,v);

  // test if they are parallel (includes either being a point)
  if (fabs(D) < EPS) // S1 and S2 are parallel
    return false;

  // the segments are skew and may intersect in a point
  // get the intersect parameter for S1
  *t1 = perp(v,w) / D;
  *t2 = perp(u,w) / D;

  return true;
}

// line intesection in 3D
// line1 = line1.a + t1*(line1.b-line1.a)
// line2 = line2.a + t2*(line2.b-line2.a)
bool lineIntersection3D(const Segment &line1, const Segment &line2, double *t1, double *t2) {
  vec3d v = line2.b - line2.a;
  vec3d u = line1.b - line1.a;
  vec3d w = line1.a - line2.a;

  // test if segments are points
  if (v.length() < EPS or u.length() < EPS)
    return false;

  double a = dot(u, u);
  double b = dot(v, u);
  double c = dot(v, v);
  double d = dot(w, u);
  double e = dot(w, v);

  double denom = a * c - b * b;
  if (fabs(denom) < EPS)
    return false;
  double numer = e * b - d * c;

  *t1 = numer / denom;
  *t2 = (e + b * (*t1)) / c;

  return true;
}

bool segmentIntersection3D(const Segment &seg1, const Segment &seg2, double *t1, double *t2) {
  if (lineIntersection3D(seg1, seg2, t1, t2) and
      ((0.0 < *t1 && *t1 < 1.0) and (0.0 < *t2 && *t2 < 1.0)))
    return true;
  return false;
}

// segmentIntersection2D(): the intersection of 2 finite 2D segments
bool segmentIntersection2D(const Segment &seg1, const Segment &seg2, double *t1, double *t2){
  if (lineIntersection2D(seg1, seg2, t1, t2) and
      (0.0 < *t1 && *t1 < 1.0) and
      (0.0 < *t2 && *t2 < 1.0))
    return true;
  return false;
}

bool intersectPlaneSegment(const vec4d& plane, const vec3d& p0, const vec3d& p1, vec3d *pt) {
  double num = -plane.w - dot(plane.xyz(), p0);
  double den = dot(plane.xyz(), p1 - p0);
  double r = num / den;
  *pt = mix(p0, p1, r);
  if (0 <= r && r <= 1)
      return true;
  return false;
}


// GIVEN 3 POINTS --> RETURN A 4D VECTOR NORMAL (PLANE EQUATION): 
vec4d makePlane(const vec3d& a, const vec3d& b, const vec3d& c) {
    vec3d normal = cross(b-a, c-a);
    normal.normalize();
    double d = dot(a, normal);
    return vec4d(normal, -d);
}


bool intersectTriangleSegment(const Segment& segment, const Triangle& tri, Point *pnt) {
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
