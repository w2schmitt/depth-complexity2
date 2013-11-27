#include "util.h"

#include <iostream>
#include <cstdio>
#include <boost/tokenizer.hpp>

// LINEAR ALGEBRA LIBRARY
#include "Eigen/Dense"

using Eigen::Affine3f;
using Eigen::Matrix3f;
using Eigen::Matrix4f;
using Eigen::Vector3f;
using Eigen::Vector4f;
using Eigen::AngleAxisf;

const double EPS = 1.0E-5;

TriMesh
loadOFFMesh(std::istream& in){
  std::clog << std::endl << "=== MODEL INFO ===" << std::endl;

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
    
    t.ca = vec4d(0.657f, 0.75f, 0.32f, 0.85f);
    t.cb = vec4d(0.657f, 0.75f, 0.32f, 0.85f);
    t.cc = vec4d(0.657f, 0.75f, 0.32f, 0.85f);
          
    //t.ca = vec4d(0.057f, 0.25f, 0.42f, 0.35f);
    //t.cb = vec4d(0.057f, 0.25f, 0.42f, 0.35f);
    //t.cc = vec4d(0.057f, 0.25f, 0.42f, 0.35f);
    if (i > 0){
      Triangle ant = mesh.faces.back();
      if (ant.a != t.a || ant.b != t.b || ant.c != t.c){
        mesh.faces.push_back(t);
      }
    } else {
      mesh.faces.push_back(t);  
    }

    for (int j=3; j<sz; ++j) {
      //std::clog << "triangulating face " << i << std::endl;

      b = c;
      in >> c;
      t.b = vertices.at(b);
      t.c = vertices.at(c);
      mesh.faces.push_back(t);
    }
  }
  mesh.vertices = vertices;
  std::clog << "  --> TYPE: " << "OFF" << std::endl;
  std::clog << "  --> VERTICES: " << mesh.vertices.size() << std::endl;
  std::clog << "  --> FACES: " << mesh.faces.size() << std::endl;
  //std::clog << "Loaded " << mesh.faces.size() << " faces" << std::endl;

  // now build the normals
  for (unsigned i=0; i<mesh.faces.size(); ++i) {
    Triangle& t = mesh.faces[i];
    const vec3d ab = t.b - t.a;
    const vec3d ac = t.c - t.a;
    vec3d n = cross(ab, ac);
    n.normalize();
    t.na = t.nb = t.nc = n;
  }
  
  
  //std::cout << "mesh vertices: " << mesh.vertices.size() << std::endl;
  
   // e Bounding Box 
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

  return mesh;
}


TriMesh
loadOBJMesh(std::istream& in) {
  std::clog << "=== MODEL INFO ===" << std::endl;

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

										//printf("ta = (%f,%f) - tb = (%f,%f) - tc = (%f,%f)\n", t.a.x, t.a.y, t.b.x, t.b.y, t.c.x, t.c.y); 

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

          //t.ca = vec4d(0.057f, 0.25f, 0.42f, 0.35f);
          //t.cb = vec4d(0.057f, 0.25f, 0.42f, 0.35f);
          //t.cc = vec4d(0.057f, 0.25f, 0.42f, 0.35f);
          
          t.ca = vec4d(0.457f, 0.45f, 0.42f, 0.85f);
          t.cb = vec4d(0.457f, 0.45f, 0.42f, 0.85f);
          t.cc = vec4d(0.457f, 0.45f, 0.42f, 0.85f);

          mesh.faces.push_back(t);			

          b = c;
          nb = nc;
        }
      }
      break;
    }
    default: {
      //std::clog << "ignored line: " << line << std::endl;
    }
    }
    mesh.vertices = vertices;
  }
  
  
  
    // e Bounding Box 
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
  std::clog << "  --> TYPE: OBJ " << std::endl;
  std::clog << "  --> VERTICES: " << mesh.vertices.size() << std::endl;
  std::clog << "  --> FACES: " << mesh.faces.size() << std::endl;

  //std::clog << "Loaded " << mesh.faces.size() << " faces" << std::endl;

  return mesh;
}


TriMesh
loadOBJMesh(std::istream& in, vec3d rotation) {
  
  std::clog << "Loading OBJ file" << std::endl;
  
  Affine3f rot = Affine3f::Identity();
  
  if (rotation.x!=0 || rotation.y !=0 || rotation.z != 0){
      Vector3f x,y,z;
      x << 1,0,0;       y << 0,1,0;     z << 0,0,1;      
      rot.rotate(AngleAxisf(rotation.x, x));   
      rot.rotate(AngleAxisf(rotation.y, y)); 
      rot.rotate(AngleAxisf(rotation.z, z)); 
  }
   

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
								//printf("px = (%f) - py = (%f) - pz = (%f)\n", p.x, p.y, p.z); 
        //Vector3f ptmp; ptmp << p.x, p.y, p.z;
        //ptmp = rot*ptmp;
        //p.x = ptmp(0); p.y=ptmp(1); p.z=ptmp(2);
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

	//printf("ta = (%f,%f) - tb = (%f,%f) - tc = (%f,%f)\n", t.a.x, t.a.y, t.b.x, t.b.y, t.c.x, t.c.y); 

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

          //t.ca = vec4d(0.057f, 0.25f, 0.42f, 0.35f);
          //t.cb = vec4d(0.057f, 0.25f, 0.42f, 0.35f);
          //t.cc = vec4d(0.057f, 0.25f, 0.42f, 0.35f);
          
          t.ca = vec4d(0.457f, 0.45f, 0.42f, 0.85f);
          t.cb = vec4d(0.457f, 0.45f, 0.42f, 0.85f);
          t.cc = vec4d(0.457f, 0.45f, 0.42f, 0.85f);

          mesh.faces.push_back(t);			

          b = c;
          nb = nc;
        }
      }
      break;
    }
    default: {
      //std::clog << "ignored line: " << line << std::endl;
    }
    }
  }
  
   // e Bounding Box 
    //mesh.aabb.merge(mesh.aabb.min - mesh.aabb.extents()/10.0);
    //mesh.aabb.merge(mesh.aabb.max + mesh.aabb.extents()/10.0);    

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

bool segmentSphereIntersection3D(const Segment &seg, const Sphere sph, vec3d & out0, vec3d & out1) {
  vec3d o = seg.a;
  vec3d d = seg.b - seg.a;

  vec3d c = sph.center;
  double r = sph.radius;

  double A = dot(d,d);

  double B = 2 * dot((o-c),d);

  double C = dot(o-c,o-c) - r*r;

  double delta = B*B - 4*A*C;

  if(delta <= 0)
    return false;

  double q = (-B - sqrt(delta))/2.;

  if(B < 0)
    q += sqrt(delta);

  double t0 = q/A;
  double t1 = C/q;

  out0 = o + t0 * d;
  out1 = o + t1 * d;

  return true;

}


vec3d cartesianToSpherical(vec3d point) {
  vec3d ans;
  double a = point.x, b = point.y, c = point.z;
  ans.x = sqrt(a*a + b*b + c*c);
  ans.y = acos(c/ans.x);
  ans.z = atan2(b,a);
  return ans;
}
