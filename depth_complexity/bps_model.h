/* 
 * File:   bps_model.h
 * Author: sombra
 *
 * Created on 12 de Julho de 2012, 17:28
 */

#ifndef BPS_MODEL_H
#define	BPS_MODEL_H

#include "flags.h"
#include "util.h"
#ifdef USE_RANDOM_DC3D
#include "dc_3d_random.h"
#else
#include "dc_3d.h"
#endif
#include "timer.h"
#include "camera/float3.h"
#include "camera/Camera.h"

#include <iostream>
#include <fstream>
#include <sstream>
#include <cstring>
/*#include <GL/glew.h>
#ifdef __APPLE__
#include <GLUT/glut.h>
#else
#include <GL/glut.h>
#endif*/
//defines usados por getPlaneTriangleIntersection para comunicar quais pontos do
//triangulo estão abaixo ou acima do plano
#define ALL_BELLOW_PLANE      -4
#define ONLY_C_ABOVE_PLANE    -3
#define ONLY_B_ABOVE_PLANE    -2
#define ONLY_A_ABOVE_PLANE    -1
#define ONLY_A_BELLOW_PLANE   1
#define ONLY_B_BELLOW_PLANE   2
#define ONLY_C_BELLOW_PLANE   3
#define ALL_ABOVE_PLANE       4  

void  subDivideModel(const vec4d plane, const TriMesh& model, TriMesh& partA,
        TriMesh& partB);

vec4d makePlane(const vec3d& normal, const vec3d& point);

std::string getExtension(const std::string& filename);

void writeOFFModel(const TriMesh& model, std::ostream& out);

std::vector<Segment> loadOFFLines(std::istream& in);

int getPlaneTriangleIntersection(const vec4d& plane, const Triangle& triangle);

Triangle cutTriangle(const vec4d& plane, const Triangle& triangle);

vec4d defineCuttingPlane(const vector<Segment> lines, const BoundingBox aabb);

Segment findSegmentBetweenLines(Segment lineA, Segment lineB);

//test if a point is inside a aabb(bounding box)
bool insideAABB(const vec3d point, const BoundingBox aabb);

//plane defined by a initial point and a normal vector
typedef struct w_plane{
  vec3d p0;
  vec3d vet_dir;
  double weight;//number of lines used to form the plan
}W_Plane;


#endif	/* BPS_MODEL_H */

