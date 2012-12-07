/* 
 * File:   RFDepthComplexity3D.cpp
 * Author: wagner
 * 
 * Created on 4 de Dezembro de 2012, 13:51
 */

#include "RFDepthComplexity3D.h"

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

void RFDepthComplexity3D::writeRays(std::ostream& out, const std::set<Segment,classcomp> & _rays, int dc) {
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

void RFDepthComplexity3D::writeMaximumRays(std::ostream& out){
    writeRays(out, _maximumRays, _maximum);
}

void RFDepthComplexity3D::writeGoodRays(std::ostream& out){
    for(unsigned i = _threshold; i <= _maximum; ++i) {
        writeRays(out, _goodRays[i], i);
    }
}

RFDepthComplexity3D::RFDepthComplexity3D(int fboWidth, int fboHeight, int discretSteps):
	_fboWidth(fboWidth),
	_fboHeight(fboHeight),
	_discretSteps(discretSteps),
	_maximum(0),
	_computeHistogram(false),
	_computeMaximumRays(false),
	_computeGoodRays(false),
	_computeRaysFromFile(false){
	
        //set up shaders
        ShaderMgr shaderMgr; 
        _shaderCountDC = shaderMgr.createShaderProgram("shader/dc.vert", "shader/dc.frag");        
        _shaderclearBuffer = shaderMgr.createShaderProgram("shader/clear.vert", "shader/clear.frag");
        
                        
        shaderMgr.linkShaderProgram(_shaderCountDC); 
        shaderMgr.linkShaderProgram(_shaderclearBuffer);

        shaderMgr.checkShaderStatus();
  
	_goodRays.resize(1);
}

void RFDepthComplexity3D::setShaderClearCounterBuffer(){
    glUseProgram(_shaderclearBuffer);
    
    
    
    glMatrixMode(GL_PROJECTION);
    glPushMatrix();
    glLoadIdentity();
    gluOrtho2D(0,1,0,1);
    
    glMatrixMode(GL_MODELVIEW);
    glPushMatrix();
    glLoadIdentity();
    
    // Pass MODELVIEW and PROJECTION matrices to Shader
    GLfloat model[16],projection[16];
    glGetFloatv(GL_MODELVIEW_MATRIX, model);
    glGetFloatv(GL_PROJECTION_MATRIX, projection);
    glProgramUniformMatrix4fv(_shaderclearBuffer, glGetUniformLocation(_shaderclearBuffer, "modelViewMat"), 1, GL_FALSE, model);
    glProgramUniformMatrix4fv(_shaderclearBuffer, glGetUniformLocation(_shaderclearBuffer, "projectionMat"), 1, GL_FALSE, projection);

    // Pass Counter Buffer Texture
    glProgramUniform1iEXT(_shaderclearBuffer, glGetUniformLocation(_shaderclearBuffer, "counterBuff"), 0);
    
    glBegin(GL_QUADS);
      glVertex2f(-1.0f , 2.0f);
      glVertex2f(2.0f, 2.0f);
      glVertex2f(2.0f, -1.0f);
      glVertex2f(-1.0f, -1.0f);
    glEnd();
    
    glMatrixMode(GL_PROJECTION);
    glPopMatrix();
    
    glMatrixMode(GL_MODELVIEW);
    glPopMatrix();
}

void RFDepthComplexity3D::setShaderCountDC(){
    // Set shader to count DC
    glUseProgram(_shaderCountDC);   
	
    // Pass MODELVIEW and PROJECTION matrices to Shader
    GLfloat model[16],projection[16];
    glGetFloatv(GL_MODELVIEW_MATRIX, model);
    glGetFloatv(GL_PROJECTION_MATRIX, projection);
    glProgramUniformMatrix4fv(_shaderCountDC, glGetUniformLocation(_shaderCountDC, "modelViewMat"), 1, GL_FALSE, model);
    glProgramUniformMatrix4fv(_shaderCountDC, glGetUniformLocation(_shaderCountDC, "projectionMat"), 1, GL_FALSE, projection);

    // Pass counter buff texture
    glProgramUniform1iEXT(_shaderCountDC, glGetUniformLocation(_shaderCountDC, "counterBuff"), 0);
}



RFDepthComplexity3D::~RFDepthComplexity3D() {
}


void RFDepthComplexity3D::setComputeHistogram(bool computeHistogram) {
	this->_computeHistogram = computeHistogram;
}

void RFDepthComplexity3D::setComputeMaximumRays(bool computeMaximumRays) {
	this->_computeMaximumRays = computeMaximumRays;
}

void RFDepthComplexity3D::setComputeGoodRays(bool computeGoodRays) {
	this->_computeGoodRays = computeGoodRays;
}

void RFDepthComplexity3D::writeRaysSpherical(std::ostream& out, int k) {
  //assert(_computeGoodRays);
  long long unsigned int total = 0;
  for(int i = 0 ; i <= k ; ++i) {
    int ind = maximum()-i;
    if (ind<0) break;
    const std::set<Segment,classcomp> & _rays = goodRays(maximum()-i);
    total += _rays.size();
  }
  out << total << "\n";

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

void RFDepthComplexity3D::setThreshold(unsigned threshold) {
	this->_threshold = threshold;
}

void RFDepthComplexity3D::writeHistogram(std::ostream& out) {
	out << "Intersections Frequency\n";
	for (unsigned i=0; i< _histogram.size(); ++i)
		out << std::left << std::setw(14) << i << _histogram[i] << "\n";
}

void RFDepthComplexity3D::writeRays(std::ostream& out) {
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

void RFDepthComplexity3D::writeRays(std::ostream& out, const std::set<Segment,classcomp> & _rays) {
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

void RFDepthComplexity3D::initTextureCounter(){
    
  //Use a texture as a COUNTER per-pixel Buffer ------------------------
  glGenTextures(1, &_counterBuffId);
  glBindTexture(GL_TEXTURE_2D, _counterBuffId);
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
  glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
  glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
  
  //bind texture to a image unit  
  glTexImage2D(GL_TEXTURE_2D, 0, GL_R32UI, _fboWidth, _fboHeight, 0,  GL_RED_INTEGER, GL_UNSIGNED_INT, 0);
  glBindImageTextureEXT(0, _counterBuffId, 0, GL_FALSE, 0,  GL_READ_WRITE, GL_R32UI);
  glBindTexture(GL_TEXTURE_2D, 0);
    
}

void RFDepthComplexity3D::erodeTriangle(vec3d &v1, vec3d &v2, vec3d &v3){
    
    double pixelSize = 1.0/512.0;
    double step = 0.5*pixelSize*sqrt(2);
    
    vec3d v1dir = (v2 + v3)*0.5 - v1;
    vec3d v2dir = (v1 + v3)*0.5 - v2;
    vec3d v3dir = (v1 + v2)*0.5 - v3;
    
    v1dir.normalize();
    v2dir.normalize();
    v3dir.normalize();
    
    v1 = v1 + step*v1dir;
    v2 = v2 + step*v2dir;
    v3 = v3 + step*v3dir;    
}

void RFDepthComplexity3D::renderScene(vec3d point){
    //glClear( GL_COLOR_BUFFER_BIT|GL_DEPTH_BUFFER_BIT );
     
    //glColorMaterial(GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE);
    //glEnable(GL_COLOR_MATERIAL);
    
    //std::clog << "sorting...";
    //if (sorted_faces.empty()) {
    //    sorted_faces = _mesh.faces;
    //}
    //glMatrixMode(GL_PROJECTION);
    //glPushMatrix();
    //glLoadIdentity();
    
    //std::sort(sorted_faces.begin(), sorted_faces.end(), ByDist(dir));
    //std::clog << "done" << std::endl;
    glPushAttrib(GL_VIEWPORT_BIT | GL_ENABLE_BIT);
    glViewport(0, 0, _fboWidth, _fboHeight);

    // desable all anti-aliasing to avoid precision error
    glDisable(GL_POINT_SMOOTH);
    glDisable(GL_LINE_SMOOTH);
    glDisable(GL_POLYGON_SMOOTH);

    glHint(GL_POLYGON_SMOOTH_HINT,GL_FASTEST);
    glHint(GL_POINT_SMOOTH, GL_FASTEST);
    glHint(GL_LINE_SMOOTH, GL_FASTEST);
    
    
    glDisable(GL_DEPTH_TEST);
    glDisable(GL_LIGHTING);
    glDisable(GL_BLEND);
    //glDisable(GL_DEPTH_TEST);
    
    glMatrixMode(GL_PROJECTION);
    glPushMatrix();
    glLoadIdentity();
    gluPerspective(50, 1024./768., 0.1, 10000);
    
    glMatrixMode(GL_MODELVIEW);
    glPushMatrix();
    glLoadIdentity();
    
    vec3d center = _mesh->aabb.center();
    gluLookAt(point.x, point.y, point.z, center.x, center.y, center.z, 0, 1, 0  );

    glActiveTexture(GL_TEXTURE0);
    glEnable(GL_TEXTURE_2D);
    glBindTexture(GL_TEXTURE_2D, _counterBuffId);  
    
     // set shader to clear the counter buffer ---------
    setShaderClearCounterBuffer();
    
    // Ensure that all texture writing is done
    glMemoryBarrierEXT(GL_FRAMEBUFFER_BARRIER_BIT_EXT);
    
    //int max = findMaxValueInCounterBuffer()-1;
    //std::cout << "initialize test: " << max << std::endl;
    
    // Set shader to write DC in the counter buffer 
    setShaderCountDC();
   
    //glLightModeli(GL_LIGHT_MODEL_TWO_SIDE, 1);    
    glColor4f(1,1,1,1);
    glBegin(GL_TRIANGLES);
    for (unsigned i=0; i< _mesh->faces.size(); i++){
        vec3d   v1(_mesh->faces[i].a.x, _mesh->faces[i].a.y, _mesh->faces[i].a.z),
                v2(_mesh->faces[i].b.x, _mesh->faces[i].b.y, _mesh->faces[i].b.z),
                v3(_mesh->faces[i].c.x, _mesh->faces[i].c.y, _mesh->faces[i].c.z);
        
        //erodeTriangle(v1,v2,v3);
        
        glVertex3f(v1.x,v1.y,v1.z);
        glVertex3f(v2.x,v2.y,v2.z);
        glVertex3f(v3.x,v3.y,v3.z);
        
        //glVertex3f(_mesh->faces[i  ].a.x, _mesh->faces[i  ].a.y, _mesh->faces[i  ].a.z );
        //glVertex3f(_mesh->faces[i+1].b.x, _mesh->faces[i+1].b.y, _mesh->faces[i+1].b.z);
        //glVertex3f(_mesh->faces[i+2].c.x, _mesh->faces[i+2].c.y, _mesh->faces[i+2].c.z);
    }
    glEnd();
    //glEnable(GL_VERTEX_ARRAY);

    //glEnableClientState(GL_VERTEX_ARRAY);    
    //glVertexPointer(3, GL_DOUBLE, 2*sizeof(vec3d)+sizeof(vec4d), &_mesh->faces[0].a.x);
    
    //glEnableClientState(GL_COLOR_ARRAY);
    //glColorPointer(4, GL_DOUBLE, 2*sizeof(vec3d)+sizeof(vec4d), &sorted_faces[0].ca.x);

    //glEnableClientState(GL_NORMAL_ARRAY);
    //glNormalPointer(GL_DOUBLE, 2*sizeof(vec3d)+sizeof(vec4d), &sorted_faces[0].na.x);

    //glDrawArrays(GL_TRIANGLES, 0, _mesh->faces.size()*3);

    //glDisableClientState(GL_VERTEX_ARRAY);
    //glDisable(GL_VERTEX_ARRAY);
    //glDisableClientState(GL_NORMAL_ARRAY);
    //glDisable(GL_NORMAL_ARRAY);
    //glDisableClientState(GL_COLOR_ARRAY);
    //glDisable(GL_COLOR_ARRAY);
    
    // Ensure that all texture writing is done
    glMemoryBarrierEXT(GL_TEXTURE_UPDATE_BARRIER_BIT_EXT);
    // Disable Shaders
    _maximum = findMaxValueInCounterBuffer()-1;
    //std::cout << "the maximum DC is: " << max << std::endl;
        
    glUseProgram(0);
    glBindTexture(GL_TEXTURE_2D, 0);   

    glMatrixMode(GL_MODELVIEW);
    glPopMatrix();
    
    glMatrixMode(GL_PROJECTION);
    glPopMatrix();
    
    glPopAttrib();
    glEnable(GL_LIGHTING);
    glEnable(GL_BLEND);
    glBlendFunc(GL_ONE, GL_ONE_MINUS_SRC_ALPHA);
    //glEnable(GL_DEPTH_TEST);
    //glDisable(GL_COLOR_MATERIAL);
    glEnable(GL_DEPTH_TEST);

}

#define CONSTANT_FACTOR 500
void RFDepthComplexity3D::process(const TriMesh &mesh) {
	 this->_mesh = &mesh;

	_usedPlanes.clear();
	_goodRays.clear();
	_goodRays.resize(1);
	_histogram.clear();
	_histogram.resize(1);
	_maximum = 0;
        
        initTextureCounter();
	
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
	
	//double r = aabb.extents().length()/2.0;
	vec3d center = aabb.center();
        vec3d radius = aabb.max - center;
        double r = radius.length();
        std::cout << r << std::endl;
	
	 /* initialize random seed: */
        srand ( time(NULL) );

	const unsigned steps = _discretSteps;
	//const unsigned Nrays = (_computeRaysFromFile? _raysFromFile.size() : steps*steps*CONSTANT_FACTOR);
        const unsigned Nrays = steps*steps*CONSTANT_FACTOR;

        std::cout << Nrays << std::endl;
        
        float aspect = 1024./768.;
        float fov = 50;
        float fovh = aspect*fov;
        
        float modelsizex = aabb.extents().x;
        float modelsizey = aabb.extents().y;
        double distance=0;
        
        if (modelsizex/aspect > modelsizey){
            distance = aabb.extents().z/2 + (modelsizex/2)/tan((fovh/2)*0.01745329238);
        } else {
            distance = aabb.extents().z/2 + (modelsizex/2)/tan((fov/2)*0.01745329238);
        }
        
        //generate random vector
        vec3d v(uniformRandom() - 0.5, uniformRandom() - 0.5, uniformRandom() - 0.5);
        v.normalize();
        v = distance*v + center;
        
        std::cout << v << std::endl;
        
        renderScene(v);
        
        /*
	for (unsigned tn = 0; tn < Nrays ; ++tn) {
		Segment randomSegment;


                
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
         */
}

unsigned int RFDepthComplexity3D::findMaxValueInCounterBuffer() {
  const int pixelNumber = _fboWidth * _fboHeight;
  unsigned int colorBuffer[pixelNumber];
  float buff[pixelNumber];
  
  glBindTexture(GL_TEXTURE_2D, _counterBuffId);
  glGetTexImage( GL_TEXTURE_2D, 0 , GL_RED_INTEGER, GL_UNSIGNED_INT, colorBuffer ); 
  glBindTexture(GL_TEXTURE_2D, 0);
  
  // show pic
  
  
  unsigned int max =  *(std::max_element(colorBuffer, colorBuffer + pixelNumber))-1;
  
  for (int i=0; i<pixelNumber; i++)
      buff[i] = (colorBuffer[i]-1)/(float)max;
      //((colorBuffer[i]-1)/((float)max-1.0f) > 0)? 1.0f: 0.0f;
      //std::cout << colorBuffer[i]-1 << std::endl;
  
  CImg<float> cb(buff,512,512,1,1);
  //cb.normalize(0,1);
  
  dualDisplay.display(cb);
  
  
  return max+1;
}

void RFDepthComplexity3D::processMeshSegment(const Segment& segment, std::vector<Point> *points) {
	assert(points);

	for (unsigned i=0; i<_mesh->faces.size(); ++i) {
		Point p;
		if (intersectTriangleSegment(segment, _mesh->faces[i], &p))
			points->push_back(p);
	}
}

bool RFDepthComplexity3D::intersectTriangleSegment(const Segment& segment, const Triangle& tri, Point *pnt) {
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

bool RFDepthComplexity3D::intersectPlaneSegment(const vec4d& plane, const vec3d& p0, const vec3d& p1, vec3d *pt) {
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

vec4d RFDepthComplexity3D::makePlane(const vec3d& a, const vec3d& b, const vec3d& c) {
		vec3d normal = cross(b-a, c-a);
		normal.normalize();
		double d = dot(a, normal);
		return vec4d(normal, -d);
}

