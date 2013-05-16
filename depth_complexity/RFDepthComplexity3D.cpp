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

Affine3f getModelview();

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

RFDepthComplexity3D::RFDepthComplexity3D(int res, int discretSteps):
    _resolution(res),
    _discretSteps(discretSteps),
    _maximum(0),
    _computeHistogram(true),
    _computeMaximumRays(true),
    _computeGoodRays(true){
    //_computeRaysFromFile(false){

    _threshold = 10;
    _limitRays = 100000;
    _computeRays = true;
    _compute3Dtexture = false;
    _isRaysLimited = false;
    //set up shaders
    ShaderMgr shaderMgr; 
    _shaderCountDC = shaderMgr.createShaderProgram("shader/dc.vert", "shader/dc.frag");        
    _shaderclearBuffer = shaderMgr.createShaderProgram("shader/clear.vert", "shader/clear.frag");
    
    _fboWidth = _fboHeight = res;

    shaderMgr.linkShaderProgram(_shaderCountDC); 
    shaderMgr.linkShaderProgram(_shaderclearBuffer);

    shaderMgr.checkShaderStatus();  
    _goodRays.resize(1);

    // initialize 3d texture (size x, size y, size z, channels, px values)
    // --> channel 1 -> DC of the ray
    // --> chanell 2 -> ray counter
    //_tex3D = CImg<unsigned int>(WIDTH,HEIGHT,DEPTH,2,0);
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
    glProgramUniform2i(_shaderclearBuffer, glGetUniformLocation(_shaderclearBuffer, "resolution"), _fboWidth, _fboHeight);
    
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
    glProgramUniform2i(_shaderCountDC, glGetUniformLocation(_shaderCountDC, "resolution"), _fboWidth, _fboHeight);

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

void RFDepthComplexity3D::setDiscreteSteps(int _dc) {
    _discretSteps = _dc;
}

void RFDepthComplexity3D::setResolution(int _res){
    _resolution = _res;
}

void RFDepthComplexity3D::writeHistogram(std::ostream& out) {
    std::cout << "-- Computing Histogram --" << std::endl;
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
  glBindImageTexture(0, _counterBuffId, 0, GL_FALSE, 0,  GL_READ_WRITE, GL_R32UI);
  glBindTexture(GL_TEXTURE_2D, 0);
    
}

/*
void RFDepthComplexity3D::erodeTriangle(vec3d &v1, vec3d &v2, vec3d &v3){
    
    double pixelSize = 1.0/_fboWidth;
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
 */

//#define CONSTANT_FACTOR 500
void RFDepthComplexity3D::process(const TriMesh &mesh) {
	this->_mesh = &mesh;
        _intersectionTriangles.clear();                 
	_usedPlanes.clear();
	_goodRays.clear();
	_goodRays.resize(1);
	_histogram.clear();
	_histogram.resize(1,0);
        _vpoints.clear();
	_maximum = 0;
        _fboWidth = _fboHeight = _resolution;
        std::cout << _fboWidth << "x" << _fboHeight << std::endl;
                
        initTextureCounter();

	BoundingBox aabb = _mesh->aabb;
	//aabb.merge(aabb.min - aabb.extents()/10.0);
	//aabb.merge(aabb.max + aabb.extents()/10.0);
        
        // initialize 3D texture
        tex3d.CreateTexture3D(512,512,1,3,0);
        tex3d.setMeshBoundingbox(aabb);
	
	vec3d center = _mesh->aabb.center();
	
	 /* initialize random seed: */
        //srand ( time(NULL) );
        
        float modelsizex = aabb.extents().x;
        float modelsizey = aabb.extents().y;
        float modelsizez = aabb.extents().z;
        double distance = sqrt(modelsizex*modelsizex + modelsizey*modelsizey + modelsizez*modelsizez);
        
        int qty = _discretSteps;        
        double inc = M_PI/60.0;
        vec3d ans;  
        
        for (float theta=0; theta <= M_PI/2.0; theta+=M_PI/(double)qty){ 
            /**/
            int length = round((double)qty*2.0 * sin(theta));
            if (length==0) length = 1;
            inc = 2*M_PI/(double)length;
            for (float phi=-M_PI; phi < M_PI; phi+=inc){
                
                ans.x = distance*sin(theta)*cos(phi); 
                ans.y = distance*sin(theta)*sin(phi);
                ans.z = distance*cos(theta);

                _vpoints.push_back(ans+center);

                unsigned int m = renderScene(ans+center);
                //std::cout << "m: " << m << std::endl;
                if (m > _maximum){
                    _maximum = m;               
                }        
            }
        }
        
        std::cout << "vpoint: " << _vpoints.size() << std::endl;
        
        // create histrogram
        //_histogram.resize(_maximum+1);
        //_histogram[_maximum] = _maximumRays.size();
        
        //for (unsigned i=0; i< _goodRays.size(); ++i){
        //    _histogram[i] = _goodRays[i].size();
        //}
        
        
        //std::cout << "fuck\n";
        if (_compute3Dtexture)
              tex3d.buildGLTexture();
        //std::cout << "fuck2\n";
        
        //tex3d.cimg2Tex(_maximum);
  
        
        // Count Intersections
        
        unsigned int max=0;
        std::vector<Point> inter;
        for (std::set<Segment,classcomp>::iterator it = _maximumRays.begin(); it!=_maximumRays.end(); ++it){
            processMeshSegment(*it,&inter);        
            if (max<inter.size()) max = inter.size();
            inter.clear();
            break;
        }
        
        std::cout << "Count intersections: " << max << std::endl;
}

unsigned int RFDepthComplexity3D::renderScene(vec3d point){

    glPushAttrib(GL_VIEWPORT_BIT | GL_ENABLE_BIT);
    glViewport(0, 0, _fboWidth, _fboHeight);

    // desable all anti-aliasing to avoid precision errors
    glDisable(GL_POINT_SMOOTH);
    glDisable(GL_LINE_SMOOTH);
    glDisable(GL_POLYGON_SMOOTH);

    glHint(GL_POLYGON_SMOOTH_HINT,GL_FASTEST);
    glHint(GL_POINT_SMOOTH, GL_FASTEST);
    glHint(GL_LINE_SMOOTH, GL_FASTEST);
    
    glDisable(GL_DEPTH_TEST);
    glDisable(GL_LIGHTING);
    glDisable(GL_BLEND);
    
    // variable computations
    vec3d up = vec3d(0,1,0);
    vec3d center = _mesh->aabb.center();
    vec3d pos;
    vec3d forward = (center-point); forward.normalize();
    vec3d side = cross(forward,up);
    
    double dist = sqrt( pow((point.x-center.x),2)+ pow((point.y-center.y),2)+ pow((point.z-center.z),2));
    pos.z = 2*dist;
        
    glMatrixMode(GL_MODELVIEW);
    glPushMatrix();
    glLoadIdentity();
    
    gluLookAt(point.x, point.y, point.z, center.x, center.y, center.z, 0, 1, 0);
    
    // get modelview matrix
    Vector4f v;
    vec3d vmin(INFINITY,INFINITY,0),vmax(-INFINITY,-INFINITY,0);
    Affine3f modelview = getModelview();
    
    for (unsigned int i=0; i<_mesh->vertices.size(); i++){
        v << _mesh->vertices[i].x, _mesh->vertices[i].y, _mesh->vertices[i].z, 1;
        v = modelview*v;
        if (vmin.x > v[0]) vmin.x = v[0];
        if (vmin.y > v[1]) vmin.y = v[1];
        if (vmax.x < v[0]) vmax.x = v[0];
        if (vmax.y < v[1]) vmax.y = v[1];      
    }

    double relaxation = 1.2; // 20 percent
    pos.x = fabs(vmax.x - vmin.x)*relaxation/2.0;
    pos.y = fabs(vmax.y - vmin.y)*relaxation/2.0;
    
    glMatrixMode(GL_PROJECTION);
    glPushMatrix();
    glLoadIdentity();  
   
    glOrtho(-pos.x, pos.x, -pos.y, pos.y ,1, pos.z);
    
    glMatrixMode(GL_MODELVIEW);
    
    
    glActiveTexture(GL_TEXTURE0);
    glEnable(GL_TEXTURE_2D);
    glBindTexture(GL_TEXTURE_2D, _counterBuffId);  
    
     // set shader to clear the counter buffer ---------
    setShaderClearCounterBuffer();
    
    std::cout << "tttt: " << findMaxValueInCounterBuffer() << std::endl;
    
    // Ensure that all texture writing is done
    glMemoryBarrierEXT(GL_FRAMEBUFFER_BARRIER_BIT_EXT);
    
    // Set shader to write DC in the counter buffer 
    setShaderCountDC();
   
    glColor4f(1,1,1,1);
    glBegin(GL_TRIANGLES);
    for (unsigned i=0; i< _mesh->faces.size(); i++){
        vec3d   v1(_mesh->faces[i].a.x, _mesh->faces[i].a.y, _mesh->faces[i].a.z),
                v2(_mesh->faces[i].b.x, _mesh->faces[i].b.y, _mesh->faces[i].b.z),
                v3(_mesh->faces[i].c.x, _mesh->faces[i].c.y, _mesh->faces[i].c.z);
        
        glVertex3f(v1.x,v1.y,v1.z);
        glVertex3f(v2.x,v2.y,v2.z);
        glVertex3f(v3.x,v3.y,v3.z);
    }
    glEnd();
        
    // Ensure that all texture writing is donetranslationY
    glMemoryBarrierEXT(GL_TEXTURE_UPDATE_BARRIER_BIT_EXT);
    
    // Disable Shaders
    unsigned int max = findMaxValueInCounterBuffer()-1;
    findMaximumRaysAndHistogram(point, forward, side, pos);
           
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
    glEnable(GL_DEPTH_TEST);
    return max;
}

Affine3f getModelview(){
    float m[16];
    glGetFloatv(GL_MODELVIEW_MATRIX, m);
    
    return Affine3f(Matrix4f(m));
}




void RFDepthComplexity3D::findMaximumRaysAndHistogram(vec3d initPos, vec3d f, vec3d s, vec3d p) {
    glPushAttrib(GL_VIEWPORT_BIT);
    glViewport(0.0, 0.0, _fboWidth, _fboHeight);
    
    unsigned int colorBuffer[_fboWidth*_fboHeight];
    
    glBindTexture(GL_TEXTURE_2D, _counterBuffId);
    glGetTexImage( GL_TEXTURE_2D, 0 , GL_RED_INTEGER, GL_UNSIGNED_INT, colorBuffer );
    glBindTexture(GL_TEXTURE_2D, 0);
    
    // compute up vector
    vec3d up = cross(s,f);
    up.normalize();
    s.normalize();
    
    for(int r=0; r<_fboHeight; ++r) {
      for(int c=0; c<_fboWidth; ++c) {
          
        unsigned int val = colorBuffer[r*_fboWidth+c]-1;
        if (val < _threshold) continue;        
                
        Segment seg;
        double pxSizeX = (2.0*p.x)/(double)_fboWidth;
        double pxSizeY = (2.0*p.y)/(double)_fboHeight;
        double translationX = (double)c*pxSizeX - p.x + pxSizeX/2.0;
        double translationY = (double)r*pxSizeY - p.y + pxSizeY/2.0;
        seg.a = initPos + translationX*s + translationY*up;
        seg.b = seg.a + p.z*f;
        
        std::cout << "antes" << std::endl;
        insertRays(val, seg);
        std::cout << "depois" << std::endl;
      }
    }
    
}

bool RFDepthComplexity3D::intersectTriangleSegment(const Segment& segment, const Triangle& tri, Point *pnt) {
	assert(pnt);

	if(!intersectPlaneSegment(makePlane(tri.a, tri.b, tri.c),segment.a,segment.b,pnt))
		return false;

	vec3d u = tri.b - tri.a;
	vec3d v = tri.c - tri.a;
	vec3d w = *pnt - tri.a;

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
   
void RFDepthComplexity3D::insertRays(unsigned int tempMax, Segment seg){
    bool includeRay = !_computeRays;
    
    if (!_computeRays) {
        _goodRays.resize(tempMax+1);
    }
     std::cout << tempMax << std::endl;
    if (tempMax == _maximum){
        if (_computeRays && ( !_isRaysLimited || _maximumRays.size()<_limitRays)) {
            _maximumRays.insert(seg);
            includeRay = true;
        }
        if (_histogram.size() <= tempMax){
            _histogram.resize(tempMax+1,0);
        }
        _histogram[tempMax] += 1;
    }
    else if (tempMax > _maximum){
         _maximum = tempMax;
         if (_computeRays){
            _goodRays.resize(tempMax+1);
            _goodRays[tempMax].insert(_maximumRays.begin(), _maximumRays.end());
            _maximumRays.clear();
            _maximumRays.insert(seg);
         }
        _histogram.resize(tempMax+1, 0);
        _histogram[tempMax] = 1;

        includeRay = true;
    } else {        
        if (_computeRays && ( !_isRaysLimited || _goodRays[tempMax].size() < _limitRays)){
            _goodRays[tempMax].insert(seg);
            includeRay = true;
        }
        _histogram[tempMax] += 1;
    }
        
    if (_compute3Dtexture && includeRay)
        tex3d.updateTexture3D(seg,tempMax);
}

unsigned int RFDepthComplexity3D::findMaxValueInCounterBuffer() {
  const int pixelNumber = _fboWidth * _fboHeight;
  unsigned int colorBuffer[pixelNumber];
  float buff[pixelNumber];
  
  glBindTexture(GL_TEXTURE_2D, _counterBuffId);
  glGetTexImage( GL_TEXTURE_2D, 0 , GL_RED_INTEGER, GL_UNSIGNED_INT, colorBuffer ); 
  glBindTexture(GL_TEXTURE_2D, 0);
  
//   for (int i=0; i<pixelNumber; i++){
//       std::cout << colorBuffer[i] << std::endl;
//       if (colorBuffer[i]==0)
//           std::cin.get();
//   }
  
  unsigned int max =  *(std::max_element(colorBuffer, colorBuffer + pixelNumber))-1;
  
  for (int i=0; i<pixelNumber; i++)
      buff[i] = (colorBuffer[i]-1)/(float)max;
  
  CImg<float> cb(buff,_fboWidth,_fboHeight,1,1);
  
  dualDisplay.resize(_fboWidth, _fboHeight, true);
  dualDisplay.display(cb);
  
  std::cin.get();
  
  return max+1;
}


void RFDepthComplexity3D::processMeshSegment(const Segment& segment, std::vector<Point> *points) {
	assert(points);

	for (unsigned i=0; i<_mesh->faces.size(); ++i) {
		Point p;
		if (intersectTriangleSegment(segment, _mesh->faces[i], &p)){
			points->push_back(p);
                        _intersectionTriangles.push_back(_mesh->faces[i]);
                        _intersectionTriangles.back().ca = vec4d(1,0,0,1);
                        _intersectionTriangles.back().cb = vec4d(1,0,0,1);
                        _intersectionTriangles.back().cc = vec4d(1,0,0,1);
                }
	}
}


bool RFDepthComplexity3D::intersectPlaneSegment(const vec4d& plane, const vec3d& p0, const vec3d& p1, vec3d *pt) {
	double num = -plane.w - dot(plane.xyz(), p0);
	double den = dot(plane.xyz(), p1 - p0);
	if (fabs(den) < EPS)
		return false;
	double r = num / den;
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


