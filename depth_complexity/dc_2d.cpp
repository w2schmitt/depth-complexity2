#include "dc_2d.h"

#include <iostream>
#include <cassert>

// Internal function to check if the framebuffer was created correctly.
bool checkFramebufferStatus();

// Draw polygons in framebuffer.
void drawQuad(const Point &p1, const Point &p2, const Point &p3, const Point &p4);
void drawTriangle(const Point &p1, const Point &p2, const Point &p3);


DepthComplexity2D::DepthComplexity2D(const int fboWidth, const int fboHeight){
  _fboWidth = fboWidth;
  _fboHeight = fboHeight;
  _computeHistogram = false;
  _computeMaximumRays = false;
  _computeGoodRays = false;
  _maximum = 0;
  _threshold = 0;
  
  fbo = new float[fboWidth*fboHeight];
  tmpfbo = new float[fboWidth*fboHeight];


  assert(initFBO());
}

bool DepthComplexity2D::initFBO() {
  // Use texture as COLOR buffer
  //glGenTextures(1, &_textureId);
  //glBindTexture(GL_TEXTURE_2D, _textureId);
  //glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
  //glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
  //glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
  //glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
  //glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA8, _fboWidth, _fboHeight, 0, GL_RGBA, GL_UNSIGNED_BYTE, 0);
  //glBindTexture(GL_TEXTURE_2D, 0);

  // User renderbuffer as DEPTH and STENCIL buffer.
  //glGenRenderbuffersEXT(1, &_rboId);
  //glBindRenderbufferEXT(GL_RENDERBUFFER_EXT, _rboId);
  //glRenderbufferStorageEXT(GL_RENDERBUFFER_EXT, GL_DEPTH_STENCIL_EXT, _fboWidth, _fboHeight);
  //glBindRenderbufferEXT(GL_RENDERBUFFER_EXT, 0);

  // Create framebuffer object
  //glGenFramebuffersEXT(1, &_fboId);
  //glBindFramebufferEXT(GL_FRAMEBUFFER_EXT, _fboId);

  // Assign COLOR, DEPTH and STENCIL buffers to the render buffer.
  //glFramebufferTexture2DEXT(GL_FRAMEBUFFER_EXT, GL_COLOR_ATTACHMENT0_EXT, GL_TEXTURE_2D, _textureId, 0);
  //glFramebufferRenderbufferEXT(GL_FRAMEBUFFER_EXT, GL_DEPTH_ATTACHMENT_EXT, GL_RENDERBUFFER_EXT, _rboId);
  //glFramebufferRenderbufferEXT(GL_FRAMEBUFFER_EXT, GL_STENCIL_ATTACHMENT_EXT, GL_RENDERBUFFER_EXT, _rboId);

  // Check if everything was ok.
  //_status = checkFramebufferStatus();

  //glBindFramebufferEXT(GL_FRAMEBUFFER_EXT, 0);
  
  return true;
}

bool checkFramebufferStatus() {
  // check FBO status
  GLenum status = glCheckFramebufferStatusEXT(GL_FRAMEBUFFER_EXT);
  switch(status) {
  case GL_FRAMEBUFFER_COMPLETE_EXT:
    std::cerr << "Framebuffer complete." << std::endl;
    return true;

  case GL_FRAMEBUFFER_INCOMPLETE_ATTACHMENT_EXT:
    std::cerr << "[ERROR] Framebuffer incomplete: Attachment is NOT complete." << std::endl;
    return false;

  case GL_FRAMEBUFFER_INCOMPLETE_MISSING_ATTACHMENT_EXT:
    std::cerr << "[ERROR] Framebuffer incomplete: No image is attached to FBO." << std::endl;
    return false;

  case GL_FRAMEBUFFER_INCOMPLETE_DIMENSIONS_EXT:
    std::cerr << "[ERROR] Framebuffer incomplete: Attached images have different dimensions." << std::endl;
    return false;

  case GL_FRAMEBUFFER_INCOMPLETE_FORMATS_EXT:
    std::cerr << "[ERROR] Framebuffer incomplete: Color attached images have different internal formats." << std::endl;
    return false;

  case GL_FRAMEBUFFER_INCOMPLETE_DRAW_BUFFER_EXT:
    std::cerr << "[ERROR] Framebuffer incomplete: Draw buffer." << std::endl;
    return false;

  case GL_FRAMEBUFFER_INCOMPLETE_READ_BUFFER_EXT:
    std::cerr << "[ERROR] Framebuffer incomplete: Read buffer." << std::endl;
    return false;

  case GL_FRAMEBUFFER_UNSUPPORTED_EXT:
    std::cerr << "[ERROR] Unsupported by FBO implementation." << std::endl;
    return false;

  default:
    std::cerr << "[ERROR] Unknown error." << std::endl;
    return false;
  }
}

void DepthComplexity2D::process(
  const Segment &from, const Segment &to, const std::vector<Segment> &segments) {
  _from = from;
  _to = to;
  _segments = &segments;
  
    glClearColor(0, 0, 0, 1.0);
    glEnable(GL_COLOR_MATERIAL);
	glColorMaterial(GL_FRONT, GL_AMBIENT_AND_DIFFUSE);
	
	glDisable(GL_LINE_SMOOTH);
	
	glEnable(GL_BLEND);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
  
  //std::cout << segments.size();
  
  /*std::cout << _textureId << " " << _fboId << " " << _rboId << " " << _status << " "
            << _computeHistogram << " " << _computeMaximumRays << " " << _computeGoodRays << " "
            << _fboWidth << " " << _fboHeight << " " << _segments->size() << " " << _threshold << std::endl;// " " << _from << " " << _to << " " */

  findDepthComplexity2D();
  //glPixelTransferf(GL_MAP_STENCIL,GL_TRUE);
  if (_computeHistogram or _computeMaximumRays or _computeGoodRays)
    findMaximumRaysAndHistogram();
    
  
  //copyStencilToColor();
}

void DepthComplexity2D::clearBuffer(float* buff, int w, int h, int ch){
	for (int i=0; i< w*h*ch; ++i){ buff[i] = 0.;}
}

// increment buff1 with values from buff2
void DepthComplexity2D::sumBuffer(float* buff1, float* buff2, int w, int h, int ch){
	for (int i=0; i< w*h*ch; ++i){ buff1[i] += buff2[i];}
}

void DepthComplexity2D::copyBuffer(float* buff1, float* buff2, int w, int h, int ch){
	for (int i=0; i< w*h*ch; ++i){ buff1[i] = buff2[i];}
}


void DepthComplexity2D::findDepthComplexity2D() {
  //glBindFramebufferEXT(GL_FRAMEBUFFER_EXT, _fboId);
  clearBuffer(fbo, _fboWidth, _fboHeight, 1);
  //clearBuffer(tmpfbo, _fboWidth, _fboHeight, 1);
  
  glPushAttrib(GL_VIEWPORT_BIT);
    glViewport(0, 0, _fboWidth, _fboHeight);

	
		
    glMatrixMode(GL_PROJECTION);
    glPushMatrix();
      glLoadIdentity();
      gluOrtho2D(0.0, 1.0, 0.0, 1.0);
      
      		glMatrixMode(GL_MODELVIEW);
		glPushMatrix();
		glLoadIdentity();
      
      // render polygons in stencil buffer
      //glEnable(GL_STENCIL_TEST);
     
      //glStencilFunc(GL_ALWAYS, 0, ~0);
      //glStencilOp(GL_INCR, GL_INCR, GL_INCR);
     // glColorMask(GL_FALSE, GL_FALSE, GL_FALSE, GL_FALSE);
      //glDepthMask(GL_FALSE);
      glDisable(GL_LINE_SMOOTH);
      glDisable(GL_DEPTH_TEST);
	  //glClearStencil(0);
	  //glStencilMask(~0);
      //glClear(GL_DEPTH_BUFFER_BIT  | GL_COLOR_BUFFER_BIT);
      
      //glMatrixMode(GL_MODELVIEW);
	  //std::cout << "Seg Size: "<< _segments->size() << "\n";
      double pixelSize = 1.0/(double)512.0; //512.0;
      double step = 0.5*pixelSize*1.41421356;
      
      for (unsigned i=0; i<_segments->size(); ++i) {
        if (!_segments->at(i).active) continue;
        Segment Tv0 = computeDualSegmentFromPoint(_segments->at(i).a);
        Segment Tv1 = computeDualSegmentFromPoint(_segments->at(i).b);
		
		//printf("Tv0: a=(%f,%f) , b=(%f,%f)\n", Tv0.a.x, Tv0.a.y, Tv0.b.x, Tv0.b.y);
		//printf("Tv1: a=(%f,%f) , b=(%f,%f)\n", Tv1.a.x, Tv1.a.y, Tv1.b.x, Tv1.b.y);
		//std::cin.get();
		
        if (Tv0.a.y < Tv1.a.y) std::swap(Tv0, Tv1);
        vec3d u = Tv0.b - Tv0.a; u.normalize();
        vec3d v = Tv1.b - Tv1.a; v.normalize();
        
        vec3d dir;

        double t1, t2;
        if (segmentIntersection2D(Tv0, Tv1, &t1, &t2) ) {
          // Shrink triangles
          dir = vec3d(u.y, -u.x);
          Segment s1(Tv0.a + dir*step, Tv0.b + dir*step);
          dir = vec3d(-v.y, v.x);
          Segment s2(Tv1.a + dir*step, Tv1.b + dir*step);
          if (s1.a.y > s2.a.y) {
            double t1, t2;
            if (lineIntersection2D(s1, s2, &t1, &t2) ) {
		      vec3d inter = s1.a*(1.f-t1) + s1.b*t1;
			  //printf("1\n");
			  //printf("%f %f %f %f %f %f\n",s1.a.x, s1.a.y, inter.x, inter.y, s2.a.x, s2.a.y);  
			  glColor3f(1,0,0);        
              drawTriangle(s1.a, inter, s2.a);
              
			  glReadPixels(0, 0, _fboWidth, _fboHeight, GL_RED, GL_FLOAT, tmpfbo);
			  sumBuffer(fbo, tmpfbo, _fboWidth, _fboHeight, 1);
			  glClear(GL_DEPTH_BUFFER_BIT  | GL_COLOR_BUFFER_BIT);
            }
          }
          

          
          // Shrink triangles
          dir = vec3d(-u.y, u.x);
          Segment s3(Tv0.a + dir*step, Tv0.b + dir*step);
          dir = vec3d(v.y, -v.x);
          Segment s4(Tv1.a + dir*step, Tv1.b + dir*step);
          if (s3.b.y < s4.b.y) {
            double t1, t2;
            if (lineIntersection2D(s3, s4, &t1, &t2) ) {
			  vec3d inter = s3.a*(1.f-t1) + s3.b*t1;
              //printf("2\n");
              //printf("%f %f %f %f %f %f\n", s3.b.x, s3.b.y, inter.x, inter.y, s4.b.x, s4.b.y);
              glColor3f(1,1,1);
              drawTriangle(s3.b, inter, s4.b);
              
              glReadPixels(0, 0, _fboWidth, _fboHeight, GL_RED, GL_FLOAT, tmpfbo);
			  sumBuffer(fbo, tmpfbo, _fboWidth, _fboHeight, 1);
			  glClear(GL_DEPTH_BUFFER_BIT  | GL_COLOR_BUFFER_BIT);
            }
          }
        } else {
          dir = vec3d(u.y, -u.x);
          Segment s1(Tv0.a + dir*step, Tv0.b + dir*step);
          dir = vec3d(-v.y, v.x);
          Segment s2(Tv1.a + dir*step, Tv1.b + dir*step);
          if (s1.a.y > s2.a.y and s1.b.y > s2.b.y) {
			 //printf("3\n");
			 //printf("%f %f %f %f %f %f %f %f\n", s1.a.x, s1.a.y, s1.b.x, s1.b.y, s2.b.x, s2.b.y, s2.a.x, s2.a.y);
			 glColor3f(1,1,1);
			 drawQuad(s1.a, s1.b, s2.b, s2.a);
			 glReadPixels(0, 0, _fboWidth, _fboHeight, GL_RED, GL_FLOAT, tmpfbo);
			 sumBuffer(fbo, tmpfbo, _fboWidth, _fboHeight, 1);
			 glClear(GL_DEPTH_BUFFER_BIT  | GL_COLOR_BUFFER_BIT);
          }
        }
      }
      //printf("4\n");

	  //printf("max:%f\n", findMaxValueInStencil());
      _maximum = findMaxValueInStencil();
	 
	  
	  //const int pixelNumber = _fboWidth * _fboHeight;
	  //float stencilBuffer[pixelNumber];
	  //glReadPixels(0, 0, _fboWidth, _fboHeight, GL_STENCIL_INDEX, GL_FLOAT, stencilBuffer);
  
	  
	  //for (int i=0; i< pixelNumber; ++i){
	//	  std::cout << (int)stencilBuffer[i] << "\n";
	 //}
	  
      //glColorMask(GL_TRUE, GL_TRUE, GL_TRUE, GL_TRUE);
      //glDepthMask(GL_TRUE);
      //glDisable(GL_STENCIL_TEST);
      //glBindFramebufferEXT(GL_FRAMEBUFFER_EXT, 0);
    glMatrixMode(GL_PROJECTION);
    glPopMatrix();
  glMatrixMode(GL_MODELVIEW);
  glPopMatrix();
  glPopAttrib();
}

Segment DepthComplexity2D::computeDualSegmentFromPoint(const Point &p) {
  Segment seg;
  double t1, t2;

  // Finding left point
  bool r1 = lineIntersection3D(_to, Segment(_from.a, p), &t1, &t2);
  assert(r1);
  seg.a = Point(0.0, t1);

  // Finding right point
  bool r2 = lineIntersection3D(_to, Segment(_from.b, p), &t1, &t2);
  assert(r2);
  seg.b = Point(1.0, t1);

  return seg;
}

void drawQuad(const Point &p1, const Point &p2, const Point &p3, const Point &p4) {
  glBegin(GL_QUADS);
    glVertex2f(p1.x, p1.y);
    glVertex2f(p2.x, p2.y);
    glVertex2f(p3.x, p3.y);
    glVertex2f(p4.x, p4.y);
  glEnd();
}

void drawTriangle(const Point &p1, const Point &p2, const Point &p3) {
  glBegin(GL_TRIANGLES);
    glVertex2f(p1.x, p1.y);
    glVertex2f(p2.x, p2.y);
    glVertex2f(p3.x, p3.y);
  glEnd();
}

float DepthComplexity2D::findMaxValueInStencil() {
  const int pixelNumber = _fboWidth * _fboHeight;
  //float stencilBuffer[pixelNumber];
  //glReadPixels(0, 0, _fboWidth, _fboHeight, GL_STENCIL_INDEX, GL_FLOAT, stencilBuffer);
  return *(std::max_element(fbo, fbo + pixelNumber));
  //std::cout << "-MAX:"<< max << "\n";
  
  
}

void DepthComplexity2D::findMaximumRaysAndHistogram() {
  _maximumRays.clear();
  _goodRays.clear();
  _goodRays.resize(_maximum + 1);
  _histogram.clear();
  _histogram.resize(_maximum + 1);
  
  //glPixelTransferi(GL_MAP_STENCIL, GL_TRUE);

  //glBindFramebufferEXT(GL_FRAMEBUFFER_EXT, _fboId);
  glPushAttrib(GL_VIEWPORT_BIT);
    glViewport(0.0, 0.0, _fboWidth, _fboHeight);
    
	
    //GLubyte stencilSave[_fboWidth*_fboHeight];
    //glReadPixels(0, 0, _fboWidth, _fboHeight, GL_STENCIL_INDEX, GL_UNSIGNED_BYTE, stencilSave);
    //printf("stencil: \n");
    //unsigned value = *(std::max_element(stencilSave, stencilSave + _fboWidth*_fboHeight));
    //std::cout << "maxxx: " << value << "\n";
	//std::cout << "height: " << _fboHeight << "\n";
    for(int r=0; r<_fboHeight; ++r) {
      for(int c=0; c<_fboWidth; ++c) {
        //GLubyte val = stencilSave[r*_fboWidth+c];
       
        //printf("val:%f\n", fbo[r*_fboWidth+c]);
        unsigned int val = fbo[r*_fboWidth+c];
        
        //printf("val: %d",val);
		//if (val==7)
		//	std::cout << "eh isso aí\n";
        if (_computeHistogram) _histogram[val]++;

        if ((_computeMaximumRays && val == _maximum) || (_computeGoodRays && val >= _threshold)) {
          Segment seg;
          double t1 = c/(double)_fboWidth;
          double t2 = r/(double)_fboHeight;
          seg.a = _from.a*(1.f-t1) + _from.b*t1;
          seg.b = _to.a*(1.f-t2) + _to.b*t2;
		  seg.sortPoints();          
          if (val == _maximum){
			 if (_maximumRays.size() < 10)
				_maximumRays.insert(seg);
            }
          else if (val >= _threshold){
           //if (_goodRays[val].size() < 200)
            _goodRays[val].insert(seg);
		  }
        }
      }
    }
  glPopAttrib();
  //glBindFramebufferEXT(GL_FRAMEBUFFER_EXT, 0);
}

float colors[10][3] = {
  {1.0f, 1.0f, 1.0f},
  {1.0f, 0.0f, 0.0f},
  {0.0f, 1.0f, 0.0f},
  {0.0f, 0.0f, 1.0f},
  {0.5f, 0.0f, 0.5f},
  {0.5f, 0.5f, 0.0f},
  {0.0f, 0.5f, 0.5f},
  {0.5f, 0.5f, 0.5f},
  {1.0f, 0.0f, 1.0f},
  {0.0f, 1.0f, 1.0f},
};

void DepthComplexity2D::copyStencilToColor() {
  glBindFramebufferEXT(GL_FRAMEBUFFER_EXT, _fboId);

  glPushAttrib(GL_VIEWPORT_BIT);
    glViewport(0.0, 0.0, _fboWidth, _fboHeight);

    GLfloat colorSave[_fboWidth*_fboHeight*3];
    GLubyte stencilSave[_fboWidth*_fboHeight];
    GLint	previousColorBuffer;
    glReadPixels(0, 0, _fboWidth, _fboHeight, GL_STENCIL_INDEX, GL_UNSIGNED_BYTE, stencilSave);

    // I'm sure this could be done much better with OpenGL
    for(register int y = 0; y < _fboHeight; ++y) {
      for(register int x = 0; x < _fboWidth; ++x) {
        int stencilValue  = stencilSave[_fboWidth * y + x];
        colorSave[(_fboWidth * y + x) * 3 + 0] = colors[stencilValue % 10][0];
        colorSave[(_fboWidth * y + x) * 3 + 1] = colors[stencilValue % 10][1];
        colorSave[(_fboWidth * y + x) * 3 + 2] = colors[stencilValue % 10][2];
      }
    }

    glMatrixMode(GL_MODELVIEW);
    glPushMatrix();
      glLoadIdentity();

      glMatrixMode(GL_PROJECTION);
      glPushMatrix();
        glLoadIdentity();
        glOrtho(0, 1, 0, 1, 0, 1);

        glDisable(GL_DEPTH_TEST);
        glDisable(GL_STENCIL_TEST);
        glColorMask(GL_TRUE, GL_TRUE, GL_TRUE, GL_TRUE);
        glGetIntegerv(GL_DRAW_BUFFER, &previousColorBuffer);
        glDrawBuffer(GL_BACK);
        glDrawPixels(_fboWidth, _fboHeight, GL_RGB, GL_FLOAT, colorSave);
        glDrawBuffer(previousColorBuffer);
        glEnable(GL_DEPTH_TEST);
      glPopMatrix();
    glMatrixMode(GL_MODELVIEW);
    glPopMatrix();

  glPopAttrib();
  glBindFramebufferEXT(GL_FRAMEBUFFER_EXT, 0);
}
