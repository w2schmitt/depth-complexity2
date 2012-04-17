#include "dc_2d.h"

#include <iostream>
#include <cassert>

//const char *DepthComplexity2D::V_PATH = "shader/dc.vert";
//const char *DepthComplexity2D::P_PATH = "shader/dc.frag";

// Internal function to check if the framebuffer was created correctly.
bool checkFramebufferStatus();

// Internal function to check if shader was compiled correctly.
GLuint createShaderProgram(const char* vFilename, const char* pFilename);
GLuint createShader(const char* filename, GLuint shaderType);
void linkShaderProgram(GLuint shaderProgramId);
char *readShaderFile(const char *filename);
bool checkShadersStatus(GLuint shaderID);

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
  //_vShaderId = _vShaderClearId = 0;
  //_pShaderId = _pShaderClearId = 0;
  _shaderProgram = 0;
  _shaderProgram2 = 0;
  _status = true;   
  //_bufferCounter = new unsigned int[_fboHeight*_fboWidth];
  
  assert(initFBO()); 

  //clear shader 
  _shaderProgram = createShaderProgram("shader/dc.vert", "shader/dc.frag");
  _shaderProgram2 = createShaderProgram("shader/clear.vert", "shader/clear.frag");
	
  linkShaderProgram(_shaderProgram);  
  linkShaderProgram(_shaderProgram2);
   
}


DepthComplexity2D::~DepthComplexity2D(){
	//delete [] _bufferCounter;
}

GLuint createShaderProgram(const char* vFilename, const char* pFilename){
	
  GLuint shaderProgram = glCreateProgram();
  
  GLuint vShaderId = createShader(vFilename, GL_VERTEX_SHADER);
  GLuint pShaderId = createShader(pFilename, GL_FRAGMENT_SHADER);
  
  glAttachShader(shaderProgram, vShaderId);
  glAttachShader(shaderProgram, pShaderId);
  
  return shaderProgram;
}

GLuint createShader(const char* filename, GLuint shaderType){
  
  GLuint shaderID = glCreateShader(shaderType);
  
  // Read shader files
	char *shaderData = readShaderFile(filename);
	
	// Pointer to shader data
	const char *ptVData = shaderData;
  
  glShaderSource(shaderID,  1, &ptVData, NULL);
  glCompileShader(shaderID);
  
  //_status = 
  checkShadersStatus(shaderID);

	free(shaderData);
  
  return shaderID;
}

void linkShaderProgram(GLuint shaderProgramId){
  glLinkProgram(shaderProgramId);
  
  GLint ok;
  glGetShaderiv(shaderProgramId, GL_LINK_STATUS, &ok);
  if (!ok){
    int ilength; char stringBuffer[STRING_BUFFER];
    glGetShaderInfoLog(shaderProgramId, STRING_BUFFER, &ilength, stringBuffer);
    std::cerr << "[ERROR] SHADER LINK ERROR : "<<stringBuffer << std::endl; 
  } else { std::cerr << "[OK] SHADER LINK OK." << std::endl; }
  
  glValidateProgram(shaderProgramId);
  
  glGetShaderiv(shaderProgramId, GL_VALIDATE_STATUS, &ok);
  if (!ok){
    int ilength; char stringBuffer[STRING_BUFFER];
    glGetShaderInfoLog(shaderProgramId, STRING_BUFFER, &ilength, stringBuffer);
    std::cerr << "[ERROR] SHADER VALIDATION ERROR : "<<stringBuffer << std::endl; 
  } else { std::cerr << "[OK] SHADER VALID OK." << std::endl; }
}


char* readShaderFile(const char *filename){
	
	FILE *fp;
	char *content = NULL;

	int count=0;

	if (filename != NULL) {
		fp = fopen(filename,"rt");
		if (fp != NULL) {      
			  fseek(fp, 0, SEEK_END);
			  count = ftell(fp);
			  rewind(fp);

				if (count > 0) {
					content = (char *)malloc(sizeof(char) * (count+1));
					count = fread(content,sizeof(char),count,fp);
					content[count] = '\0';
			}
			fclose(fp);
		} else {
			std::cerr << "[ERROR] Cannot load Vertex or Fragment Shader File! " << std::endl;
		}
	}
	
	return content;
}

bool checkShadersStatus(GLuint shaderID){
	
	GLint ok;
	glGetShaderiv(shaderID, GL_COMPILE_STATUS, &ok);
	if (!ok){
		int ilength; char stringBuffer[STRING_BUFFER];
		glGetShaderInfoLog(shaderID, STRING_BUFFER, &ilength, stringBuffer);
		std::cerr<<"[ERROR] Shader Compilation Error : "<<stringBuffer<< std::endl; 
		return false;
	}
	
  /*
	glGetShaderiv(pShaderID, GL_COMPILE_STATUS, &ok);
	if (!ok){
		int ilength; char stringBuffer[STRING_BUFFER];
		glGetShaderInfoLog(pShaderID, STRING_BUFFER, &ilength, stringBuffer);
		std::cerr << "[ERROR] Shader Compilation Error ("<<DepthComplexity2D::P_PATH<<") : "<<stringBuffer << std::endl; 
		return false;
	}
	*/
  
	std::cerr << "[OK] Shaders OK." << std::endl;
	return true;
}

bool DepthComplexity2D::initFBO() {
  
  // Use a texture as COLOR buffer for FBO -----------------------------
  glGenTextures(1, &_cboTexId);
  glBindTexture(GL_TEXTURE_2D, _cboTexId);
  // Set texture filters
  glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
  glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
  glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
  glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
  
  glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA8, _fboWidth, _fboHeight, 0, GL_RGBA, GL_UNSIGNED_BYTE , NULL);
  glBindTexture(GL_TEXTURE_2D, 0);

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

  // User renderbuffer as DEPTH buffer ---------------------------------
  glGenRenderbuffersEXT(1, &_rboId);
  glBindRenderbufferEXT(GL_RENDERBUFFER_EXT, _rboId);
  glRenderbufferStorageEXT(GL_RENDERBUFFER_EXT, GL_DEPTH_COMPONENT, _fboWidth, _fboHeight);
  glBindRenderbufferEXT(GL_RENDERBUFFER_EXT, 0);

  // Create framebuffer object -----------------------------------------
  glGenFramebuffersEXT(1, &_fboId);
  glBindFramebufferEXT(GL_FRAMEBUFFER_EXT, _fboId);

  // Assign COLOR1, COLOR2 and DEPTH buffers to the render buffer.
  glFramebufferTexture2DEXT(GL_FRAMEBUFFER_EXT, GL_COLOR_ATTACHMENT0_EXT, GL_TEXTURE_2D, _cboTexId, 0);
  glFramebufferRenderbufferEXT(GL_FRAMEBUFFER_EXT, GL_DEPTH_ATTACHMENT_EXT, GL_RENDERBUFFER_EXT, _rboId);

  // Check if everything is ok.
  _status = checkFramebufferStatus();

  glBindFramebufferEXT(GL_FRAMEBUFFER_EXT, 0);
  
  return true;
}

bool checkFramebufferStatus() {
  // check FBO status
  GLenum status = glCheckFramebufferStatusEXT(GL_FRAMEBUFFER_EXT);
  switch(status) {
  case GL_FRAMEBUFFER_COMPLETE_EXT:
    std::cerr << "[OK] FBO Complete." << std::endl;
    return true;

  case GL_FRAMEBUFFER_INCOMPLETE_ATTACHMENT_EXT:
    std::cerr << "[ERROR] FBO incomplete: Attachment is NOT complete." << std::endl;
    return false;

  case GL_FRAMEBUFFER_INCOMPLETE_MISSING_ATTACHMENT_EXT:
    std::cerr << "[ERROR] FBO incomplete: No image is attached to FBO." << std::endl;
    return false;

  case GL_FRAMEBUFFER_INCOMPLETE_DIMENSIONS_EXT:
    std::cerr << "[ERROR] FBO incomplete: Attached images have different dimensions." << std::endl;
    return false;

  case GL_FRAMEBUFFER_INCOMPLETE_FORMATS_EXT:
    std::cerr << "[ERROR] FBO incomplete: Color attached images have different internal formats." << std::endl;
    return false;

  case GL_FRAMEBUFFER_INCOMPLETE_DRAW_BUFFER_EXT:
    std::cerr << "[ERROR] FBO incomplete: Draw buffer." << std::endl;
    return false;

  case GL_FRAMEBUFFER_INCOMPLETE_READ_BUFFER_EXT:
    std::cerr << "[ERROR] FBO incomplete: Read buffer." << std::endl;
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
  
  if (!_status){
	  std::cerr << "[ERROR] Check FBO or SHADERS." << std::endl;
	  return;
  }
  //for (int i=0; i<_fboHeight*_fboWidth; ++i) { _bufferCounter[i] = 0;}
  //glPixelStorei(GL_UNPACK_ALIGNMENT, 1);  
  
  findDepthComplexity2D();
  
  if (_computeHistogram or _computeMaximumRays or _computeGoodRays)
    findMaximumRaysAndHistogram();
  
  //copyStencilToColor();
}


void DepthComplexity2D::findDepthComplexity2D() {
  glBindFramebufferEXT(GL_FRAMEBUFFER_EXT, _fboId);  
  
  glPushAttrib(GL_VIEWPORT_BIT | GL_ENABLE_BIT);
    glViewport(0, 0, _fboWidth, _fboHeight);
	


				glDisable(GL_POINT_SMOOTH);
	   glDisable(GL_LINE_SMOOTH);
				glDisable(GL_POLYGON_SMOOTH);
				
				glHint(GL_POLYGON_SMOOTH_HINT,GL_FASTEST);
				glHint(GL_POINT_SMOOTH, GL_FASTEST);
				glHint(GL_LINE_SMOOTH, GL_FASTEST);
				
	   glDisable(GL_DEPTH_TEST);
		
    glMatrixMode(GL_PROJECTION);
    glPushMatrix();
      glLoadIdentity();
      gluOrtho2D(0.0, 1.0, 0.0, 1.0);
      
      glMatrixMode(GL_MODELVIEW);
		glPushMatrix();
		glLoadIdentity();

		glActiveTexture(GL_TEXTURE0);
		glEnable(GL_TEXTURE_2D);
		glBindTexture(GL_TEXTURE_2D, _counterBuffId);   
		
    
  // CLEAR ======================================================================================================================
		//Use Shader
		glUseProgram(_shaderProgram2);
    
		//Pass MODELVIEW and PROJECTION matrices to Shader
		GLfloat model[16],projection[16];
		glGetFloatv(GL_MODELVIEW_MATRIX, model);
		glGetFloatv(GL_PROJECTION_MATRIX, projection);
		glProgramUniformMatrix4fv(_shaderProgram2, glGetUniformLocation(_shaderProgram2, "modelViewMat"), 1, GL_FALSE, model);
		glProgramUniformMatrix4fv(_shaderProgram2, glGetUniformLocation(_shaderProgram2, "projectionMat"), 1, GL_FALSE, projection);
    
    //Pass counter buff texture
		glProgramUniform1iEXT(_shaderProgram2, glGetUniformLocation(_shaderProgram2, "counterBuff"), 0);
		  
		glClearColor(0.f, 0.f, 0.f, 0.f);    
		glClear(GL_COLOR_BUFFER_BIT);    //GL_DEPTH_BUFFER_BIT  | 

    
    float myVertexArray[8] = {0.0f, 1.0f, 
                              1.0f, 1.0f,
                              1.0f, 0.0f,
                              0.0f, 0.0f};
                              
    // first get the attribute-location
    GLint vertexLocation = glGetAttribLocation(_shaderProgram2, "vertexPos");
    // enable an array for the attribute
    glEnableVertexAttribArray(vertexLocation);
    // set attribute pointer
    glVertexAttribPointer(vertexLocation, 2, GL_FLOAT,GL_FALSE, 0, myVertexArray);
    // Draw ("render") the triangle
    glDrawArrays(GL_QUADS, 0, 4);
    // Done with rendering. Disable vertex attribute array
    glDisableVertexAttribArray(vertexLocation);
    
    
     glMemoryBarrierEXT(GL_TEXTURE_UPDATE_BARRIER_BIT_EXT);
    // ======================================================================================================================
    glUseProgram(_shaderProgram);
   
	
		//Pass MODELVIEW and PROJECTION matrices to Shader
		//GLfloat model[16],projection[16];
		//glGetFloatv(GL_MODELVIEW_MATRIX, model);
		//glGetFloatv(GL_PROJECTION_MATRIX, projection);
		glProgramUniformMatrix4fv(_shaderProgram, glGetUniformLocation(_shaderProgram, "modelViewMat"), 1, GL_FALSE, model);
		glProgramUniformMatrix4fv(_shaderProgram, glGetUniformLocation(_shaderProgram, "projectionMat"), 1, GL_FALSE, projection);
      
    //std::cout << "glProjection:\n" << projection[0] << " - " << projection[5] << " - " << projection[10] << std::endl;
    
    //Pass counter buff texture
		glProgramUniform1iEXT(_shaderProgram, glGetUniformLocation(_shaderProgram, "counterBuff"), 0);
		   
		//glClearColor(0.f, 0.f, 0.f, 0.f);    
		//glClear(GL_COLOR_BUFFER_BIT);    //GL_DEPTH_BUFFER_BIT  | 

    /*
    float myVertexArray2[2] = {0.5f,0.5f};
    // first get the attribute-location
    vertexLocation = glGetAttribLocation(_shaderProgram, "vertexPos");
    // enable an array for the attribute
    glEnableVertexAttribArray(vertexLocation);
    // set attribute pointer
    glVertexAttribPointer(vertexLocation, 2, GL_FLOAT,GL_FALSE, 0, myVertexArray2);
    // Draw ("render") the triangle
    glDrawArrays(GL_POINTS, 0, 1);
    // Done with rendering. Disable vertex attribute array
    glDisableVertexAttribArray(vertexLocation);
    
    */
  
    //glBegin(GL_POINTS);
    //  glVertex2f(0.5f, 0.5f);
    //glEnd();
    /*
    glBegin(GL_QUADS);
      glVertex2f(-1.0f, 1.0f);
      glVertex2f(1.0f, 1.0f);
      glVertex2f(1.0f, -1.0f);
      glVertex2f(-1.0f, -1.0f);
    glEnd();
    */
    
		double pixelSize = 1.0/512.0;
		double step = 0.5*pixelSize*sqrt(2); //1.41421356;
		for (unsigned i=0; i<_segments->size(); ++i) {
		  if (!_segments->at(i).active) continue;
		
		  Segment Tv0 = computeDualSegmentFromPoint(_segments->at(i).a);
		  Segment Tv1 = computeDualSegmentFromPoint(_segments->at(i).b);

		  if (Tv0.a.y < Tv1.a.y) std::swap(Tv0, Tv1);

		  double t1, t2;
		  if (segmentIntersection2D(Tv0, Tv1, &t1, &t2) ) {			
						double t1, t2;
						
						// Triangle 1
						if (lineIntersection2D(Tv0.a, Tv1.a, &t1, &t2)) {
								// find intersection point
								vec3d inter = Tv0.a*(1.0-t1) + Tv0.b*t1;

								// find median points of the triangle sides
								vec3d pm1 = (Tv0.a + inter)*0.5; // for Tv1
								vec3d pm2 = (inter + Tv1.a)*0.5; // for Tv0
								vec3d pm3 = (Tv1.a + Tv0.a)*0.5; // for inter
								
								// find inward directions for each vertex, to shrink triangle;
								vec3d u = pm1 - Tv1.a; u.normalize();
								vec3d v = pm2 - Tv0.a; v.normalize();
								vec3d t = pm3 - inter; t.normalize();
	
								// Triangle 1
								Triangle tri;
								tri.a = Tv0.a + v*step;
								tri.b = Tv1.a + u*step;
								tri.c = inter + t*step;

								if (tri.a.y > tri.b.y) {
										drawTriangle(tri.a, tri.c, tri.b);
								}
						}
						
						// Triangle 2
						if (lineIntersection2D(Tv0.b, Tv1.b, &t1, &t2) ) {
								vec3d inter = Tv1.a*(1.f-t1) + Tv1.b*t1;
						
								// find median points of the triangle sides
								vec3d pm1 = (Tv0.b + inter)*0.5; // for Tv1
								vec3d pm2 = (inter + Tv1.b)*0.5; // for Tv0
								vec3d pm3 = (Tv1.b + Tv0.b)*0.5; // for inter

							 // find inward directions for each vertex, to shrink triangle;
								vec3d u = pm1 - Tv1.b; u.normalize();
								vec3d v = pm2 - Tv0.b; v.normalize();
								vec3d t = pm3 - inter; t.normalise();

								// Triangle 1
								Triangle tri;
								tri.a = Tv0.b + v*step;
								tri.b = Tv1.b + u*step;
								tri.c = inter + t*step;
									
								if (tri.a.y < tri.b.y) {								
											drawTriangle(tri.a, tri.c, tri.b);
						  }   
						}
   } else {
												// u,v inward directions to shrink trapezoid
												vec3d u = Tv1.b - Tv0.a; u.normalize();
												vec3d v = Tv1.a - Tv0.b; v.normalize();
												
												// increment vertices of trapezoid inward
            Segment s1(Tv0.a + u*step, Tv0.b + v*step);
            Segment s2(Tv1.a - v*step, Tv1.b - u*step);

            if (s1.a.y > s2.a.y and s1.b.y > s2.b.y) {
              drawQuad(s1.a, s1.b, s2.b, s2.a);
            }
          }
        }
		    
        //Ensure that all texture writing is done
        glMemoryBarrierEXT(GL_TEXTURE_UPDATE_BARRIER_BIT_EXT);
        glUseProgram(0);
	
        _maximum = findMaxValueInCounterBuffer()-1;
        glBindTexture(GL_TEXTURE_2D, 0);   

	    glMatrixMode(GL_MODELVIEW);
	    glPopMatrix();
     
	  glMatrixMode(GL_PROJECTION);
	glPopMatrix();
	  
	  glEnable(GL_DEPTH_TEST);
    glPopAttrib();
  glBindFramebufferEXT(GL_FRAMEBUFFER_EXT, 0);
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

unsigned int DepthComplexity2D::findMaxValueInCounterBuffer() {
  const int pixelNumber = _fboWidth * _fboHeight;
  unsigned int colorBuffer[pixelNumber];
  
  glBindTexture(GL_TEXTURE_2D, _counterBuffId);
  glGetTexImage( GL_TEXTURE_2D, 0 , GL_RED_INTEGER, GL_UNSIGNED_INT, colorBuffer ); 
  glBindTexture(GL_TEXTURE_2D, 0);
  
  for (int i=0; i<pixelNumber; ++i){
    std::cout << colorBuffer[i]-1 << std::endl;
  }
  
  return *(std::max_element(colorBuffer, colorBuffer + pixelNumber));
}

void DepthComplexity2D::findMaximumRaysAndHistogram() {
  _maximumRays.clear();
  _goodRays.clear();
  _goodRays.resize(_maximum + 1);
  _histogram.clear();
  _histogram.resize(_maximum + 1);

    glBindFramebufferEXT(GL_FRAMEBUFFER_EXT, _fboId);
    glPushAttrib(GL_VIEWPORT_BIT);
    glViewport(0.0, 0.0, _fboWidth, _fboHeight);
    
    unsigned int colorBuffer[_fboWidth*_fboHeight];
    
    glBindTexture(GL_TEXTURE_2D, _counterBuffId);
    glGetTexImage( GL_TEXTURE_2D, 0 , GL_RED_INTEGER, GL_UNSIGNED_INT, colorBuffer ); 
    glBindTexture(GL_TEXTURE_2D, 0);
    
    for(int r=0; r<_fboHeight; ++r) {
      for(int c=0; c<_fboWidth; ++c) {
        unsigned int val = colorBuffer[r*_fboWidth+c]-1;

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
			   
          if (_goodRays[val].size() < 200)
            _goodRays[val].insert(seg);
            
		  }
        }
      }
    }
  glPopAttrib();
  glBindFramebufferEXT(GL_FRAMEBUFFER_EXT, 0);
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

  glPushAttrib(GL_VIEWPORT_BIT | GL_ENABLE_BIT);
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
