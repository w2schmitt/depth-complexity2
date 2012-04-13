#include "dc_2d.h"

#include <iostream>
#include <cassert>

const char *DepthComplexity2D::V_PATH = "shader/dc.vert";
const char *DepthComplexity2D::P_PATH = "shader/dc.frag";

// Internal function to check if the framebuffer was created correctly.
bool checkFramebufferStatus();
// Internal function to check if shader was compiled correctly.
bool checkShadersStatus(GLuint vShaderID, GLuint pShaderID);

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
  _vShaderId = -1;
  _pShaderId = -1;
  _shaderProgram = -1;
  _status = true;

  assert(initFBO());
  createShader();
  attachShader();  
}

void DepthComplexity2D::createShader(){
	
	_vShaderId = glCreateShader(GL_VERTEX_SHADER);
	_pShaderId = glCreateShader(GL_FRAGMENT_SHADER);
	
	// Read shader files
	char *vShaderData = readShaderFile(V_PATH);
	char *pShaderData = readShaderFile(P_PATH);
	
	// Pointer to shader data
	const char *ptVData = vShaderData;
	const char *ptPData = pShaderData;

	glShaderSource(_vShaderId,  1, &ptVData, NULL);
	glShaderSource(_pShaderId, 1, &ptPData, NULL);
	
	glCompileShader(_vShaderId);
	glCompileShader(_pShaderId);
	
	//Check compilation status
	_status = checkShadersStatus(_vShaderId, _pShaderId);

	free(vShaderData);
	free(pShaderData);
}

void DepthComplexity2D::attachShader(){
	
	_shaderProgram = glCreateProgram();
	
	glAttachShader(_shaderProgram, _vShaderId);
	glAttachShader(_shaderProgram, _pShaderId);
	
	glLinkProgram(_shaderProgram);
}


char* DepthComplexity2D::readShaderFile(const char *filename){
	
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

bool checkShadersStatus(GLuint vShaderID, GLuint pShaderID){
	
	GLint ok;
	glGetShaderiv(vShaderID, GL_COMPILE_STATUS, &ok);
	if (!ok){
		int ilength; char stringBuffer[STRING_BUFFER];
		glGetShaderInfoLog(vShaderID, STRING_BUFFER, &ilength, stringBuffer);
		std::cerr<<"[ERROR] Compilation error ("<<DepthComplexity2D::V_PATH<<") : "<<stringBuffer<< std::endl; 
		return false;
	}
	
	glGetShaderiv(pShaderID, GL_COMPILE_STATUS, &ok);
	if (!ok){
		int ilength; char stringBuffer[STRING_BUFFER];
		glGetShaderInfoLog(pShaderID, STRING_BUFFER, &ilength, stringBuffer);
		std::cerr << "[ERROR] Compilation error ("<<DepthComplexity2D::P_PATH<<") : "<<stringBuffer << std::endl; 
		return false;
	}
	
	std::cerr << "[OK] Shaders OK." << std::endl;
	return true;
}

bool DepthComplexity2D::initFBO() {
  // Use texture 1 as COLOR buffer for FBO
  glGenTextures(1, &_textureId);
  glBindTexture(GL_TEXTURE_2D, _textureId);
  glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
  glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
  glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
  glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
  
  glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA8, _fboWidth, _fboHeight, 0, GL_RGBA, GL_UNSIGNED_BYTE, NULL);
  glBindTexture(GL_TEXTURE_2D, 0);
    
  //Use texture2 as a counter per-pixel Buffer//
  glGenTextures(1, &_textureId2);
  glBindTexture(GL_TEXTURE_2D, _textureId2);
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
  glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
  glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);

  //Uses GL_R32F instead of GL_R32I that is not working in R257.15
  glTexImage2D(GL_TEXTURE_2D, 0, GL_R32F, _fboWidth, _fboHeight, 0,  GL_RED, GL_FLOAT, 0);
  glBindImageTextureEXT(1, _textureId2, 0, false, 0,  GL_READ_WRITE, GL_R32UI);
  glBindTexture(GL_TEXTURE_2D, 0);
  

  /*
  glGenTextures(1, &_textureId2);
  glBindTexture(GL_TEXTURE_2D, _textureId2);
  glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
  glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
  glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
  glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
  glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA16, _fboWidth, _fboHeight, 0, GL_RGBA, GL_UNSIGNED_SHORT, NULL);
  glBindTexture(GL_TEXTURE_2D, 0);
	*/
  // User renderbuffer as DEPTH buffer
  glGenRenderbuffersEXT(1, &_rboId);
  glBindRenderbufferEXT(GL_RENDERBUFFER_EXT, _rboId);
  glRenderbufferStorageEXT(GL_RENDERBUFFER_EXT, GL_DEPTH_COMPONENT, _fboWidth, _fboHeight);
  glBindRenderbufferEXT(GL_RENDERBUFFER_EXT, 0);

  // Create framebuffer object
  glGenFramebuffersEXT(1, &_fboId);
  glBindFramebufferEXT(GL_FRAMEBUFFER_EXT, _fboId);

  // Assign COLOR1, COLOR2 and DEPTH buffers to the render buffer.
  glFramebufferTexture2DEXT(GL_FRAMEBUFFER_EXT, GL_COLOR_ATTACHMENT0_EXT, GL_TEXTURE_2D, _textureId, 0);
  //glFramebufferTexture2DEXT(GL_FRAMEBUFFER_EXT, GL_COLOR_ATTACHMENT1_EXT, GL_TEXTURE_2D, _textureId2, 0);
  glFramebufferRenderbufferEXT(GL_FRAMEBUFFER_EXT, GL_DEPTH_ATTACHMENT_EXT, GL_RENDERBUFFER_EXT, _rboId);
  //glFramebufferRenderbufferEXT(GL_FRAMEBUFFER_EXT, GL_STENCIL_ATTACHMENT_EXT, GL_RENDERBUFFER_EXT, _rboId);

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
    std::cerr << "[OK] Framebuffer complete." << std::endl;
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
  
  std::cout << "--------------------------------------------------------< VERSION >\n";
  std::cout << "OGL VERSION: " << glGetString(GL_VERSION) << std::endl;
  std::cout << "SHADER VERSION: " << glGetString(GL_SHADING_LANGUAGE_VERSION) << std::endl;
  std::cout << "--------------------------------------------------------- \n";
  
  if (!_status){
	  std::cerr << "[ERROR] Check FBO or SHADERS." << std::endl;
	  return;
  }
  
  glPixelStorei(GL_UNPACK_ALIGNMENT, 1);  
  findDepthComplexity2D();
  
  if (_computeHistogram or _computeMaximumRays or _computeGoodRays)
    findMaximumRaysAndHistogram();
  
  //copyStencilToColor();
}


void DepthComplexity2D::findDepthComplexity2D() {
  glBindFramebufferEXT(GL_FRAMEBUFFER_EXT, _fboId);
  
  
  // Disable LIGHT, BLENDING and TEXTURE
  glDisable(GL_LIGHT0);
  glDisable(GL_LIGHTING);
  glDisable(GL_BLEND);
  //glDisable(GL_TEXTURE_2D);
  glDisable(GL_LINE_SMOOTH);
  glDisable(GL_DEPTH_TEST);
    
  glPushAttrib(GL_VIEWPORT_BIT);
    glViewport(0, 0, _fboWidth, _fboHeight);
		
    glMatrixMode(GL_PROJECTION);
    glPushMatrix();
      glLoadIdentity();
      gluOrtho2D(0., 1.0, 0., 1.0);
      
      glMatrixMode(GL_MODELVIEW);
	  glPushMatrix();
		glLoadIdentity();
		
		
		  //Pass MODELVIEW and PROJECTION matrices to Shader
		  GLfloat model[16],projection[16];
		  glGetFloatv(GL_MODELVIEW_MATRIX, model);
		  glGetFloatv(GL_PROJECTION_MATRIX, projection);
		  glProgramUniformMatrix4fvEXT(_shaderProgram, glGetUniformLocation(_shaderProgram, "modelViewMat"), 1, GL_FALSE, model);
		  glProgramUniformMatrix4fvEXT(_shaderProgram, glGetUniformLocation(_shaderProgram, "projectionMat"), 1, GL_FALSE, projection);
		  
		  
		  //Pass texture2
		  glProgramUniform1iEXT(_shaderProgram, glGetUniformLocation(_shaderProgram, "counterBuff"), 0);
		  //Use Shader
		  glUseProgram(_shaderProgram);
		  // Select both color buffer to clear
		  //GLenum buffers[] = {GL_COLOR_ATTACHMENT0_EXT, GL_COLOR_ATTACHMENT1_EXT};
		  //glDrawBuffers(2, buffers);
		  
		   glClearColor(0.f, 0.f, 0.f, 0.f);    
		   glClear(GL_DEPTH_BUFFER_BIT  | GL_COLOR_BUFFER_BIT); 

		   glDrawBuffer(GL_COLOR_ATTACHMENT0_EXT);
	
		   

		   glColor4f(0.5f, 0.2f, 0.1f, 0.15f);
		   //
		   glBegin(GL_QUADS);
			 glVertex2f(0.0f, 0.0f);
			 glVertex2f(1.0f, 0.0f);
			 glVertex2f(1.0f, 1.0f);
			 glVertex2f(0.0f, 1.0f);
		   glEnd();
		   
		   //glFinish();
		   		   
		   //glDrawBuffer(GL_COLOR_ATTACHMENT1_EXT);
		  
		   //glActiveTexture(GL_TEXTURE0);		   
		   //glBindTexture(GL_TEXTURE_2D, _textureId2);
		   //glEnable(GL_TEXTURE_2D);
		   //glUniform1i(glGetUniformLocation(shaderProgram, "dcCounter"), 0);
		    
		    //
			glFramebufferTexture2DEXT(GL_FRAMEBUFFER_EXT, GL_COLOR_ATTACHMENT0_EXT, GL_TEXTURE_2D, _textureId2, 0);		
	
			glReadBuffer(GL_COLOR_ATTACHMENT0_EXT);
			
		    const int pixelNumber = _fboWidth * _fboHeight;
		    //unsigned short int colorBuffer[pixelNumber];	
		    GLushort colorBuffer[pixelNumber];	  
	
  
			glReadPixels(0, 0, _fboWidth, _fboHeight, GL_RED, GL_UNSIGNED_SHORT, colorBuffer);

			std::cout << "-MAX:"<< (int)*(std::max_element(colorBuffer, colorBuffer + pixelNumber)) << "\n";

		  
	  
      

	  glMatrixMode(GL_MODELVIEW);
	  glPopMatrix();
     
	  glMatrixMode(GL_PROJECTION);
	  glPopMatrix();
	  
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

unsigned short int DepthComplexity2D::findMaxValueInStencil() {
  const int pixelNumber = _fboWidth * _fboHeight;
  GLushort colorBuffer[pixelNumber];
  

  //glTexImage2D(GL_TEXTURE_2D, 0, GL_LUMINANCE32UI_EXT, _fboWidth, _fboHeight, 0, GL_LUMINANCE_INTEGER_EXT, GL_UNSIGNED_INT, NULL);
  glReadPixels(0, 0, _fboWidth, _fboHeight, GL_LUMINANCE, GL_UNSIGNED_SHORT, colorBuffer);
  unsigned short int val;
  val = *(std::max_element(colorBuffer, colorBuffer + pixelNumber));
  std::cout << "-MAX:"<< val << "\n";
  return val;
 
  
  
}

void DepthComplexity2D::findMaximumRaysAndHistogram() {
	
  _maximumRays.clear();
  _goodRays.clear();
  _goodRays.resize(_maximum + 1);
  _histogram.clear();
  _histogram.resize(_maximum + 1);
  
  //glPixelTransferi(GL_MAP_STENCIL, GL_TRUE);

  glBindFramebufferEXT(GL_FRAMEBUFFER_EXT, _fboId);
  glPushAttrib(GL_VIEWPORT_BIT);
    glViewport(0.0, 0.0, _fboWidth, _fboHeight);
    
	
    GLushort colorBuffer[_fboWidth*_fboHeight];
    glReadPixels(0, 0, _fboWidth, _fboHeight, GL_RED, GL_UNSIGNED_SHORT, colorBuffer);
    
    //printf("stencil: \n");
    //unsigned value = *(std::max_element(stencilSave, stencilSave + _fboWidth*_fboHeight));
    //std::cout << "maxxx: " << value << "\n";
	//std::cout << "height: " << _fboHeight << "\n";
    for(int r=0; r<_fboHeight; ++r) {
      for(int c=0; c<_fboWidth; ++c) {
        //GLubyte val = stencilSave[r*_fboWidth+c];
       
        //printf("val:%u\n", colorBuffer[r*_fboWidth+c]);
        unsigned val = colorBuffer[r*_fboWidth+c];
        //std::cout << "cr: "<< r*_fboWidth+c << "\n";
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
		  //std::cout << "fuck2\n";
        }
      }
    }
  std::cout << "fuck yeah\n";
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
