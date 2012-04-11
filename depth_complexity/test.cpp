
#include <GL/glew.h>
#include <GL/glut.h>


#include <iostream>
#include <algorithm>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <vector>
#include <string.h>

using namespace std;

/*
#define ERR 1.E-2

struct Point { float x,y; };
struct Triangle { Point a,b,c; };
struct Trapezoid { Point a,b,c,d; };

bool myfn(float i, float j) { return i<j; }

//---------------------------------- GLOBALS
int Nintersections = 0;

int fbo_w = 512;
int fbo_h = 512;
int fbo_c = 1;
float *fbo,*tmpFBO;

std::vector< Triangle > vecTri1;
std::vector< Triangle > vecTri2;
std::vector< Trapezoid > vecTrap;

bool normalVisualization = true;
bool displayTri1 = true;
bool displayTri2 = true;
bool displayTrap = true;

bool readFbo = false;




void readInput(){
	
	int t=0;
	bool terminate = true;
	Triangle tri;
	Trapezoid trap;
	
	while (terminate){
		cin >> t;

		switch (t){
			case 1:
				cin >> tri.a.x >> tri.a.y >> tri.b.x >> tri.b.y >> tri.c.x >> tri.c.y;
				vecTri1.push_back(tri);
				break;
			case 2:
				cin >> tri.a.x >> tri.a.y >> tri.b.x >> tri.b.y >> tri.c.x >> tri.c.y;
				vecTri2.push_back(tri);
				break;
			case 3:
				cin >> trap.a.x >> trap.a.y >> trap.b.x >> trap.b.y >> trap.c.x >> trap.c.y >> trap.d.x >> trap.d.y;
				vecTrap.push_back(trap);
				break;
			case 4:
			   //terminate = false;
			   break;
			case 5:
				terminate = false;
				break;
		}		
	}
}

void readFboInput(){
	float t;
	int i=0;
	
	while (cin >> t){
		fbo[i] = t;		
		//if (fbo[i]!=0)
		//	std::cout<<fbo[i]<<std::endl;
		i++;
	}
}


float abs1(float x){
	if (signbit(x)!= 0){
		return x*-1.;
	}
	return x;
	
}

void intersectionBuffer(float*  buff, int n, int w, int h, int ch){
	for (int i=0; i< w*h*ch; ++i){
		if (abs1(float(buff[i] - (float)n) ) <= ERR)
			buff[i] = 1.;
		else
			buff[i] = 0.;
	}
}

void clearBuffer(float* buff, int w, int h, int ch){
	for (int i=0; i< w*h*ch; ++i){ buff[i] = 0.;}
}

// increment buff1 with values from buff2
void sumBuffer(float* buff1, float* buff2, int w, int h, int ch){
	for (int i=0; i< w*h*ch; ++i){ buff1[i] += buff2[i];}
}

void copyBuffer(float* buff1, float* buff2, int w, int h, int ch){
	for (int i=0; i< w*h*ch; ++i){ buff1[i] = buff2[i];}
}


void setupRC(){
	glEnable(GL_COLOR_MATERIAL);
	glColorMaterial(GL_FRONT, GL_AMBIENT_AND_DIFFUSE);
	glClearColor(0, 0, 0, 1.0);
	glDisable(GL_LINE_SMOOTH);
	
	glEnable(GL_BLEND);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
    
    //initialize fbo
    if ((fbo = new float[fbo_w*fbo_h])){
		clearBuffer(fbo,fbo_w,fbo_h,fbo_c);
	}
	else
		std::clog << "ERROR: cannot initilize FBO!";
		
	if ((tmpFBO = new float[fbo_w*fbo_h])){
		clearBuffer(tmpFBO,fbo_w,fbo_h,fbo_c);	
	}
	else
		std::clog << "ERROR: cannot initilize tmpFBO!";
		
	//std::cout << ERR << "\n";
}

void reshape(int w, int h) {
  //	float aspect = w/h;	
  if(h==0) h = 1;
  glMatrixMode(GL_PROJECTION);
  glLoadIdentity();	
  glViewport(0, 0, w, h);
  //glOrtho(-1, 2, -1, 2, -1.0, 1.0);
  //glOrtho(-0.6, 1.6, -0.4, 2.2, -1.0, 1.0);
  gluOrtho2D(0, 1, 0, 1);
  glMatrixMode(GL_MODELVIEW);
}

void display(){

	
	
	if (readFbo){
		glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
		
		copyBuffer(tmpFBO, fbo, fbo_w, fbo_h, fbo_c);
		intersectionBuffer(tmpFBO, Nintersections, fbo_w, fbo_h, fbo_c);
		
		//for (int i=0; i<fbo_w*fbo_h; ++i) std::cout << fbo[i] << std::endl;
		
		glDrawPixels(fbo_w,fbo_h, GL_LUMINANCE, GL_FLOAT, tmpFBO );	
	}	
	else if (normalVisualization){
		clearBuffer(fbo, fbo_w, fbo_h, fbo_c);
		glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
		
		if (displayTri1){
			glColor3f(1,0,0);
			glBegin(GL_TRIANGLES);	
			for (unsigned int i=0; i < vecTri1.size(); ++i)
			{
				glVertex2f(vecTri1[i].a.x,vecTri1[i].a.y);
				glVertex2f(vecTri1[i].b.x,vecTri1[i].b.y);
				glVertex2f(vecTri1[i].c.x,vecTri1[i].c.y);
			}
			glEnd();
		}
		
		if (displayTri2){	
			glColor3f(0,1,0);
			glBegin(GL_TRIANGLES);	
			for (unsigned int i=0; i < vecTri2.size(); ++i)
			{
				glVertex2f(vecTri2[i].a.x,vecTri2[i].a.y);
				glVertex2f(vecTri2[i].b.x,vecTri2[i].b.y);
				glVertex2f(vecTri2[i].c.x,vecTri2[i].c.y);		
			}		
			glEnd();
		}
		
		if (displayTrap){
			glColor3f(0,0,1);
			glBegin(GL_QUADS);
			for (unsigned int i=0; i < vecTrap.size(); ++i)
			{
				glVertex2f(vecTrap[i].a.x,vecTrap[i].a.y);
				glVertex2f(vecTrap[i].b.x,vecTrap[i].b.y);
				glVertex2f(vecTrap[i].c.x,vecTrap[i].c.y);
				glVertex2f(vecTrap[i].d.x,vecTrap[i].d.y);				
			}
			glEnd();
		}
		
					
	}
	else {
		clearBuffer(fbo, fbo_w, fbo_h, fbo_c);
		if (displayTri1){
			glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
			glColor3f(1,0,0);
			for (unsigned int i=0; i < vecTri1.size(); ++i){
				glBegin(GL_TRIANGLES);			
					glVertex2f(vecTri1[i].a.x,vecTri1[i].a.y);
					glVertex2f(vecTri1[i].b.x,vecTri1[i].b.y);
					glVertex2f(vecTri1[i].c.x,vecTri1[i].c.y);
				glEnd();
			}		
			glReadPixels(0,0,fbo_w,fbo_h, GL_LUMINANCE, GL_FLOAT, tmpFBO);
			sumBuffer(fbo,tmpFBO,fbo_w,fbo_h,1);
		}
	
		if (displayTri2){
			glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
			glColor3f(0,1,0);
			for (unsigned int i=0; i < vecTri2.size(); ++i){
				glBegin(GL_TRIANGLES);				
					glVertex2f(vecTri2[i].a.x,vecTri2[i].a.y);
					glVertex2f(vecTri2[i].b.x,vecTri2[i].b.y);
					glVertex2f(vecTri2[i].c.x,vecTri2[i].c.y);				
				glEnd();
			}
			glReadPixels(0,0,fbo_w,fbo_h, GL_LUMINANCE, GL_FLOAT, tmpFBO);
			sumBuffer(fbo,tmpFBO,fbo_w,fbo_h,1);
		}
				
		if (displayTrap){
			glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
			glColor3f(0,0,1);
			for (unsigned int i=0; i < vecTrap.size(); ++i){
				glBegin(GL_QUADS);			
					glVertex2f(vecTrap[i].a.x,vecTrap[i].a.y);
					glVertex2f(vecTrap[i].b.x,vecTrap[i].b.y);
					glVertex2f(vecTrap[i].c.x,vecTrap[i].c.y);
					glVertex2f(vecTrap[i].d.x,vecTrap[i].d.y);			
				glEnd();
			}
			glReadPixels(0,0,fbo_w,fbo_h, GL_LUMINANCE, GL_FLOAT, tmpFBO);
			sumBuffer(fbo,tmpFBO,fbo_w,fbo_h,1);
		}
		
		glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
		//std::cout << "max depth: " << *(std::max_element(fbo, fbo + fbo_w*fbo_h*fbo_c)) << "\n";
		
		intersectionBuffer(fbo, Nintersections, fbo_w, fbo_h, fbo_c);	
		glDrawPixels(fbo_w,fbo_h, GL_LUMINANCE, GL_FLOAT, fbo );
	}
	
	glutSwapBuffers();
	glFinish();
}




void keyboard(unsigned char key, int x, int y){
	
	if (key=='a')
		displayTri1 = !displayTri1;
	else if (key=='s')
		displayTri2 = !displayTri2;
	else if (key=='d')
		displayTrap = !displayTrap;
	else if (key=='q')
		normalVisualization = !normalVisualization;
	else if (key >= '0' && key <= '0'+9)
		Nintersections = (int)(key-'0');
		
	glutPostRedisplay();
}


int main (int argc, char **argv) {
	//std::cout << argv[1] << "\n";
	if (argc==2){
		if (strcmp(argv[1],"-sbo")==0)
			readFbo=true;
	}
	glutInit(&argc, argv);
	glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGBA | GLUT_DEPTH);
	glutInitWindowSize(512, 512);
	glutInitWindowPosition(300,300);
	glutCreateWindow("Dual Ray Implementation");
	glutDisplayFunc(display);	

	setupRC();
	
	if (readFbo)
		readFboInput();
	else
		readInput();
	
	
	glutReshapeFunc(reshape);
	glutDisplayFunc(display);
	glutKeyboardFunc( keyboard );
	//glutSpecialFunc( special );

	glutMainLoop();
	
	return 0;
}
* 
*/

  GLuint                                		_textureId;
  GLuint                                		_fboId;
  GLuint                                		_rboId;

int _fboWidth = 512;
int _fboHeight = 512;

bool initFBO() {
  // Use texture as COLOR buffer
  glGenTextures(1, &_textureId);
  glBindTexture(GL_TEXTURE_2D, _textureId);
  glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
  glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
  glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
  glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
  //glTexImage2D(GL_TEXTURE_2D, 0, GL_LUMINANCE16, _fboWidth, _fboHeight, 0, GL_LUMINANCE, GL_UNSIGNED_SHORT, NULL); 
  glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA16, _fboWidth, _fboHeight, 0, GL_RGBA, GL_FLOAT, NULL);
	
	GLint actual_format;
	glGetTexLevelParameteriv(GL_TEXTURE_2D, 0, GL_TEXTURE_INTERNAL_FORMAT, &actual_format);
	if (actual_format == GL_RGBA16){
	  std::cout << "format OK\n";
	}
	else {
	  std::cout << "teste: "<< actual_format << "  -  " << GL_RGBA16 << "\n";
	}
  
  glBindTexture(GL_TEXTURE_2D, 0);

  // User renderbuffer as DEPTH and STENCIL buffer.
  glGenRenderbuffersEXT(1, &_rboId);
  glBindRenderbufferEXT(GL_RENDERBUFFER_EXT, _rboId);
  glRenderbufferStorageEXT(GL_RENDERBUFFER_EXT, GL_DEPTH_COMPONENT, _fboWidth, _fboHeight);
  glBindRenderbufferEXT(GL_RENDERBUFFER_EXT, 0);

  // Create framebuffer object
  glGenFramebuffersEXT(1, &_fboId);
  glBindFramebufferEXT(GL_FRAMEBUFFER_EXT, _fboId);

  // Assign COLOR, DEPTH and STENCIL buffers to the render buffer.
  glFramebufferTexture2DEXT(GL_FRAMEBUFFER_EXT, GL_COLOR_ATTACHMENT0_EXT, GL_TEXTURE_2D, _textureId, 0);
  glFramebufferRenderbufferEXT(GL_FRAMEBUFFER_EXT, GL_DEPTH_ATTACHMENT_EXT, GL_RENDERBUFFER_EXT, _rboId);
  //glFramebufferRenderbufferEXT(GL_FRAMEBUFFER_EXT, GL_STENCIL_ATTACHMENT_EXT, GL_RENDERBUFFER_EXT, _rboId);

  // Check if everything was ok.
  //bool _status = true; //checkFramebufferStatus();

  glBindFramebufferEXT(GL_FRAMEBUFFER_EXT, 0);
  
  return true;
}



void setupRC(){
	//glEnable(GL_COLOR_MATERIAL);
	//glColorMaterial(GL_FRONT, GL_AMBIENT_AND_DIFFUSE);
	glClearColor(0, 0, 0, 1.0);
	
	glDisable(GL_LINE_SMOOTH);
	glDisable(GL_LIGHTING);
	glDisable(GL_BLEND);
	initFBO();
}


void display(){
	glBindFramebufferEXT(GL_FRAMEBUFFER_EXT, _fboId);
	
	glPushAttrib(GL_VIEWPORT_BIT);
    glViewport(0, 0, _fboWidth, _fboHeight);

	    glMatrixMode(GL_MODELVIEW);
		glPushMatrix();
		glLoadIdentity();
		
		glDrawBuffer(GL_COLOR_ATTACHMENT0_EXT);
		glClearColor(0.0f, 0.0f, 0.0f, 0.00000f);
		glClear(GL_DEPTH_BUFFER_BIT  | GL_COLOR_BUFFER_BIT);
		
		glColor4f( 0.1f, 0.1f, 0.1f, 0.005f);
		glBegin(GL_QUADS);
			//glColor3f(0.f, 0.f, 0.5f);
			glVertex2f(0.f, 0.f);
			//glColor3f(0.5f, 0.f, 0.f);
			glVertex2f(1.f, 0.f);
			//glColor3f(0.f, 0.5f, 0.f); 
			glVertex2f(1.f, 1.f);
			glVertex2f(0.f, 1.f);
		glEnd();
		
		const int pixelNumber = _fboWidth * _fboHeight;
		float colorBuffer[pixelNumber];		  

		//GAGA
		glReadPixels(0, 0, _fboWidth, _fboHeight, GL_ALPHA, GL_FLOAT, colorBuffer);

		for (int i=0; i<pixelNumber; i++){
		  //printf("val: %hd\n", colorBuffer[i]);
		  std::cout << colorBuffer[i] << std::endl;
		}
		std::cout << "-MAX:"<< *(std::max_element(colorBuffer, colorBuffer + pixelNumber)) << "\n";
		
		
	glutSwapBuffers();
	glFinish();	
}



void reshape(int w, int h) {
  //	float aspect = w/h;	
  if(h==0) h = 1;
  glMatrixMode(GL_PROJECTION);
  glLoadIdentity();	
  glViewport(0, 0, w, h);
  gluOrtho2D(0, 1, 0, 1);
  glMatrixMode(GL_MODELVIEW);
}


void keyboard(unsigned char key, int x, int y){
	
}

int main (int argc, char **argv) {
	//std::cout << argv[1] << "\n";

	glutInit(&argc, argv);
	glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGBA | GLUT_DEPTH);
	glutInitWindowSize(512, 512);
	glutInitWindowPosition(300,300);
	glutCreateWindow("Dual Ray Implementation");
	glutDisplayFunc(display);	
	glutReshapeFunc(reshape);
	glutKeyboardFunc( keyboard );

	glewInit();
	  
	setupRC();
	

	//glutSpecialFunc( special );

	glutMainLoop();
	
	return 0;
}




