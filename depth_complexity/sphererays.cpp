//****************************************************************************
// File: sphere-rays.cpp
//
// Author: Joao Comba (comba@inf.ufrgs.br)
//
// Description: Draw lines connecting the projected spherical coordinates
//
//****************************************************************************
//
//      SPHERICAL:(R, THETA, PHI)
//      THETA [0,PI] 
//      PHI   [-PI,PI]
//
//      INPUT: THETA1, PHI1, THETA2, PHI2, DC
//      AXIS : X --> THETA, Y --> PHI
//****************************************************************************

#ifdef __APPLE__
#include <GLUT/glut.h>
#else
#include <GL/glut.h>
#endif

#ifdef WIN32
#  include <windows.h>
#endif

#include<iostream>
#include<vector>
#include<math.h>

//#include "vector.hpp"
#include "camera/Camera.h"
#include "Eigen/Dense"

using Eigen::Affine3f;
using Eigen::Matrix3f;
using Eigen::Matrix4f;
using Eigen::Vector3f;
using Eigen::Vector4f;
using Eigen::AngleAxisf;

using namespace std;

#define PI 3.1415926

typedef struct {
  float theta1, phi1, theta2, phi2;
  int dc;
} DCLine;

typedef struct  {
    int x,y;
    int Dx, Dy;
    bool buttons[3];
    float sensitivity;
} Mouse;

struct SegmentCart{
    Vector3f a;
    Vector3f b;
};

struct SegmentSph{
    Vector3f sph_a;
    Vector3f sph_b;
};

Camera camera;
Mouse mouse;

Vector3f rotationAxis;

int thetaDom1 = 0; 
int thetaDom2 = PI;
int phiDom1   = -PI;
int phiDom2   = PI;


int fps = 30;

int number_of_lines, maxDC;
vector<DCLine> DClines;
SegmentSph mainSegment[3];
double s_x[4], s_y[4];
int win_id[2] = {0,0};
bool updateAngle = true;



//-----------------------------------------------------------------------------
// conversion
//-----------------------------------------------------------------------------
Vector3f cartesianToSpherical(Vector3f point) {
  Vector3f ans;
  double a = point(0), b = point(1), c = point(2);
  ans(0) = sqrt(a*a + b*b + c*c);
  ans(1) = acos(c/ans(0));
  ans(2) = atan2(b,a);
  return ans;
}

Vector3f sphericalToCartesian(Vector3f point) {
  Vector3f ans;
  double r = point(0), theta = point(1), phi = point(2);
  ans(0) = r*sin(theta)*cos(phi); 
  ans(1) = r*sin(theta)*sin(phi);
  ans(2) = r*cos(theta);
  return ans;
}



//-----------------------------------------------------------------------------
// displayBox
//-----------------------------------------------------------------------------
void
displayBox() {
  glColor3f(0, 0, 0);
  glLineWidth(2.0);
  glPolygonMode(GL_FRONT, GL_LINE);
  glBegin(GL_POLYGON);
  glVertex3f(s_x[0], s_y[0], 0.0);
  glVertex3f(s_x[1], s_y[1], 0.0);
  glVertex3f(s_x[3], s_y[3], 0.0);
  glVertex3f(s_x[2], s_y[2], 0.0);
  glEnd();
  glLineWidth(1.0);
}

float normalizey(float angle) {
  return ((angle/3.1415926+1.0f)/2.0f);
}

float normalizex(float angle) {
  return (angle/3.1415926);
}

//-----------------------------------------------------------------------------
// display
//-----------------------------------------------------------------------------
void 
display(void) 
{
  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
  glEnable (GL_BLEND);
  glBlendFunc (GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
  
    glLineWidth(1.0);
    glColor4f(0,0,0,1);
    glBegin(GL_LINES);
    for (unsigned i=0; i<3; i++){
        switch(i) {
        case 0: glColor4f(1, 0, 0, 1); break;
        case 1: glColor4f(0, 1, 0, 1); break;
        case 2: glColor4f(0, 0, 1, 1); break;
        }
        glVertex3f(mainSegment[i].sph_a(1), mainSegment[i].sph_a(2), 0.0);
        glVertex3f(mainSegment[i].sph_b(1), mainSegment[i].sph_b(2), 0.0);
    }
    glEnd();

  glFlush();
  glutSwapBuffers();
  glFinish();
  
  glutPostRedisplay();
}

float angle=0;

void idle(){
    
}

#include<unistd.h>

void
display2(){
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();
    
    sleep(0.25);
    
    camera.update();    
    camera.lookAt();

    // Draw Wire Sphere
    glLineWidth(1.0f);
    glColor4f(0,0,0,0.2);
    glutWireSphere(2,12,12);
    
    glLineWidth(1.0f);
    glColor4f(0.0, 1.0, 0.3, 1.0);
    SegmentCart d[3];
    
    for (unsigned i=0; i<3; i++){
        if (updateAngle){        

            //rotT.rotate();
            //d[i].a = sphericalToCartesian(mainSegment[i].sph_a);
            //d[i].b = sphericalToCartesian(mainSegment[i].sph_b);

            d[i].a = AngleAxisf(angle, rotationAxis)*sphericalToCartesian(mainSegment[i].sph_a);
            d[i].b = AngleAxisf(angle, rotationAxis)*sphericalToCartesian(mainSegment[i].sph_b);

            //std::cout << d[i].a << "\n" << d[i].b << "\n";
            mainSegment[i].sph_a = cartesianToSpherical(d[i].a);
            mainSegment[i].sph_b = cartesianToSpherical(d[i].b);


            
        }else {
            d[i].a = sphericalToCartesian(mainSegment[i].sph_a);
            d[i].b = sphericalToCartesian(mainSegment[i].sph_b);
        }
    }
    updateAngle = false;
    
    glLineWidth(3.0f);
    glBegin(GL_LINES);
    for (unsigned i=0; i<3; i++){
        switch(i) {
            case 0: glColor4f(1, 0, 0, 1); break;
            case 1: glColor4f(0, 1, 0, 1); break;
            case 2: glColor4f(0, 0, 1, 1); break;
        }
        glVertex3fv(d[i].a.data());
        glVertex3fv(d[i].b.data());
    }
    glEnd();
    
    glPointSize(8.0f);
    glBegin(GL_POINTS);
    for (unsigned i=0; i<3; i++){
        switch(i) {
            case 0: glColor4f(1, 0, 0, 1); break;
            case 1: glColor4f(0, 1, 0, 1); break;
            case 2: glColor4f(0, 0, 1, 1); break;
        }
        glVertex3f(d[i].a(0), d[i].a(1), d[i].a(2));
        glVertex3f(d[i].b(0), d[i].b(1), d[i].b(2));
    }
    glEnd();
    
    Vector3f rot;
    rot << 2*rotationAxis;
    glLineWidth(6.0f);
    glColor4f(0.0, 0.0, 0.0, 1.0);
    glBegin(GL_LINES);
    glVertex3f(0,0,0); glVertex3fv(rot.data());
    glEnd();
    
    // Draw Axis
    glLineWidth(1.0f);
    glColor4f(0.0, 0.0, 0.0, 1.0);
    glBegin(GL_LINES);
    glVertex3f(0,0,0); glVertex3f(100,0,0);
    glVertex3f(0,0,0); glVertex3f(0,100,0);
    glVertex3f(0,0,0); glVertex3f(0,0,100);
    glEnd();
    
    

    
    glFlush();
    glutSwapBuffers();
    glFinish();  
    glutPostRedisplay();
}

//-----------------------------------------------------------------------------
// reshape
//-----------------------------------------------------------------------------
void 
reshape(int w, int h) {
  if(h==0) h = 1;
  glMatrixMode(GL_PROJECTION);
  glLoadIdentity();	
  glViewport(0, 0, w, h);
  glOrtho(thetaDom1, thetaDom2, phiDom1, phiDom2, -1.0, 1.0);
  glMatrixMode(GL_MODELVIEW);

}

void 
reshape2(int w, int h) {
  //	float aspect = w/h;	
  if(h==0) h = 1;
  glMatrixMode(GL_PROJECTION);
  glLoadIdentity();	
  glViewport(0, 0, w, h);
  camera.setPerspec(60, (double)w/h, 0.1, 100);

  glMatrixMode(GL_MODELVIEW);
  camera.update();	
  camera.lookAt();
}

//-----------------------------------------------------------------------------
// setupRC
//-----------------------------------------------------------------------------
void 
setupRC(void) {
    glEnable(GL_COLOR_MATERIAL);
    glColorMaterial(GL_FRONT, GL_AMBIENT_AND_DIFFUSE);
    glClearColor(1.0, 1.0, 1.0, 1.0);
  
    
}

float lightAmbient[4] = {0.2, 0.2, 0.2, 1.0};
float lightDiffuse[4] = {0.2, 0.2, 0.2, 1.0};
float lightDir[4] = {0.2, -0.6, 0.2, 0.0};

void 
setupRC2(void) {
    glClearColor(1.0, 1.0, 1.0, 1.0);
    
    glEnable(GL_COLOR_MATERIAL);
    glColorMaterial(GL_FRONT, GL_AMBIENT_AND_DIFFUSE);   
    
    glEnable (GL_BLEND);
    glBlendFunc (GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
  
    // mouse
    mouse.sensitivity = 0.0025;
    mouse.x = mouse.y = 0;
    
    rotationAxis << 1,3,-0.82;
    rotationAxis.normalize();
    std::cout << "rotationAxis:\n"<< rotationAxis << std::endl;
}

//-----------------------------------------------------------------------------
// keyboard
//-----------------------------------------------------------------------------
void keyboard( unsigned char key, int x, int y ) {
}

//-----------------------------------------------------------------------------
// special
//-----------------------------------------------------------------------------
void special1( int key, int x, int y )
{
    glutPostRedisplay();
}

void special2( int key, int x, int y )
{
    for (unsigned i=0; i<3; i++){
        if (key == GLUT_KEY_UP){   
            mainSegment[i].sph_a(1) += 0.05;
            mainSegment[i].sph_b(1) -= 0.05;
        }
        if (key == GLUT_KEY_DOWN){
            mainSegment[i].sph_a(1) -= 0.05;
            mainSegment[i].sph_b(1) += 0.05;
        }
        if (key == GLUT_KEY_LEFT){
            mainSegment[i].sph_a(2) += 0.05;
            mainSegment[i].sph_b(2) += 0.05;
        }
        if (key == GLUT_KEY_RIGHT){
            mainSegment[i].sph_a(2) -= 0.05;
            mainSegment[i].sph_b(2) -= 0.05;
        }

        if (mainSegment[i].sph_a(1) > PI)  mainSegment[i].sph_a(1) -= PI;
        if (mainSegment[i].sph_a(1) < 0 )  mainSegment[i].sph_a(1) += PI;
        if (mainSegment[i].sph_a(2) > PI)  mainSegment[i].sph_a(2) -= 2*PI;
        if (mainSegment[i].sph_a(2) <= -PI) mainSegment[i].sph_a(2) += 2*PI;

        if (mainSegment[i].sph_b(1) > PI)  mainSegment[i].sph_b(1) -= PI;
        if (mainSegment[i].sph_b(1) < 0 )  mainSegment[i].sph_b(1) += PI;
        if (mainSegment[i].sph_b(2) > PI)  mainSegment[i].sph_b(2) -= 2*PI;
        if (mainSegment[i].sph_b(2) <= -PI) mainSegment[i].sph_b(2) += 2*PI;
    }
    
  glutPostRedisplay();
}

//-----------------------------------------------------------------------------
// mouse motion
//-----------------------------------------------------------------------------
void mouseMotion2(int x, int y){
    
    mouse.Dx = mouse.x - x;
    mouse.Dy = mouse.y - y;
    
    // rotate camera 
    if (mouse.buttons[0]){        
        camera.lookLefRigObj(mouse.Dx*mouse.sensitivity);
        camera.lookUpDownObj(mouse.Dy*mouse.sensitivity);
    }
    // zoom camera
    else if (mouse.buttons[2]){
        camera.MoveFrente(mouse.Dy*mouse.sensitivity);
    }
    
    // update mouse coords
    mouse.x = x;
    mouse.y = y;
    
    glutPostRedisplay();
}


void mouseFunc2(int button, int state, int x, int y){    
    // update mouse buttons and coordinates
    mouse.buttons[button] = !mouse.buttons[button];
    mouse.x = x;
    mouse.y = y;
}

//-----------------------------------------------------------------------------
// timer motion
//-----------------------------------------------------------------------------
void update(int v){
    angle = 0.02;
    glutTimerFunc(1000/fps, update, 0);
    updateAngle = true;
}


//-----------------------------------------------------------------------------
// main
//-----------------------------------------------------------------------------
int 
main (int argc, char **argv)
{		
  glutInit(&argc, argv);
  glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGBA | GLUT_DEPTH);
  // ************ WINDOW 1
  glutInitWindowSize(900, 900);
  win_id[0] = glutCreateWindow("Dual Ray Implementation");
  glutInitWindowPosition (10, 100);
  setupRC();
  glutDisplayFunc(display);
  glutReshapeFunc(reshape);
  glutKeyboardFunc( keyboard );
  glutSpecialFunc( special1 );
  
  
  // ************ WINDOW 2
  win_id[1] = glutCreateWindow("Spherical Coordiantes");
  glutInitWindowPosition (500, 100);
  setupRC2();  
  glutDisplayFunc(display2);
  glutReshapeFunc(reshape2);
  glutMouseFunc(mouseFunc2);
  glutMotionFunc(mouseMotion2);
  //glutKeyboardFunc( keyboard );
  glutSpecialFunc( special2 );
  glutTimerFunc(1000/fps,update,0);
  glutIdleFunc(idle);
  
  
  //* -----------------------------------------------
  s_x[0] = 0.0; s_y[0] = -PI;  
  s_x[1] = PI; s_y[1] = -PI;  
  s_x[2] = PI; s_y[2] = PI;  
  s_x[3] = 0.0; s_y[3] = -PI;	
  
  Vector3f p1,p2;
  
  p1 << 2,PI/2.0f, -PI;
  p2 << 2,PI/2.0f, 0;

  mainSegment[0].sph_a = p1;
  mainSegment[0].sph_b = p2;
  
  p1 << 2, 0, -PI/2.0f;
  p2 << 2, PI, 0;
  
  mainSegment[1].sph_a = p1;
  mainSegment[1].sph_b = p2;
  
  p1 << 2, PI/2, PI/2;
  p2 << 2, +PI/2, -PI/2;
  
  mainSegment[2].sph_a = p1;
  mainSegment[2].sph_b = p2;
  
  glutMainLoop();

  return 0;
}
