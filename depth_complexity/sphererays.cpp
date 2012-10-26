//****************************************************************************
// File: sphere-rays.cpp
//
// Author: Joao Comba (comba@inf.ufrgs.br)
//
// Description: Draw lines connecting the projected spherical coordinates
//
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
using namespace std;

typedef struct {
  float alpha1, omega1, alpha2, omega2;
  int dc;
} DCLine;

int number_of_lines, maxDC;
vector<DCLine> DClines;
double s_x[4], s_y[4];

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
  displayBox();
  float x1, y1, x2, y2; 

  glEnable (GL_BLEND);
  glBlendFunc (GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

  for (int i=number_of_lines-1; i >= 0; i--) {
    switch(maxDC-DClines[i].dc) {
    case 0: glColor4f(1, 0, 0, 1); break;
    case 1: glColor4f(0, 1, 0, 0.5); break;
    case 2: glColor4f(0, 0, 1, 0.1); break;
    case 3: glColor4f(1, 0, 1, 0.1); break;
    case 4: glColor4f(0, 1, 1, 0.1); break;
    case 5: glColor4f(1, 1, 0, 0.1); break;
    }
    // Perform the proper scaling - alpha 1 ranges from -PI to PI
    //    cout << "\n" << DClines[i].alpha1 << "," << DClines[i].omega1 << " -  " << DClines[i].alpha2 << "," << DClines[i].omega2;
    x1 = normalizex(DClines[i].alpha1);
    y1 = normalizey(DClines[i].omega1);
    x2 = normalizex(DClines[i].alpha2);
    y2 = normalizey(DClines[i].omega2);
    //cout << "\n" << x1 << "," << y1 << " -  " << x2 << "," << y2;
    glLineWidth(1.0);
    glBegin(GL_LINES);
    glVertex3f(x1, y1, 0.0);
    glVertex3f(x2, y2, 0.0);
    glEnd();
  }
  glFlush();
  glutSwapBuffers();
  glFinish();
}

//-----------------------------------------------------------------------------
// reshape
//-----------------------------------------------------------------------------
void 
reshape(int w, int h) {
  //	float aspect = w/h;	
  if(h==0) h = 1;
  glMatrixMode(GL_PROJECTION);
  glLoadIdentity();	
  glViewport(0, 0, w, h);
  glOrtho(0, 1, 0, 1, -1.0, 1.0);
  glMatrixMode(GL_MODELVIEW);

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

//-----------------------------------------------------------------------------
// keyboard
//-----------------------------------------------------------------------------
void keyboard( unsigned char key, int x, int y ) {
}

//-----------------------------------------------------------------------------
// special
//-----------------------------------------------------------------------------
void special( int key, int x, int y )
{
  glutPostRedisplay();
}

//-----------------------------------------------------------------------------
// main
//-----------------------------------------------------------------------------
int 
main (int argc, char **argv)
{		
  glutInit(&argc, argv);
  glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGBA | GLUT_DEPTH);
  glutInitWindowSize(1200, 1200);
  glutCreateWindow("Dual Ray Implementation");
  glutDisplayFunc(display);

  s_x[0] = 0.0; s_y[0] = 0.0;
  s_x[1] = 1.0; s_y[1] = 0.0;
  s_x[2] = 0.0; s_y[2] = 1.0;
  s_x[3] = 1.0; s_y[3] = 1.0;	

  setupRC();
  glutReshapeFunc(reshape);
  glutKeyboardFunc( keyboard );
  glutSpecialFunc( special );

  int first = 1;

  cin >> number_of_lines ;

  DClines.clear(); DClines.resize(number_of_lines);

  //float alpha1, omega1, alpha2, omega2;

  for (int i=0; i < number_of_lines; i++) {
    DCLine line;
    cin >> line.alpha1 >> line.omega1 >> line.alpha2 >> line.omega2 >> line.dc;
    if (first) {
      first = 0;
      maxDC = line.dc;
    }
    //cout << "\n" << line.alpha1 << "," << line.omega1 << " -  " << line.alpha2 << "," << line.omega2;
    DClines.push_back(line);
    DClines[i].alpha1 = line.alpha1;
    DClines[i].alpha2 = line.alpha2;
    DClines[i].omega1 = line.omega1;
    DClines[i].omega2 = line.omega2;
    DClines[i].dc = line.dc;
    //cout << DClines[i].alpha1 << "," << DClines[i].omega1 << " -  " << DClines[i].alpha2 << "," << DClines[i].omega2;
  }

  glutMainLoop();

  return 0;
}
