
#include <GL/glut.h>

#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <vector>

struct Point { float x,y; };
struct Triangle { Point a,b,c; };
struct Trapezoid { Point a,b,c,d; };

std::vector< Triangle > vecTri1;
std::vector< Triangle > vecTri2;
std::vector< Trapezoid > vecTrap;

using namespace std;

void displayQuad();
void displayTri1();
void displayTri2();

void readInput(){
	
	int t=0;
	bool terminate = true;
	Triangle tri;
	Trapezoid trap;
	
	while (terminate){
		cin >> t;
		cout << t << "\n";
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
			   terminate = false;
			   break;
		}		
	}
}

void setupRC(){
	glEnable(GL_COLOR_MATERIAL);
	glColorMaterial(GL_FRONT, GL_AMBIENT_AND_DIFFUSE);
	glClearColor(0.76, 0.56, 1.0, 1.0);
	
	glEnable(GL_BLEND);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
}

void reshape(int w, int h) {
  //	float aspect = w/h;	
  if(h==0) h = 1;
  glMatrixMode(GL_PROJECTION);
  glLoadIdentity();	
  glViewport(0, 0, w, h);
  //glOrtho(-1, 2, -1, 2, -1.0, 1.0);
  //glOrtho(-0.6, 1.6, -0.4, 2.2, -1.0, 1.0);
  gluOrtho2D(-0.1, 1.24, -2.1, 2.7);
  glMatrixMode(GL_MODELVIEW);
}

void display(){
	
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	
	
	glBegin(GL_TRIANGLES);
		for (unsigned int i=0; i < vecTri1.size(); ++i){
			glVertex2f(vecTri1[i].a.x,vecTri1[i].a.y);
			glVertex2f(vecTri1[i].b.x,vecTri1[i].b.y);
			glVertex2f(vecTri1[i].c.x,vecTri1[i].c.y);
		}
	glEnd();
	
	glBegin(GL_TRIANGLES);
		for (unsigned int i=0; i < vecTri2.size(); ++i){
			glVertex2f(vecTri2[i].a.x,vecTri2[i].a.y);
			glVertex2f(vecTri2[i].b.x,vecTri2[i].b.y);
			glVertex2f(vecTri2[i].c.x,vecTri2[i].c.y);
		}
	glEnd();
	
	glBegin(GL_QUADS);
		for (unsigned int i=0; i < vecTrap.size(); ++i){
			glVertex2f(vecTrap[i].a.x,vecTrap[i].a.y);
			glVertex2f(vecTrap[i].b.x,vecTrap[i].b.y);
			glVertex2f(vecTrap[i].c.x,vecTrap[i].c.y);
			glVertex2f(vecTrap[i].d.x,vecTrap[i].d.y);
		}
	glEnd();
	
	//displayTri2();
	//displayTri1();
	
	//displayQuad();
	
	glFlush();
	glutSwapBuffers();
	glFinish();
}


int main (int argc, char **argv) {

	glutInit(&argc, argv);
	glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGBA | GLUT_DEPTH);
	glutInitWindowSize(512, 512);
	glutInitWindowPosition(300,300);
	glutCreateWindow("Dual Ray Implementation");
	glutDisplayFunc(display);	

	readInput();
	setupRC();
	glutReshapeFunc(reshape);
	glutDisplayFunc(display);
	//glutKeyboardFunc( keyboard );
	//glutSpecialFunc( special );

	glutMainLoop();
	
	return 0;
}



void displayTri1(){
	glColor3f(0.0, 0.75, 0.2);
	

	glEnd();

}

void displayTri2(){
			glColor3f(0.5, 0.25, 0);
	
	glEnd();
}


void displayQuad(){
	glColor3f(0.5, 0, 0.7);
	

	glEnd();
}
