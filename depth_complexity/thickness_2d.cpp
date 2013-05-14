#include <iostream>

#include <GL/glut.h>
#include <vector>
#include <math.h>
//#include <GL/glew.h>

#include "vector.hpp"

#define NUMBER_OF_RAYS 100000
#define OBJECT_TYPE OBJ_TRAPEZOID
#define OUTPUT_FORMAT TWO_COLUMNS_WITH_LABEL

enum OutputFormat{ ONLY_ONE_COLUMN, TWO_COLUMNS, TWO_COLUMNS_WITH_LABEL };
enum Object{ OBJ_SQUARE, OBJ_TRIANGLE, OBJ_CIRCLE, OBJ_RECTANGLE, OBJ_TRAPEZOID};

const double EPS = 1.0E-5;
#define perp(u,v)  ((u).x * (v).y - (u).y * (v).x)  // perp product (2D)


//STRUCTS ---------------------------------------------------------------------
struct Segment{
    vec3d a,b;    
};

// geometry
struct Mesh {
    std::vector<vec3d> vlist;   
};

//struct Cube {
//    vec3d v[4];
//};

//GLOBALS
float histMin, histMax;
int intervals;
float boundingRadius;
Mesh obj;

std::vector<Segment> lines;
std::vector<vec3d> intersection;
std::vector<double> thickness;
std::vector<long long unsigned int> histogram;



// AUXILIAR COLLISION

double uniformRandom() { return rand()/double(RAND_MAX); }

// lineIntersect in 2D
bool lineIntersection2D(const Segment &line1, const Segment &line2, double *t1, double *t2) {
  vec3d u = line1.b - line1.a;
  vec3d v = line2.b - line2.a;
  vec3d w = line1.a - line2.a;
  // just X and Y values are used
  double D = perp(u,v);

  // test if they are parallel (includes either being a point)
  if (fabs(D) < EPS) // S1 and S2 are parallel
    return false;

  // the segments are skew and may intersect in a point
  // get the intersect parameter for S1
  *t1 = perp(v,w) / D;
  *t2 = perp(u,w) / D;

  return true;
}

// segmentIntersection2D(): the intersection of 2 finite 2D segments
bool segmentIntersection2D(const Segment &seg1, const Segment &seg2, double *t1, double *t2){
  if (lineIntersection2D(seg1, seg2, t1, t2) and
      (0.0 < *t1 && *t1 < 1.0) and
      (0.0 < *t2 && *t2 < 1.0))
    return true;
  return false;
}


void createRays(){
    
    //std::cout << uniformRandom() << std::endl;
    for (int i=0; i<NUMBER_OF_RAYS; i++){
        Segment s1;
        float val1 = 2*M_PI*uniformRandom();
        float val2 = 2*M_PI*uniformRandom();
        s1.a = vec3d(boundingRadius*sin(val1), boundingRadius*cos(val1), 0);
        s1.b = vec3d(boundingRadius*sin(val2), boundingRadius*cos(val2), 0);
        lines.push_back(s1);
        
        std::vector<vec3d> distanceThick;
        
        // check intersection points
        for (unsigned int j=0; j<obj.vlist.size(); j++){
            double t1,t2; //intersection info
            Segment s2;
            s2.a = obj.vlist[j];
            s2.b = obj.vlist[(j+1)%obj.vlist.size()];
            if (segmentIntersection2D(s1,s2,&t1,&t2)){
                vec3d ipoint = s1.a*(1.0f-t1) + t1*s1.b;
                intersection.push_back(ipoint); 
                distanceThick.push_back(ipoint);
            }
        }
        if (distanceThick.size()==2){
            double thicknessValue = sqrt(pow(distanceThick[1].x - distanceThick[0].x, 2) + 
                                    pow(distanceThick[1].y - distanceThick[0].y, 2));
            
            thickness.push_back(thicknessValue);
        }
    }
    
    //for (unsigned int i=0; i<thickness.size(); i++){
    //    std::cout << thickness[i] << std::endl;
        
    //}
    
    for (unsigned int i=0; i<thickness.size(); i++){
        double factor = static_cast<double>(histMax)/intervals;
        histogram[ thickness[i] / factor ] += 1;   
    }    
    
    for (unsigned int i=0; i<histogram.size(); i++){
       
        switch(OUTPUT_FORMAT){
            case ONLY_ONE_COLUMN:
                std::cout << histogram[i] << std::endl;
            break;
            
            case TWO_COLUMNS:
                std::cout << i << "\t" << histogram[i] << std::endl;
            break;
                
            case TWO_COLUMNS_WITH_LABEL:
                if (i==0)
                    std::cout << "Thickness\tFrequency" << std::endl;
                
                std::cout << i << "\t" << histogram[i] << std::endl;
            break;
        }
        
    }
    
    
}

void initialize(){
    glClearColor(1.0f, 1.0f, 1.0f, 1.0f);
    
    switch(OBJECT_TYPE){
        case OBJ_SQUARE:
            obj.vlist.push_back(vec3d(-1,-1,0));
            obj.vlist.push_back(vec3d(-1,1,0));
            obj.vlist.push_back(vec3d(1,1,0));
            obj.vlist.push_back(vec3d(1,-1,0));
        break;
        
        case OBJ_TRIANGLE:
            obj.vlist.push_back(vec3d(-1,-1,0));
            obj.vlist.push_back(vec3d(0,2*0.73,0));
            obj.vlist.push_back(vec3d(1,-1,0));
        break;
        
        case OBJ_CIRCLE:
            for (float i=0; i<2*M_PI; i+=2*M_PI/16.0 )
                obj.vlist.push_back(vec3d(sqrt(2)*sin(i),sqrt(2)*cos(i),0));
        break;
        
        case OBJ_RECTANGLE:
            obj.vlist.push_back(vec3d(-0.5,-1,0));
            obj.vlist.push_back(vec3d(-0.5,1,0));
            obj.vlist.push_back(vec3d(0.5,1,0));
            obj.vlist.push_back(vec3d(0.5,-1,0));
        break;
        
        case OBJ_TRAPEZOID:
            obj.vlist.push_back(vec3d(-1,-1,0));
            obj.vlist.push_back(vec3d(-1,0.25,0));
            obj.vlist.push_back(vec3d(1,1,0));
            obj.vlist.push_back(vec3d(1,-1,0));
        break;
    }
    
    
    
    intervals = 100;
    histogram.resize(intervals,0);
    histMin = 0;
    histMax = 2*sqrt(2)*1.1;
    boundingRadius = sqrt(2)*1.1;

    
    createRays();
}



void display(){
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    
        //draw axis   
    glBegin(GL_LINES);
        glColor3f(1.0,0,0);
        glVertex2i(-10,0);
        glVertex2i(10,0);
        glColor3f(0.0,1.0,0);
        glVertex2i(0,-10);
        glVertex2i(0,10);        
    glEnd();
    
    
    glColor3f(0,0,0);
    
    
    for (unsigned int i=0; i<obj.vlist.size(); i++){
        glBegin(GL_LINE_STRIP);
            glVertex2f(obj.vlist[i].x, obj.vlist[i].y);
            glVertex2f(obj.vlist[(i+1)%obj.vlist.size()].x, obj.vlist[(i+1)%obj.vlist.size()].y);

            //glVertex2f(tri.v[1].x, tri.v[1].y);
            //glVertex2f(tri.v[2].x, tri.v[2].y);

            //glVertex2f(tri.v[2].x, tri.v[2].y);
            //glVertex2f(tri.v[0].x, tri.v[0].y);
        glEnd();
    }
    
    
    /*
    for (unsigned int i=0; i<lines.size(); i++){        
         glBegin(GL_LINES);
         glVertex2f(lines[i].a.x, lines[i].a.y);  
         glVertex2f(lines[i].b.x, lines[i].b.y);     
         
         glEnd();
    }
     */
    
    
    glLineWidth(2.0);
    //draw bounding circle
    glBegin(GL_LINE_STRIP);
    for (float i=0; i<2*M_PI; i+=2*M_PI/30.0){
        glVertex2f(boundingRadius*sin(i),boundingRadius*cos(i));
    }
    glEnd();
    
    /*
    glColor3f(1,0,0);
    glPointSize(5.0f);
    glBegin(GL_POINTS);
    for (unsigned int i=0; i<intersection.size(); i++){
        glVertex2f(intersection[i].x,intersection[i].y);
    }
    glEnd();
    */

    
    glutSwapBuffers();
}

void reshape(int h, int w){
     // Prevent a divide by zero, when window is too short
    if(h == 0) h = 1; 
    //mainWin.x = w;
    //mainWin.y = h;
    // Use the Projection Matrix
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();

    // Set the viewport to be the entire window
    glViewport(0, 0, 600, 600);

    //gluPerspective(45, (double)w/h, 0.5, 1000);
    gluOrtho2D(-3,3,-3,3);
    //camera.setPerspective(45, (double)w/h, 0.1, 250);

    // Get Back to the Modelview
    glMatrixMode(GL_MODELVIEW);     
}

int main(int argc, char** argv){
    /* initialize random seed: */
    srand (time(NULL));
    
    glutInit(&argc, argv);
    glutInitDisplayMode(GLUT_RGB | GLUT_DOUBLE | GLUT_DEPTH );
    glutInitWindowPosition(300, 300);
    glutInitWindowSize(600, 600);
    glutCreateWindow("Thickness 2D");    
    
    // callbacks
    glutDisplayFunc(display);
    glutReshapeFunc(reshape);
    
    initialize();
    glutMainLoop();
    return 0;   
}
