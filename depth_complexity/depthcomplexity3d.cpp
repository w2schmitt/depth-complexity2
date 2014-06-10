#include <iostream>
#include <algorithm>
#include <istream>
#include <cmath>
#include <string>
#include <sstream>
#include <cstdlib>
#include <vector>
#include <assert.h>
#include <fstream>
#include <map>
#include <list>
#include <set>



#include <GL/glew.h>
#include <GLFW/glfw3.h>
#include <AntTweakBar.h>
#ifdef __APPLE__
#include <GLUT/glut.h>
#else
#include <GL/glut.h>
#endif

#include "flags.h"
#include "vector.hpp"
#include "camera/float3.h"
#include "camera/Camera.h"
#include "util.h"
#include "RFDepthComplexity3D.h" 
#include "timer.h"
#include "CImg.h"

using namespace cimg_library;

float radius = 0.01;
bool showPlanes = false;
bool doGoodRays = true;
bool showMaxRays = true;
bool showGoodRays = false;
bool thresholdRays = false;
bool highlightTris = false;
bool discretePts = true;
bool showScale = false;
vec4f sphereColor(0.0, 1.0, 0.0, 1.0);
unsigned int showRayIndex = 0;
unsigned int showRayThicknessIndex = 0;
GLFWwindow* window;

// for 3d textue display
CImgDisplay main_disp;
CImg<float> cscale;

// DEPTH COMPLEXITY 3D ---
RFDepthComplexity3D *dc3d; // = new RFDepthComplexity3D(512, 50);

//#define COMPUTE_OFFLINE

//#ifdef USE_RANDOM_DC3D
//RFDepthComplexity3D *dc3d;
//#else
//DepthComplexity3D *dc3d;
//#endif

template<class T>
T clamp(const T& min, const T& x, const T& max)
{
    return x < min ? min : (x > max ? max : x);
}

void drawPlane(const vec4d& plane)
{
    const vec3d& n = plane.xyz();
    
    // axis = axis of min rotation from z to n
    const vec3d& axis = cross(vec3d(0,0,1), n);
    
    double degrees = std::asin(clamp(-1.0, axis.length(), 1.0)) * 180 / M_PI;
    
    glPushMatrix();
      glRotatef(degrees, axis.x, axis.y, axis.z);
      glTranslatef(0, 0, -plane.w);
    
      glEnable(GL_BLEND);
      glBlendFunc(GL_ONE, GL_ONE_MINUS_SRC_ALPHA);

      glBegin(GL_QUADS);
        glColor4f(0, .2, .5, 0.3);
        glTexCoord2f(0.0, 0.0); glVertex2f(-2, -2);
        glTexCoord2f(1.0, 0.0); glVertex2f(2, -2);
        glTexCoord2f(1.0, 1.0); glVertex2f(2, 2);
        glTexCoord2f(0.0, 1.0); glVertex2f(-2, 2);
      glEnd();
    glPopMatrix();
    glDisable(GL_BLEND);
}


void drawQuadTex(const Point &p1, const Point &p2, const Point &p3, const Point &p4) {
  glColor4f(1.f, 0.f, 0.f, 1.0);
  glBegin(GL_QUADS);
    glTexCoord2f(0.0, 0.0); glVertex2f(p1.x, p1.y);
    glTexCoord2f(1.0, 0.0); glVertex2f(p2.x, p2.y);
    glTexCoord2f(1.0, 1.0); glVertex2f(p3.x, p3.y);
    glTexCoord2f(0.0, 1.0); glVertex2f(p4.x, p4.y);
  glEnd();
}

static int winWidth, winHeight;


// Callback function called by GLFW when window size changes
void WindowSizeCB(GLFWwindow* w, int width, int height)
{
    //std::cout << width << " " << height << std::endl;
    winWidth = width;
    winHeight = height;
    // Send the new window size to AntTweakBar
    TwWindowSize(width, height);
}

vec3d triCenter(const Triangle& t)
{
    return (t.a + t.b + t.c) / 3.0;
}

template<class T>
vec3<T> sorted(const T& a, const T& b, const T& c)
{
    vec3<T> r(a, b, c);
    if (r.y < r.x)
        std::swap(r.x, r.y);
    if (r.z < r.x)
        std::swap(r.x, r.z);
    if (r.z < r.y)
        std::swap(r.y, r.z);
    return r;
}

struct ByDist {
    const vec3d dir;
    ByDist(const vec3f& d) : dir(d) {}
    
    bool operator()(const Triangle& t1, const Triangle& t2) const
    {
#if 1
        vec3d c1 = triCenter(t1);
        vec3d c2 = triCenter(t2);
        return dot(c1, dir) > dot(c2, dir);
#else
        double da1 = dot(dir, t1.a);
        double db1 = dot(dir, t1.b);
        double dc1 = dot(dir, t1.c);

        double da2 = dot(dir, t2.a);
        double db2 = dot(dir, t2.b);
        double dc2 = dot(dir, t2.c);

        vec3d s1 = sorted(da1, db1, dc1);
        vec3d s2 = sorted(da2, db2, dc2);

        if (s1.x > s2.x)
            return true;
        if (s1.x < s2.x)
            return false;
        if (s1.y > s2.y)
            return true;
        if (s1.y < s2.y)
            return false;
        return s1.z > s2.z;
#endif
    }
};

std::vector<Triangle> sorted_faces;

void drawMesh(const TriMesh& mesh, const vec3f& dir)
{
 
    glEnable(GL_TEXTURE_3D);
    glBindTexture(GL_TEXTURE_3D, dc3d->getTextureID());
    glColorMaterial(GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE);
    glEnable(GL_COLOR_MATERIAL);
    
    // set texture 3d shader
    dc3d->setShaderTex3d();
    
    //std::clog << "sorting...";
    if (sorted_faces.empty()) {
        sorted_faces = mesh.faces;
    }
    
    std::sort(sorted_faces.begin(), sorted_faces.end(), ByDist(dir));
    //std::clog << "done" << std::endl;
    
    glEnable(GL_BLEND);
    glBlendFuncSeparate(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA, GL_ONE, GL_ZERO);
    //glBlendFunc(GL_ONE, GL_ONE_MINUS_SRC_ALPHA);
    //glDisable(GL_DEPTH_TEST);
    glEnable(GL_DEPTH_TEST);
    glEnable(GL_LIGHTING);
    glLightModeli(GL_LIGHT_MODEL_TWO_SIDE, 1);    

    glEnable(GL_VERTEX_ARRAY);
    //glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);

    glEnableClientState(GL_VERTEX_ARRAY);    
    glVertexPointer(3, GL_DOUBLE, 3*sizeof(vec3d)+sizeof(vec4d), &sorted_faces[0].a.x);
    
    glEnableClientState (GL_TEXTURE_COORD_ARRAY);
    glTexCoordPointer(3, GL_DOUBLE, 3*sizeof(vec3d)+sizeof(vec4d), &sorted_faces[0].tca.x);
    
    glEnableClientState(GL_COLOR_ARRAY);
    glColorPointer(4, GL_DOUBLE, 3*sizeof(vec3d)+sizeof(vec4d), &sorted_faces[0].ca.x);

    glEnableClientState(GL_NORMAL_ARRAY);
    glNormalPointer(GL_DOUBLE, 3*sizeof(vec3d)+sizeof(vec4d), &sorted_faces[0].na.x);

    glDrawArrays(GL_TRIANGLES, 0, sorted_faces.size()*3);
   
    
    glUseProgram(0);

    //glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
    
    glDisableClientState(GL_VERTEX_ARRAY);
    glDisable(GL_VERTEX_ARRAY);
    glDisableClientState(GL_TEXTURE_COORD_ARRAY);
    glDisable(GL_TEXTURE_COORD_ARRAY);
    glDisableClientState(GL_NORMAL_ARRAY);
    glDisable(GL_NORMAL_ARRAY);
    glDisableClientState(GL_COLOR_ARRAY);
    glDisable(GL_COLOR_ARRAY);

    glDisable(GL_BLEND);
    glEnable(GL_DEPTH_TEST);
    glDisable(GL_COLOR_MATERIAL);
    
    glDisable(GL_TEXTURE_3D);
    glBindTexture(GL_TEXTURE_3D, 0);    
}

//void drawPlaneIntersection(const std::vector<Segment>& segments)
//{
//    // draw the plane-mesh intersection
//    glDisable(GL_LIGHTING);
//    glLineWidth(2);
//    glBegin(GL_LINES);
//    glColor3f(1, 0, 0);
//    for (unsigned i=0; i<intersectionVectors.size(); ++i) {
//        const Segment & s = intersectionVectors[i];
//        glVertex3f(s.a.x, s.a.y, s.a.z);
//        glVertex3f(s.b.x, s.b.y, s.b.z);
//    }
//    glEnd();
//}

unsigned rayIndexT = 0;
    
void drawRays()
{
    /*
    const std::vector<Segment> &s = dc3d->cameraRays();
    
    if (s.size() > 0){
        glLineWidth(2.5);
        glBegin(GL_LINES);
        glColor3f(0.5, 0.0, 0.5);
        for (unsigned int i=0; i<s.size(); i++){
            glVertex3f(s[i].a.x, s[i].a.y, s[i].a.z);
            glVertex3f(s[i].b.x, s[i].b.y, s[i].b.z);
        }
        glEnd();
    }
    /**/ 


    
    if (showMaxRays){
        const std::set<Segment,classcomp>& rays = dc3d->maximumRays();
        std::set<Segment,classcomp>::const_iterator it = rays.begin();
        if (rays.size()>0){
            // draw rays
            glLineWidth(2.5);
            glBegin(GL_LINES);
            glColor3f(0.5, 0.0, 0.5);
              for (; it!=rays.end(); ++it) {
                const Segment &r = *it;
                
                //glVertex3f(0,0,0);
                //glVertex3f(50,50,50);

                glVertex3f(100*r.a.x, 100*r.a.y, 100*r.a.z);
                glVertex3f(100*r.b.x, 100*r.b.y, 100*r.b.z);
              }
            glEnd();
        }
    }
        
        // HIGHLIGHT TRIANGLES 
        /*
        if (highlightTris){
            const std::vector<Triangle> &interTris = dc3d->intersectionTriangles();
            for  (std::vector<Triangle>::const_iterator it=interTris.begin(); it!=interTris.end(); ++it){ // interTris.size()>0){

                glBegin(GL_TRIANGLES);
                glColor4f(it->ca.x,it->ca.y,it->ca.z,it->ca.w);
                glNormal3f(it->na.x,it->na.y,it->na.z);
                glVertex3f(it->a.x,it->a.y,it->a.z);   

                glColor4f(it->cb.x,it->cb.y,it->cb.z,it->cb.w);
                glNormal3f(it->nb.x,it->nb.y,it->nb.z);
                glVertex3f(it->b.x,it->b.y,it->b.z);  

                glColor4f(it->cc.x,it->cc.y,it->cc.z,it->cc.w);
                glNormal3f(it->nc.x,it->nc.y,it->nc.z);
                glVertex3f(it->c.x,it->c.y,it->c.z);  

                glEnd();

                glLineWidth(3.0);

                glBegin(GL_LINE_LOOP);
                glColor4f(0,0,0,1);
                glNormal3f(it->na.x,it->na.y,it->na.z);
                glVertex3f(it->a.x,it->a.y,it->a.z);   

                glColor4f(0,0,0,1);
                glNormal3f(it->nb.x,it->nb.y,it->nb.z);
                glVertex3f(it->b.x,it->b.y,it->b.z);  

                glColor4f(0,0,0,1);
                glNormal3f(it->nb.x,it->nb.y,it->nb.z);
                glVertex3f(it->b.x,it->b.y,it->b.z);  

                glColor4f(0,0,0,1);
                glNormal3f(it->nc.x,it->nc.y,it->nc.z);
                glVertex3f(it->c.x,it->c.y,it->c.z);  

                glEnd();        
                glLineWidth(1.0);
            }
        }
    }
    */

    
    const std::vector<Segment>& thicknessRays = dc3d->thicknessRays(showRayThicknessIndex);
    std::vector<Segment>::const_iterator it = thicknessRays.begin();
    
    glBegin(GL_LINES);
    glColor3f(1, 0, 0);
        for (; it!=thicknessRays.end(); ++it) {
          const Segment &r = *it;
          glVertex3f(r.a.x, r.a.y, r.a.z);
          glVertex3f(r.b.x, r.b.y, r.b.z);
        }
    glEnd();

    
    if (showGoodRays){
        for(unsigned i = dc3d->_threshold ; i < dc3d->_maximum ; ++i) {

          const std::set<Segment,classcomp>& gRays = dc3d->goodRays(i);
          std::set<Segment,classcomp>::const_iterator it = gRays.begin();
          // draw rays

          //double color = ((dc3d->_maximum - i)*(0.5))/(dc3d->_maximum-dc3d->_threshold);
          vec3d color(0,0,0);
          if (i == dc3d->maximum()){
              glLineWidth(2.5);
              color = vec3d(0,1,0);
          } else if (i == dc3d->maximum()-1){
              glLineWidth(1.5);
              color = vec3d(1,0,0);
          } else if (i == dc3d->maximum()-2){
              glLineWidth(1.0);
              color = vec3d(0.5,0,1);
          } else {
              glLineWidth(0.5);
              color = vec3d(0,0,1);
          }
          //glLineWidth( (i - dc3d->_threshold)*3 + 1 );
          glBegin(GL_LINES);
          glColor3f(color.x, color.y, color.z);
          unsigned value=0; 
          std::vector<Point> points;
          for (; it!=gRays.end(); ++it) {
                if (value++ == rayIndexT){                    
                    const Segment &r = *it;
                    dc3d->processMeshSegment(r, &points);
                    //if (points.size() != dc3d->_threshold){
                        glVertex3f(r.a.x, r.a.y, r.a.z);
                        glVertex3f(r.b.x, r.b.y, r.b.z);
                     //   std::cout << points.size() << std::endl;
                   // }
                }
          }
          glEnd();
          for (unsigned int t = 0; t < points.size(); t++){
            double inc = (double)t/(double)points.size();
            glColor3f(1.f-inc,0.f+inc,0);
            glPointSize(20.f);
              glBegin(GL_POINTS);
              glVertex3f(points[t].x, points[t].y, points[t].z);
              printf("[%f, %f, %f]\n",points[t].x, points[t].y, points[t].z);
              glEnd();
          }
          std::cout << std::endl;
          points.clear();
          if (thresholdRays) break;
        }
    }
    
     
    if (discretePts){
        glColor3f(0,0,0);
        glPointSize(3.0);
        const std::vector<vec3d> &vpoints = dc3d->visualizationPoints();
        //vec3d center = mesh.aabb.center();
        glBegin(GL_POINTS);
        for(unsigned i=0; i<dc3d->visualizationPoints().size(); i++){
            glVertex3f(vpoints[i].x,vpoints[i].y,vpoints[i].z );
        }
        glEnd();
    }
    
    glBegin(GL_LINES);
    //for(unsigned i=0; i<dc3d->visualizationPoints().size(); i++){
        glVertex3f(0,0,0);
        glVertex3f(1000,0,0);
        
        glVertex3f(0,0,0);
        glVertex3f(0,1000,0);
        
        glVertex3f(0,0,0);
        glVertex3f(0,0,1000);
    //    glVertex3f( (center.x - vpoints[i].x) , (center.y-vpoints[i].y), (center.z-vpoints[i].z) );
    //}
   glEnd();
    
    
    if(showPlanes) {
      const std::vector<Segment>& bounds = dc3d->usedPlanes();
      // draw planes
      glLineWidth(1);
      glBegin(GL_LINES);
      glColor3f(0.5, 0.5, 0.5);
        for (unsigned i=0; i<bounds.size(); ++i) {
          const Segment &r = bounds[i];
          glVertex3f(r.a.x, r.a.y, r.a.z);
          glVertex3f(r.b.x, r.b.y, r.b.z);
        }
      glEnd();
    }

    const std::vector<Point>& points = dc3d->intersectionPoints();
    glEnable(GL_LIGHTING);
    glMaterialfv(GL_FRONT_AND_BACK, GL_DIFFUSE, &sphereColor.x);
    for (unsigned i=0; i<points.size(); ++i)  {
      const Point &p = points[i];
      glPushMatrix();
        glTranslatef(p.x, p.y, p.z);
        glutSolidSphere(radius, 10, 10);
      glPopMatrix();
    }
}

void setupCamera(Camera& camera)
{
    glViewport(0, 0, winWidth, winHeight);
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    camera.setPerspec(50, (double)winWidth/winHeight, 0.01, 500);
   // gluPerspective(50, (double)winWidth/winHeight, 0.1, 1000);

    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();
    //cam.applyTransform();
    camera.update();	
    camera.lookAt();
}


const unsigned int colorTableSize=6;

float ColorTable[colorTableSize][3] = {
    {0.0, 0.0, 1.0}, // Azul
    {1.0, 0.0, 1.0}, // Roxo/Rosa
    {1.0, 0.5, 0.0}, // Laranja
    {1.0, 1.0, 0.0}, // Amarelo
    {0.0, 1.0, 0.0}, // Verde
    {0.0, 1.0, 0.0} // Verde
};

CImg<float> color_scale(int w, int h){
    
    float parts = h/(float)(colorTableSize-2);
   
 
    CImg<float> teste(w+30,h+3,1,3, 0.5);
    
    for (unsigned int i=1; i<colorTableSize-1; i++){
        vec3d inc((ColorTable[i][0] - ColorTable[i-1][0])/parts,
                  (ColorTable[i][1] - ColorTable[i-1][1])/parts,
                  (ColorTable[i][2] - ColorTable[i-1][2])/parts);
        float color[3] = {ColorTable[i-1][0], ColorTable[i-1][1], ColorTable[i-1][2]};
        
        for (int y=(i-1)*parts; y<(i-1)*parts+2; y++){
         for (int x=0; x<w+30; x++){
                teste(x,y,0,0) = 0.2;
                teste(x,y,0,1) = 0.2;
                teste(x,y,0,2) = 0.2;
            }
        }
        
        for (int y=(i-1)*parts+2; y<(i*parts); y++){            
            for (int x=0; x<w+30; x++){
                teste(x,y,0,0) = color[0];
                teste(x,y,0,1) = color[1];
                teste(x,y,0,2) = color[2];
            }
            color[0] += inc.x; color[1] += inc.y; color[2] += inc.z;
        }
    }
    
    for (int y=h; y<h+2; y++){
        for (int x=0; x<w+30; x++){
            teste(x,y,0,0) = 0.2;
            teste(x,y,0,1) = 0.2;
            teste(x,y,0,2) = 0.2;
        }
    }
    //const char *t = "teste";
    //teste.draw_text(10,10,t,0,0);
    return teste;
    
}

void Tex3DTweakBar(){
    // initialize
    dc3d->tex3d.min = dc3d->_threshold;
    dc3d->tex3d.max = dc3d->_maximum;
    dc3d->tex3d.colorAlpha[0] = 0.1f;
    dc3d->tex3d.colorAlpha[1] = 1.0f;
    dc3d->tex3d.colorAlpha[2] = 1.0f;
    dc3d->tex3d.colorAlpha[3] = 1.0f;
    dc3d->tex3d.colorAlpha[4] = 1.0f;
    dc3d->tex3d.colorAlpha[5] = 1.0f;
    
    
    TwBar *bar = TwNewBar("mybar");
    TwDefine(" mybar label='texture3D' position='50 400' size='160 100'");
      
    TwAddVarRW(bar, "Min", TW_TYPE_INT32,  &(dc3d->tex3d.min), "");
    TwSetParam(bar, "Min", "max", TW_PARAM_INT32, 1, &dc3d->_maximum);
    TwSetParam(bar, "Min", "min", TW_PARAM_INT32, 1, &dc3d->_threshold);
    
    TwAddVarRW(bar, "Max", TW_TYPE_INT32,  &(dc3d->tex3d.max), "");
    TwSetParam(bar, "Max", "max", TW_PARAM_INT32, 1, &dc3d->_maximum);
    TwSetParam(bar, "Max", "min", TW_PARAM_INT32, 1, &dc3d->_threshold);
    
    TwAddVarRW(bar, "alpha1", TW_TYPE_FLOAT, &(dc3d->tex3d.colorAlpha[0]), "min=0 max=1 step=0.1");
    TwAddVarRW(bar, "alpha2", TW_TYPE_FLOAT, &(dc3d->tex3d.colorAlpha[1]), "min=0 max=1 step=0.1");
    TwAddVarRW(bar, "alpha3", TW_TYPE_FLOAT, &(dc3d->tex3d.colorAlpha[2]), "min=0 max=1 step=0.1");
    TwAddVarRW(bar, "alpha4", TW_TYPE_FLOAT, &(dc3d->tex3d.colorAlpha[3]), "min=0 max=1 step=0.1");
    TwAddVarRW(bar, "alpha5", TW_TYPE_FLOAT, &(dc3d->tex3d.colorAlpha[4]), "min=0 max=1 step=0.1");
    TwAddVarRW(bar, "alpha6", TW_TYPE_FLOAT, &(dc3d->tex3d.colorAlpha[5]), "min=0 max=1 step=0.1");
    

}


void recompute(void *data)
{
    //const TriMesh* mesh = reinterpret_cast<const TriMesh*>(data);
    TriMesh* mesh = reinterpret_cast<TriMesh*>(data);
    //if(doGoodRays)
    //  dc3d->setComputeGoodRays(true);
    //else
    //  dc3d->setComputeGoodRays(false);
    //std::cout << "mehs faces: " << mesh->faces.size() << std::endl;
    
    
    tic();    
    dc3d->process(*mesh);
    toc("Depth Complexity");
    
    std::clog << std::endl << "=== RESULTS ===" << std::endl;
    std::clog << "  --> NUM SAMPLES: " << dc3d->visualizationPoints().size() << std::endl;
    std::clog << "  --> DC MAX: " << dc3d->maximum() << std::endl;
    std::clog << "  --> NUM MAXIMUM RAYS: " << dc3d->maximumRays().size() << std::endl;
    //if(doGoodRays) {
    unsigned numRays = 0;
    for(unsigned i = dc3d->getThreshold() ; i <= dc3d->maximum() ; ++i)
        numRays += dc3d->goodRays(i).size();
    std::clog << "  --> NUM GOOD RAYS: " << numRays << std::endl;
    
    std::cout << std::endl;
    std::cout << "tri " << mesh->faces.size() << "\n";
    std::cout << "vert " << mesh->vertices.size(); 
    
    //}
    
    // get 2D histogram
    //dc3d->display2dHistogram();
    
      
      /*
    // create color scale
    unsigned int val = dc3d->maximum();
    
    if (showScale){
        cscale = color_scale(40,500); 

        float color[3] = {0.0, 0.0, 0.0};
        cscale.draw_text(3, 2  , "%d", color, 0, 1, 18, val*0/4);
        cscale.draw_text(3, 127, "%d", color, 0, 1, 18, val*1/4);
        cscale.draw_text(3, 252, "%d", color, 0, 1, 18, val*2/4);
        cscale.draw_text(3, 377, "%d", color, 0, 1, 18, val*3/4);
        cscale.draw_text(3, 482, "%d", color, 0, 1, 18, val*4/4);

        main_disp.display(cscale);
    }
    */
    
    
    //Tex3DTweakBar();
    
    // check interception
    //const std::list<unsigned int>& ilist = dc3d->intersectionTris();
    //std::list<unsigned int>::const_iterator it = ilist.begin();
    //for (;it!=ilist.end(); ++it){
    //    mesh->faces[*it].intercepted = true;
    //}
}

void
drawBackground(const vec3f& top, const vec3f& mid, const vec3f& bot)
{
    glDisable(GL_DEPTH_TEST);
    //glDepthMask(GL_FALSE);
    glDisable(GL_LIGHTING);
    glDisable(GL_BLEND);
    
    glMatrixMode(GL_MODELVIEW);
    glPushMatrix();
    glLoadIdentity();
    glMatrixMode(GL_PROJECTION);
    glPushMatrix();
    glLoadIdentity();
    glOrtho(-1, 1, -1, 1, -1, 1);

    glBegin(GL_QUAD_STRIP);
    glColor3fv(&bot.x);
    glVertex2f(-1, -1);
    glVertex2f( 1, -1);
    glColor3fv(&mid.x);
    glVertex2f(-1, 0);
    glVertex2f( 1, 0);
    glColor3fv(&top.x);
    glVertex2f(-1, 1);
    glVertex2f( 1, 1);
    glEnd();

    glPopMatrix();
    glMatrixMode(GL_MODELVIEW);
    glPopMatrix();

    glDepthMask(GL_TRUE);
    glEnable(GL_DEPTH_TEST);
}

//________________________________________________ CAMERA
Camera camera;

//______________ MOUSE
double mDx = 0;
double mDy = 0;
bool mclicked;
int mbutton;

//_______________________________________________________________ Called when a mouse button is pressed or released
void mouse_click(GLFWwindow* w,int button, int action, int mod){
	if(TwEventMouseButtonGLFW(button,action))
		return;
	glfwGetCursorPos(window, &mDx,&mDy);
	mclicked=action==GLFW_PRESS;
	mbutton=button;
}

//_______________________________________________________________ Called when a mouse move and a button is pressed
void mouse_motion(GLFWwindow* w, double x, double y){
	if(TwEventMousePosGLFW(x,y))
		return;
	if( mclicked ){
		float dx = mDx - x;
		float dy = mDy - y;
		float d = sqrt(dx*dx + dy*dy);
		float delta = 0.0015;

		switch(mbutton){
			case GLFW_MOUSE_BUTTON_RIGHT : //look
				//if( glutGetModifiers() == GLUT_ACTIVE_SHIFT){
					camera.lookLefRigObj(dx * delta);
					camera.lookUpDownObj(dy * delta);
				/*} else { 
					camera.lookLefRig(dx * delta);
					camera.lookUpDown(dy * delta);
				}*/
			break;
			case GLFW_MOUSE_BUTTON_LEFT:			
				//if( glutGetModifiers() != GLUT_ACTIVE_SHIFT){
					camera.MoveLado(dx*delta);
					camera.MoveCimaBaixo(-dy*delta);
				/*}else{
					camera.MoveLadoObj(dx*delta);
					camera.MoveCimaBaixoObj(-dy*delta);				
				//}*/
			break;
                    case GLFW_MOUSE_BUTTON_MIDDLE :
				if( dy>0) d = -d;
				//if( glutGetModifiers() != GLUT_ACTIVE_SHIFT)
					camera.MoveFrente(d*delta);
				/*else
					camera.MoveFrenteObj(d*delta);*/
			break;
		}//end switch
		mDx = x;	mDy = y;
	}//end if
}

void initializeGraphics(){
    //std::cout << dc3d << std::endl;
    glfwInit();
    std::atexit(glfwTerminate);
    
    //GLFWvidmode mode;
    
    //glfwGetDesktopMode(&mode);
#ifdef COMPUTE_OFFLINE
    glfwWindowHint(GLFW_VISIBLE, GL_FALSE);
#endif
    window = glfwCreateWindow(1024, 768, "Shape Retrieval", NULL, NULL);
    
    if(!window)
    {
        std::cerr << "failed to open window!" << std::endl;
        glfwTerminate();
        return;
    }
    
    /* Make the window's context current */
    glfwMakeContextCurrent(window);
    
    // initialize GLEW
    GLenum glewStatus = glewInit();
    if (glewStatus != GLEW_OK){
      std::cerr << "[ERROR] "<< glewGetErrorString(glewStatus)<<std::endl;
    }
    
    // Print OPENGL, SHADER and GLEW versions
    std::clog << std::endl << "=== SYSTEM INFO ===" << std::endl;
    std::clog << "  --> OPENGL VERSION: " << glGetString(GL_VERSION) << std::endl;
    std::clog << "  --> SHADER VERSION: " << glGetString(GL_SHADING_LANGUAGE_VERSION) << std::endl;
    std::cerr << "  --> GLEW VERSION: "<<glewGetString(GLEW_VERSION)<<std::endl;
    //std::clog << "----------------" std::endl;
 
    std::clog << std::endl << "=== DEBUG ===" << std::endl;

    //glfwEnable(GLFW_MOUSE_CURSOR);
    //glfwEnable(GLFW_KEY_REPEAT);
    //glfwSetWindowTitle("Plane-Triangle intersection test");

	
    // Initialize AntTweakBar
    if( !TwInit(TW_OPENGL, NULL) )
    {
        std::cerr << "AntTweakBar initialization failed: " << TwGetLastError() << std::endl;
        return;
    }
    TwWindowSize(1024, 768);
    
    // Set GLFW event callbacks
    glfwSetWindowSizeCallback(window, WindowSizeCB);

    glfwSetMouseButtonCallback(window,mouse_click);
    glfwSetCursorPosCallback(window,mouse_motion);
    glfwSetScrollCallback(window,(GLFWscrollfun)TwEventMouseWheelGLFW);
    glfwSetKeyCallback(window, (GLFWkeyfun)TwEventKeyGLFW);
    glfwSetCharCallback(window, (GLFWcharfun)TwEventCharGLFW);
    
}

int doInteractive(TriMesh& mesh)
{
    
    // setup opengl, glut, glew stuff.
    //initializeGraphics();    

    TwBar *bar = TwNewBar("Controls");
    TwDefine(" GLOBAL ");
	
    vec3f top(0.25, 0.25, .5), mid(0.75, 0.75, .85), bot(1, 1, 1);

    vec4f objdiff(0.55, 0.5, 0, 0.5), objspec(.75, .75, .75, .2);
    GLfloat shine = 50;
    bool showObj = true;
    
    TwAddVarRW(bar, "showPlanes", TW_TYPE_BOOLCPP, &showPlanes, " label='show discret. planes' ");

    TwAddVarRW(bar, "goodRays", TW_TYPE_BOOLCPP, &doGoodRays, " label='show more rays' ");
    TwAddVarRW(bar, "discrete pts", TW_TYPE_BOOLCPP, &discretePts, " label='disc. pts' ");

    TwAddVarRW(bar, "discretSteps", TW_TYPE_UINT32, &dc3d->_discretSteps, " label='discret. steps' min=2 step=5 ");
    TwAddVarRW(bar, "resolution", TW_TYPE_UINT32, &dc3d->_resolution, " label='resolution' min=30 max=700");
    TwAddVarRO(bar, "maxDepth", TW_TYPE_UINT32, &dc3d->_maximum, " label='Max. depth' ");
    
    
    TwAddButton(bar, "recompute", recompute, (void*)&mesh, " label='Recompute' ");
    
    TwAddVarRW(bar, "Compute Rays", TW_TYPE_BOOLCPP, &dc3d->_computeRays, "group='Rays' label='Compute Rays'");
    TwAddVarRW(bar, "goodThreshold", TW_TYPE_UINT32, &dc3d->_threshold, "group='Rays' label='threshold' min=0 ");
    TwAddVarRW(bar, "limit rays", TW_TYPE_UINT32, &dc3d->_limitRays, " group='Rays' label='Limit Rays'");
    TwAddVarRW(bar, "Show Max", TW_TYPE_BOOLCPP, &showMaxRays, " group='Rays' label='Show Max'");
    TwAddVarRW(bar, "Show Good", TW_TYPE_BOOLCPP, &showGoodRays, " group='Rays' label='Show Good'");
    TwAddVarRW(bar, "Threshold Rays", TW_TYPE_BOOLCPP, &thresholdRays, " group='Rays' label='Threshold Rays'");
    TwAddVarRW(bar, "Good Rays", TW_TYPE_UINT32, &rayIndexT, " group='Rays' label='Good rays Index'");

    
    TwAddVarRW(bar, "highlight Tris", TW_TYPE_BOOLCPP, &highlightTris, " group='Rays' label='Highlight Tris'");
    
    TwAddVarRW(bar, "thickness0", TW_TYPE_BOOLCPP, &dc3d->_computeThicknessFlag, "group='Thickness' label='Compute Thickness'");
    TwAddVarRW(bar, "thickness1", TW_TYPE_UINT32, &showRayThicknessIndex, " min=0 max=999 group='Thickness' label='thickness index'");    

    TwAddVarRW(bar, "top", TW_TYPE_COLOR3F, &top.x, " group='background' ");
    TwAddVarRW(bar, "mid", TW_TYPE_COLOR3F, &mid.x, " group='background' ");
    TwAddVarRW(bar, "bot", TW_TYPE_COLOR3F, &bot.x, " group='background' ");

    TwAddVarRW(bar, "Diff", TW_TYPE_COLOR4F, &objdiff.x, " group='Object' ");
    TwAddVarRW(bar, "Spec", TW_TYPE_COLOR4F, &objspec.x, " group='Object' ");
    TwAddVarRW(bar, "Shin", TW_TYPE_FLOAT, &shine, " group='Object' min='1' max='128' ");
    TwAddVarRW(bar, "Show", TW_TYPE_BOOLCPP, &showObj, " group='Object' ");
		
    TwAddVarRW(bar, "radius", TW_TYPE_FLOAT, &radius, "  group='Sphere' label='radius' min=0.0 step=0.001 max=2.0");
    TwAddVarRW(bar, "Color", TW_TYPE_COLOR4F, &sphereColor.x, " group='Sphere' ");
    //TwAddVarRW(bar, "Show Ray", TW_TYPE_UINT32, &showRayIndex, " group='rays' min=0");

    glEnable(GL_LIGHTING);
    glEnable(GL_LIGHT0);
    glDisable(GL_CULL_FACE);
    glEnable(GL_DEPTH_TEST);

    //Camera cam;
    BoundingBox aabb = mesh.aabb;
    //std::cout << "aabb: " << aabb.center() << std::endl;
    
    camera.bbox(float3(aabb.min.x,aabb.min.y,aabb.min.z), float3(aabb.max.x,aabb.max.y,aabb.max.z), true );
		camera.front();
    
    //cam.target = aabb.center();
    //cam.up = vec3f(0, 1, 0);
    //cam.pos = cam.target + vec3f(0, 0, 2*aabb.extents().z);
    while( !glfwWindowShouldClose(window) && !glfwGetKey(window,GLFW_KEY_ESCAPE) ) {
        glClear( GL_COLOR_BUFFER_BIT|GL_DEPTH_BUFFER_BIT );

        drawBackground(top, mid, bot);

        setupCamera(camera);
        
        #ifdef COMPUTE_OFFLINE
            //std::cout<<"OFFLINE RENDERING --" << std::endl;
            //std::cout << "before" << std::endl;
            recompute((void*)&mesh);
            //std::cout << "finalize" << std::endl;
            break;
            
        #endif
       
        /*camera.update();	
	    camera.lookAt();*/

        /*GLfloat lpos[4] = { camera.GetEye().x, camera.GetEye().y, camera.GetEye().z, 1 };
        glLightfv(GL_LIGHT0, GL_POSITION, lpos);*/
       
        //std::cout << "not\n";
        //drawPlaneIntersection(intersectionVectors);
        drawRays();
        



        //glMaterialfv(GL_FRONT_AND_BACK, GL_DIFFUSE, &objdiff.x);
        glMaterialfv(GL_FRONT_AND_BACK, GL_SPECULAR, &objspec.x);
        glMaterialf(GL_FRONT_AND_BACK, GL_SHININESS, shine);
        if (showObj) drawMesh(mesh, vec3f(camera.GetDir().x,camera.GetDir().y,camera.GetDir().z));

        Triangle *t = &mesh.faces[550];
        t->ca = vec4d(1.0f, 0.0f, 0.0f, 0.7f);
        t->cb = vec4d(1.0f, 0.0f, 0.0f, 0.7f);
        t->cc = vec4d(1.0f, 0.0f, 0.0f, 0.7f);

        t = &mesh.faces[608];
        t->ca = vec4d(1.0f, 0.0f, 0.0f, 0.7f);
        t->cb = vec4d(1.0f, 0.0f, 0.0f, 0.7f);
        t->cc = vec4d(1.0f, 0.0f, 0.0f, 0.7f);
        
        

        // Draw tweak bars
        TwDraw();

        // Present frame buffer
        glfwSwapBuffers(window);
        glfwPollEvents();
        
        
    }
                

    glfwTerminate(); 
    return 0;
}

std::string getExtension(const std::string& filename)
{
    std::string::size_type dotpos = filename.rfind(".");
    if (dotpos != std::string::npos)
        return filename.substr(dotpos+1);
    return "";
}

#ifndef COMPUTE_OFFLINE
int main(int argc, char **argv)
{
    std::clog << "*----- ONLINE RENDERING -----*" << std::endl;
    if (argc == 1) {
        std::cerr << "[ERROR] Missing Input File!" << std::endl;
        return 1;
    }
   
 
    glutInit(&argc,argv);
    
    try {
        std::ifstream file(argv[1]);
        std::string ext = getExtension(argv[1]);

        TriMesh mesh;

       if (ext == "off" || ext == "OFF")
            mesh = loadOFFMesh(file);
        else if (ext == "obj" || ext == "OBJ")
            mesh = loadOBJMesh(file);
        else
            throw "Unknown file type!";
       file.close();
        
        // setup opengl, glut, glew stuff.
        initializeGraphics();       
        dc3d = new RFDepthComplexity3D(128, 50);
        
        if (doInteractive(mesh))
            return 1;

    }
    catch (const char* msg)  {
        std::cerr << "Failed: " << msg << std::endl;
        return 1;
    }
    catch (std::string msg) {
        std::cerr << "Failed: " << msg << std::endl;
    }
}

#else

int main (int argc, char **argv) {
  cmd_usage("Program to find depth complexity out (offline mode)");
  
  const char *filename = cmd_option("-f", "models/suzanne.obj", "Model in OBJ or OFF format.");
  const char *file_basename = cmd_option("-fb", "hist2D", "the basename of the filename");
  //const int fboWidth  = cmd_option("-fboWidth",  512, "Framebuffer width.");
  //const int fboHeight = cmd_option("-fboHeight", 512, "Framebuffer height.");
  const int resolution = cmd_option("-res", 512, "Framebuffer Resolution.");  
  const int discretSteps = cmd_option("-dsteps", 50, "Discrete steps.");
  const char *filenameParallelData = cmd_option("-fpr", "", "Save a *.txt containing rays informaton");
  const char *filenameHistogram = cmd_option("-fh", "", "Save a *.txt file with histogram information");
  const char *filenameThicknessHistogram = cmd_option("-fth", "", "Save a *.txt file with thickness histogram information");
  const char *filename2DHistogram = cmd_option("-f2dh", "", "Save a *.js file with a 2D histogram");
  //const char *filenameRays = cmd_option("-fr", "", "Save a *.off file with rays in ");
  const char *filenameRaysSpherical = cmd_option("-frs", "", "Save a *.txt file with rays in spherical coordinates");
  const int sphericalThreshold = cmd_option("-k",  0, "Spherical coordinates are calculated for rays with DC between MDC-k and MDC");
  //const int rotateX = cmd_option("-x",  0, "Rotate model around x-axis x degrees");
  //const int rotateY = cmd_option("-y",  0, "Rotate model around y-axis x degrees");
  //const int rotateZ = cmd_option("-z",  0, "Rotate model around z-axis x degrees");
  const char *filenameRays = cmd_option("-fr", "", "Save a *.off file with rays");
  const bool computeMoreRays = cmd_option("-cmr", strcmp(filenameRaysSpherical, "")!=0, "Whether rays above the threshold of intersections should be output");
  //const bool computeMoreRays = cmd_option("-cmr", false, "Whether rays above the threshold of intersections should be output");
  const int intersectionThreshold  = cmd_option("-it",  0, "Threshold of intersections");

  //std::cout << "Threshold: " << intersectionThreshold << std::endl;
  try {
    glutInit(&argc, argv);
    //tic();
    std::ifstream file(filename);
    std::string ext = getExtension(filename);
    TriMesh mesh;
    if (ext == "off" || ext == "OFF")
      mesh = loadOFFMesh(file);
    else if (ext == "obj" || ext == "OBJ")
      mesh = loadOBJMesh(file);
    else
      throw "Unknown file type!";
    file.close();
    
    //rotatePoint(mesh.aabb.min ,rotateX,rotateY,rotateZ);
    //rotatePoint(mesh.aabb.max ,rotateX,rotateY,rotateZ); 
    //toc("Loading Mesh");

    // TODO(jpocom) Look for another way to call glewInit() and work in offscreen.
    // found and working =)
    
    
    std::cout << "*---- RENDERING OFFLINE ----*" << std::endl;
    // setup opengl, glut, glew stuff.
    initializeGraphics(); 
    dc3d = new RFDepthComplexity3D(512, 50);    

    dc3d->setDiscreteSteps(discretSteps);
    dc3d->setResolution(resolution);
    dc3d->setComputeHistogram(strcmp(filenameHistogram, "")!=0);
    dc3d->setComputeMaximumRays(strcmp(filenameRays, "")!=0);
    dc3d->setComputeGoodRays(false);
    dc3d->setThreshold(intersectionThreshold);
    
    //std::cout << "parameters: " << resolution << << std::endl;

    //tic();
    doInteractive(mesh);
    //toc("Computing Depth Complexity");

    //std::cout << "Maximum: " << dc3d->maximum() << std::endl;
    
    // Saving Histogram file
    if (strcmp(filenameHistogram, "")!=0) {
      std::string extTxt = getExtension(filenameHistogram);
      if (extTxt == "js" || extTxt == "JS") {
        std::ofstream fileHistogram(filenameHistogram);
        dc3d->writeHistogram(fileHistogram, file_basename);
        fileHistogram.close();
      } else throw "Histogram's file should be *.txt!";
    }
    
    
     // Saving Thickness Histogram file
    if (strcmp(filenameThicknessHistogram, "")!=0) {
      std::string extTxt = getExtension(filenameThicknessHistogram);
      if (extTxt == "js" || extTxt == "JS") {
        std::ofstream fileThicknessHistogram(filenameThicknessHistogram);
        dc3d->writeThicknessHistogram(fileThicknessHistogram,file_basename);
        fileThicknessHistogram.close();
      } else throw "Histogram's file should be *.txt!";
    }

    // Saving 2D Histograma JS
    if (strcmp(filename2DHistogram, "")!=0) {
      //std::string extTxt = getExtension(filename2DHistogram);
      //if (extTxt == "js" || extTxt == "JS") {
        std::string basename(filename2DHistogram);
        basename+=".js";
        std::ofstream file2DHistogram(basename.c_str());
        dc3d->write2dHistogram(file2DHistogram, file_basename);
        file2DHistogram.close();
      //} else throw "Histogram's file should be *.js!";
    }

    // Saving Rays file
    if (strcmp(filenameRays, "")!=0) {
      std::string extOff = getExtension(filenameRays);
      if (extOff == "off" || extOff == "OFF") {
        std::ofstream fileRays(filenameRays);
        dc3d->writeRays(fileRays);
        fileRays.close();
      } else throw "Ray's file should be *.off!";
    }
   
    // Saving MoreRays file
    if(computeMoreRays) {
      unsigned numRays = 0;
      for(unsigned i = dc3d->getThreshold() ; i <= dc3d->maximum() ; ++i)
        numRays += dc3d->goodRays(i).size();
      std::cout << "Number of good rays: " << numRays << std::endl;
      for(unsigned i=intersectionThreshold;i<=dc3d->maximum();++i) {
        std::ostringstream filenameMoreRays;
        filenameMoreRays << "rays" << i << ".off";
        std::ofstream fileRays(filenameMoreRays.str().c_str());
        dc3d->writeRays(fileRays,dc3d->goodRays(i),i);
        fileRays.close();
      }
    }
    
     // Saving RaysSpherical file
    if (strcmp(filenameRaysSpherical, "")!=0) {
      std::string extTxt = getExtension(filenameRaysSpherical);
      if (extTxt == "txt" || extTxt == "TXT") {
        std::ofstream fileRaysSpherical(filenameRaysSpherical);
        dc3d->writeRaysSpherical(fileRaysSpherical,sphericalThreshold);
        fileRaysSpherical.close();
      } else throw "Spherical Rays' file should be *.txt!";
    }
    
    if (strcmp(filenameParallelData, "")!=0) {
        std::string extTxt = getExtension(filenameRaysSpherical);
        if (extTxt == "txt" || extTxt == "TXT") {
            std::ofstream fileRaysSpherical(filenameRaysSpherical);
            //dc3d->writeRaysSpherical(fileRaysSpherical,sphericalThreshold);
            fileRaysSpherical.close();        
        }
    }
     
      
  } catch (const char* msg)  {
    std::cerr << "Failed: " << msg << std::endl;
    return 1;
  } catch (std::string msg) {
    std::cerr << "Failed: " << msg << std::endl;
    return 1;
  }

  return 0;
}

#endif
