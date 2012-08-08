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
#include <GL/glfw.h>
#include <AntTweakBar.h>
#ifdef __APPLE__
#include <GLUT/glut.h>
#else
#include <GL/glut.h>
#endif


#include "CImg.h"

#include "flags.h"
#include "vector.hpp"
#include "camera/float3.h"
//#include "camera.hpp"
#include "camera/Camera.h"
#include "util.h"
#ifdef USE_RANDOM_DC3D
#include "dc_3d_random.h"
#else
#include "dc_3d.h"
#endif
#include "timer.h"

using namespace cimg_library;

float radius = 0.01;
bool showPlanes = false;
bool doGoodRays = true;
vec4f sphereColor(0.0, 1.0, 0.0, 1.0);
unsigned int showRayIndex = 0;

CImgDisplay main_disp;
CImg<float> cscale;//color_scale(40,500); 

#ifdef USE_RANDOM_DC3D
RDepthComplexity3D *dc3d;
const char *filenameRays;
#else
DepthComplexity3D *dc3d;
#endif

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
void GLFWCALL WindowSizeCB(int width, int height)
{
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
#ifndef USE_RANDOM_DC3D
    glBindTexture(GL_TEXTURE_3D, dc3d->getTextureID());
#endif
    glColorMaterial(GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE);
    glEnable(GL_COLOR_MATERIAL);
    
    //glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
    //std::clog << "sorting...";
    if (sorted_faces.empty()) {
        sorted_faces = mesh.faces;
    }
    
    std::sort(sorted_faces.begin(), sorted_faces.end(), ByDist(dir));
    //std::clog << "done" << std::endl;
    
    glEnable(GL_BLEND);
    glBlendFunc(GL_ONE, GL_ONE_MINUS_SRC_ALPHA);
    //glDisable(GL_DEPTH_TEST);
    glEnable(GL_DEPTH_TEST);
    glEnable(GL_LIGHTING);
    glLightModeli(GL_LIGHT_MODEL_TWO_SIDE, 1);

    glEnable(GL_VERTEX_ARRAY);

    glEnableClientState(GL_VERTEX_ARRAY);    
    glVertexPointer(3, GL_DOUBLE, 3*sizeof(vec3d)+sizeof(vec4d), &sorted_faces[0].a.x);
    
    glEnableClientState (GL_TEXTURE_COORD_ARRAY);
    glTexCoordPointer(3, GL_DOUBLE, 3*sizeof(vec3d)+sizeof(vec4d), &sorted_faces[0].tca.x);
    
    glEnableClientState(GL_COLOR_ARRAY);
    glColorPointer(4, GL_DOUBLE, 3*sizeof(vec3d)+sizeof(vec4d), &sorted_faces[0].ca.x);

    glEnableClientState(GL_NORMAL_ARRAY);
    glNormalPointer(GL_DOUBLE, 3*sizeof(vec3d)+sizeof(vec4d), &sorted_faces[0].na.x);

    glDrawArrays(GL_TRIANGLES, 0, sorted_faces.size()*3);

    glDisableClientState(GL_VERTEX_ARRAY);
    glDisable(GL_VERTEX_ARRAY);
    glDisableClientState(GL_NORMAL_ARRAY);
    glDisable(GL_NORMAL_ARRAY);
    glDisableClientState(GL_COLOR_ARRAY);
    glDisable(GL_COLOR_ARRAY);

    glDisable(GL_BLEND);
    glEnable(GL_DEPTH_TEST);
    //glDisable(GL_COLOR_MATERIAL);
    
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
    
void drawRays()
{
    if(!doGoodRays) {
        const std::set<Segment,classcomp>& rays = dc3d->maximumRays();
        std::set<Segment,classcomp>::const_iterator it = rays.begin();
        // draw rays
        glLineWidth(1);
        glBegin(GL_LINES);
        glColor3f(0.5, 0.0, 0.5);
          for (; it!=rays.end(); ++it) {
            const Segment &r = *it;
            glVertex3f(r.a.x, r.a.y, r.a.z);
            glVertex3f(r.b.x, r.b.y, r.b.z);
          }
        glEnd();
    }
    
    for(unsigned i = dc3d->_threshold ; i <= dc3d->_maximum && doGoodRays ; ++i) {
      const std::set<Segment,classcomp>& gRays = dc3d->goodRays(i);
      std::set<Segment,classcomp>::const_iterator it = gRays.begin();
      // draw rays
      double color = ((dc3d->_maximum - i)*(0.5))/(dc3d->_maximum-dc3d->_threshold);
      glLineWidth( (i - dc3d->_threshold)*3 + 1 );
      glBegin(GL_LINES);
      glColor3f(1.0 - color, color, color);
        for (; it!=gRays.end(); ++it) {
          const Segment &r = *it;
          glVertex3f(r.a.x, r.a.y, r.a.z);
          glVertex3f(r.b.x, r.b.y, r.b.z);
        }
      glEnd();
    }
#ifndef USE_RANDOM_DC3D
    if(showPlanes) {
      const std::vector<Plane>& bounds = dc3d->usedPlanes();
      // draw planes
      glEnable(GL_BLEND);
      glBlendFunc(GL_ONE, GL_ONE_MINUS_SRC_ALPHA);

      glColor4f(0.5f, 0.45f, 0.3f, 0.34f);
      for (unsigned int i=0; i<bounds.size(); i++){
          const Plane &p = bounds[i];
          glBegin(GL_QUADS);
                glVertex3f(p.a.x, p.a.y, p.a.z);
                glVertex3f(p.b.x, p.b.y, p.b.z);
                glVertex3f(p.c.x, p.c.y, p.c.z);
                glVertex3f(p.d.x, p.d.y, p.d.z);
          glEnd();
      }
      glDisable(GL_BLEND);

      glLineWidth(4);
      const std::vector<Segment>& bounds2 = dc3d->usedPlanes2();
      
      glBegin(GL_LINES);
      //glColor3f(0.5, 0.5, 0.5);
        for (unsigned i=0; i<bounds2.size(); ++i) {
          const Segment &r = bounds2[i];
          if (i%2 == 0 )  glColor4f(0.2, 0.9, 0.1, 1.0);
          else            glColor4f(0.9, 0.3, 0.1, 1.0);
         
          glVertex3f(r.a.x, r.a.y, r.a.z);         
          glVertex3f(r.b.x, r.b.y, r.b.z);
        }
      glEnd();

     }
#endif

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
    camera.setPerspec(50, (double)winWidth/winHeight, 0.1, 10000);
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

void recompute(void *data)
{
    const TriMesh* mesh = reinterpret_cast<const TriMesh*>(data);
    if(doGoodRays)
      dc3d->setComputeGoodRays(true);
    else
      dc3d->setComputeGoodRays(false);
    tic();
    dc3d->process(*mesh);
    toc("Depth Complexity");
    std::clog << "Maximum: " << dc3d->maximum() << std::endl;
    if(doGoodRays) {
      unsigned numRays = 0;
      for(unsigned i = dc3d->getThreshold() ; i <= dc3d->maximum() ; ++i)
        numRays += dc3d->goodRays(i).size();
      std::clog << "Number of good rays: " << numRays << std::endl;
    }
    
    // create color scale
    unsigned int val = dc3d->maximum();
    cscale = color_scale(40,500); 

    float color[3] = {0.0, 0.0, 0.0};
    cscale.draw_text(3, 2  , "%d", color, 0, 1, 18, val*0/4);
    cscale.draw_text(3, 127, "%d", color, 0, 1, 18, val*1/4);
    cscale.draw_text(3, 252, "%d", color, 0, 1, 18, val*2/4);
    cscale.draw_text(3, 377, "%d", color, 0, 1, 18, val*3/4);
    cscale.draw_text(3, 482, "%d", color, 0, 1, 18, val*4/4);
    
    main_disp.display(cscale);
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
int mDx = 0;
int mDy = 0;
bool mclicked;
int mbutton;

//_______________________________________________________________ Called when a mouse button is pressed or released
void GLFWCALL mouse_click(int button, int action){
	if(TwEventMouseButtonGLFW(button,action))
		return;
	glfwGetMousePos(&mDx,&mDy);
	mclicked=action==GLFW_PRESS;
	mbutton=button;
}

//_______________________________________________________________ Called when a mouse move and a button is pressed
void GLFWCALL mouse_motion(int x, int y){
	if(TwEventMousePosGLFW(x,y))
		return;
	if( mclicked ){
		float dx = mDx - x;
		float dy = mDy - y;
		float d = sqrt(dx*dx + dy*dy);
		float delta = 0.001;

		switch(mbutton){
			case GLFW_MOUSE_BUTTON_MIDDLE : //look
				//if( glutGetModifiers() == GLUT_ACTIVE_SHIFT){
					camera.lookLefRigObj(dx * delta );
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
			case GLFW_MOUSE_BUTTON_RIGHT: 
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




int doInteractive(TriMesh& mesh)
{
   
    glfwInit();
    std::atexit(glfwTerminate);
    
    GLFWvidmode mode;
    
    glfwGetDesktopMode(&mode);
    if( !glfwOpenWindow(1024, 768, mode.RedBits, mode.GreenBits, mode.BlueBits, 
                        0, 16, 0, GLFW_WINDOW) )
    {
        std::cerr << "failed to open window!" << std::endl;
        return -1;
    }
    
    // initialize GLEW
    GLenum glewStatus = glewInit();
    if (glewStatus != GLEW_OK){
      std::cerr << "[ERROR] "<< glewGetErrorString(glewStatus)<<std::endl;
    }
    
    // Print OPENGL, SHADER and GLEW versions
    std::clog << "----- << VERSION >>\n";
    std::clog << "OPENGL VERSION: " << glGetString(GL_VERSION) << std::endl;
    std::clog << "SHADER VERSION: " << glGetString(GL_SHADING_LANGUAGE_VERSION) << std::endl;
    std::clog << "GLEW VERSION: "<<glewGetString(GLEW_VERSION)<<std::endl;
    std::clog << "---------------- \n";
 

    glfwEnable(GLFW_MOUSE_CURSOR);
    glfwEnable(GLFW_KEY_REPEAT);
    glfwSetWindowTitle("Plane-Triangle intersection test");

	
    // Initialize AntTweakBar
    if( !TwInit(TW_OPENGL, NULL) )
    {
        std::cerr << "AntTweakBar initialization failed: " << TwGetLastError() << std::endl;
        return 1;
    }
    // Set GLFW event callbacks
    glfwSetWindowSizeCallback(WindowSizeCB);

    glfwSetMouseButtonCallback(mouse_click);
    glfwSetMousePosCallback(mouse_motion);
    glfwSetMouseWheelCallback((GLFWmousewheelfun)TwEventMouseWheelGLFW);
    glfwSetKeyCallback((GLFWkeyfun)TwEventKeyGLFW);
    glfwSetCharCallback((GLFWcharfun)TwEventCharGLFW);

    TwBar *bar = TwNewBar("Controls");
    TwDefine(" GLOBAL ");
	
    vec3f top(0.25, 0.25, .5), mid(0.75, 0.75, .85), bot(1, 1, 1);

    vec4f objdiff(0.55, 0.5, 0, 0.5), objspec(.75, .75, .75, .2);
    GLfloat shine = 50;
    bool showObj = true;
    
    #ifdef USE_RANDOM_DC3D
		if (strcmp(filenameRays, "")!=0)
			dc3d = new RDepthComplexity3D(512, 512, 2, filenameRays);
		else
			dc3d = new RDepthComplexity3D(512, 512, 2);
    #else
    dc3d = new DepthComplexity3D(512, 512, 2);
    #endif
    dc3d->setComputeMaximumRays(true);
    dc3d->setComputeHistogram(true);
    dc3d->setThreshold(10);
    
    TwAddVarRW(bar, "showPlanes", TW_TYPE_BOOLCPP, &showPlanes, " label='show discret. planes' ");

    TwAddVarRW(bar, "goodRays", TW_TYPE_BOOLCPP, &doGoodRays, " label='show more rays' ");

    TwAddVarRW(bar, "goodThreshold", TW_TYPE_UINT32, &dc3d->_threshold, " label='intersection threshold' min=0 ");

    TwAddVarRW(bar, "discretSteps", TW_TYPE_UINT32, &dc3d->_discretSteps, " label='discret. steps' min=2 ");

    TwAddButton(bar, "recompute", recompute, (void*)&mesh, " label='Recompute' ");

    TwAddVarRO(bar, "maxDepth", TW_TYPE_UINT32, &dc3d->_maximum, " label='Max. depth' ");

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
    
    camera.bbox(float3(aabb.min.x,aabb.min.y,aabb.min.z), float3(aabb.max.x,aabb.max.y,aabb.max.z), true );
		camera.front();
    
    //cam.target = aabb.center();
    //cam.up = vec3f(0, 1, 0);
    //cam.pos = cam.target + vec3f(0, 0, 2*aabb.extents().z);
                

    


    while( glfwGetWindowParam(GLFW_OPENED) && !glfwGetKey(GLFW_KEY_ESC) ) {
        glClear( GL_COLOR_BUFFER_BIT|GL_DEPTH_BUFFER_BIT );

        drawBackground(top, mid, bot);

        setupCamera(camera);
        
        /*camera.update();	
	    camera.lookAt();*/

        /*GLfloat lpos[4] = { camera.GetEye().x, camera.GetEye().y, camera.GetEye().z, 1 };
        glLightfv(GL_LIGHT0, GL_POSITION, lpos);*/
       

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
        glfwSwapBuffers();
    }

    return 0;
}

std::string getExtension(const std::string& filename)
{
    std::string::size_type dotpos = filename.rfind(".");
    if (dotpos != std::string::npos)
        return filename.substr(dotpos+1);
    return "";
}


int main(int argc, char **argv)
{
    if (argc == 1) {
        std::cerr << "[ERROR] Missing Input File!" << std::endl;
        return 1;
    }
 
    glutInit(&argc,argv);
    
    try {
        std::ifstream file(argv[1]);
        std::string ext = getExtension(argv[1]);

#ifdef USE_RANDOM_DC3D
				filenameRays = cmd_option("-fr", "", "load text file with rays");
				//std::clog << "[RAYS FILES] : " << filenameRays << std::endl;
#endif
        TriMesh mesh;

       if (ext == "off" || ext == "OFF")
            mesh = loadOFFMesh(file);
        else if (ext == "obj" || ext == "OBJ")
            mesh = loadOBJMesh(file);
        else
            throw "Unknown file type!";
        
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
