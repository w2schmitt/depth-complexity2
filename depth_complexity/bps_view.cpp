/* 
 * File:   bps_view.cpp
 * Author: sombra
 *
 * Created on 6 de Agosto de 2012, 16:03
 */

#include <cstdlib>
# include <GL/glew.h>
# include <GL/glfw.h>
# include <AntTweakBar.h>
# ifdef __APPLE__
#   include <GLUT/glut.h>
# else
#   include <GL/glut.h>
# endif

# include "dc_3d_random.h"
# include "dc_3d.h"
# include "camera/Camera.h"
#include "bps_model.h"

using namespace std;

// parametros para o setup do depthcompleity
TwBar *mainBar;
TwBar *treeBar;

struct MeshTree{
  MeshTree *front;
  MeshTree *back;
  bool leaf;
  TriMesh mesh;
  
  //visual configs
  vec4d color;
  bool show;
  
  MeshTree(TriMesh newMesh):front(NULL),back(NULL),mesh(newMesh){}
  ~MeshTree()
  { delete front;
    delete back;
    delete mesh;
  }
  TriMesh makeMesh();
  TriMesh recompute(int maxDepth);
  void putMeshColor();
  vec4d MeshTree::computeLinesRandom(int const threshold);
  vec4d MeshTree::computeLinesNormal(int const threshold);
};

MeshTree tree;

struct _barParam{
  bool doGoodRays = true;//retirar
  bool maxRays = true;//retirar
  bool doHist = true;//retirar
  unsigned int showRayIndex = 0;
  int discretSteps = 2;
  unsigned int maxDepth = 0;
  unsigned int threshold = 10;
  int fboWidth = 512;
  int fboHeight = 512;
  vec3f top(0.25, 0.25, .5), mid(0.75, 0.75, .85), bot(1, 1, 1);
  vec4f objdiff(0.55, 0.5, 0, 0.5), objspec(.75, .75, .75, .2);
  GLfloat shine = 50;
  //vec4f sphereColor(0.0, 1.0, 0.0, 1.0);
  //float sphereRadius = 0.01;
  bool showObj = true;
  bool normal_dc3d = false; // depth complexity normal ou random
  bool showOriginal = true;
  
  void TW_CALL setShowOriginal(const void *Valor, void *clientData);
  void TW_CALL getShowOriginal(void *Valor, void *clientData);
} Config;

bool input_model = true; // se um modelo foi fornecido na chamada do programa
TriMesh originalMesh; // modelo fornecido para ser particionado
TriMesh dividedMesh = NULL;    // mesh after the divisions
std::vector<Triangle> sorted_faces;

DepthComplexity3D *dc3d = NULL;
RDepthComplexity3D *dc3dr = NULL;
//________________________________________________ CAMERA
Camera camera;
//______________ MOUSE
int mDx = 0;
int mDy = 0;
bool mclicked;
int mbutton;
int mwhell;
//______________ WINDOW
static int winWidth, winHeight;

//subroutines
int doInteractive();
void setupDepth3d();
void setupRDepth3d();
void setupTweakBar(TwBar* bar);
void drawMesh(const TriMesh& mesh, const vec3f& dir);
TriMesh recompute(TriMesh& mesh, int threshold);

void GLFWCALL mouse_motion(int x, int y);
void GLFWCALL mouse_click(int button, int action);
void GLFWCALL mouse_whell(int pos);
void GLFWCALL WindowSizeCB(int width, int height);

int main(int argc, char** argv) {
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


//_______________________________________________________________ Called when a mouse button is pressed or released
void GLFWCALL mouse_click(int button, int action){
	if(TwEventMouseButtonGLFW(button,action))
		return;
	glfwGetMousePos(&mDx,&mDy);
	mclicked=action==GLFW_PRESS;
	mbutton=button;
}

//_______________________________________________________________ Called when the mouse whell move
void GLFWCALL mouse_whell(int pos){
  if(TwEventMouseWheelGLFW(pos)){
    mwhell = pos;
    return;
  }
  float Dw = pos - mwhell;
  float delta = 0.1;
  
  camera.MoveFrente(Dw * delta);
  
  mwhell = pos;
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
			case GLFW_MOUSE_BUTTON_LEFT : //look
				//if( glutGetModifiers() == GLUT_ACTIVE_SHIFT){
					camera.lookLefRigObj(dx * delta );
					camera.lookUpDownObj(dy * delta);
				/*} else { 
					camera.lookLefRig(dx * delta);
					camera.lookUpDown(dy * delta);
				}*/
			break;
			case GLFW_MOUSE_BUTTON_RIGHT:			
				//if( glutGetModifiers() != GLUT_ACTIVE_SHIFT){
					camera.MoveLado(dx*delta);
					camera.MoveCimaBaixo(-dy*delta);
				/*}else{
					camera.MoveLadoObj(dx*delta);
					camera.MoveCimaBaixoObj(-dy*delta);				
				//}*/
			break;
			case GLFW_MOUSE_BUTTON_MIDDLE: 
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

// Callback function called by GLFW when window size changes
void GLFWCALL WindowSizeCB(int width, int height)
{
    winWidth = width;
    winHeight = height;
    // Send the new window size to AntTweakBar
    TwWindowSize(width, height);
}

void TW_CALL _barParam::getShowOriginal(void* Valor, void* clientData){
  *(bool*)Valor = showOriginal;
}

void TW_CALL _barParam::setShowOriginal(const void* Valor, void* clientData){
  showOriginal = *(const bool *)Valor;
  
  if(showOriginal){
    sorted_faces = originalMesh.faces;
  }else{
    sorted_faces = dividedMesh.faces;
  }
}
int doInteractive()
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
    std::cerr << "GLEW VERSION: "<<glewGetString(GLEW_VERSION)<<std::endl;
    std::clog << "---------------- \n";
 

    glfwEnable(GLFW_MOUSE_CURSOR);
    glfwEnable(GLFW_KEY_REPEAT);
    glfwSetWindowTitle("Minimal intersection division");

	
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
    glfwSetMouseWheelCallback(mouse_whell);
    glfwSetKeyCallback((GLFWkeyfun)TwEventKeyGLFW);
    glfwSetCharCallback((GLFWcharfun)TwEventCharGLFW);

    mwhell = glfwGetMouseWheel();
    
    mainBar = TwNewBar("Controls");

    glEnable(GL_LIGHTING);
    glEnable(GL_LIGHT0);
    glDisable(GL_CULL_FACE);
    glEnable(GL_DEPTH_TEST);

    //Camera cam;
    BoundingBox aabb = originalMesh.aabb;
    
    camera.bbox(float3(aabb.min.x,aabb.min.y,aabb.min.z), float3(aabb.max.x,aabb.max.y,aabb.max.z), true );
		camera.front();
    
    //cam.target = aabb.center();
    //cam.up = vec3f(0, 1, 0);
    //cam.pos = cam.target + vec3f(0, 0, 2*aabb.extents().z);

    while( glfwGetWindowParam(GLFW_OPENED) && !glfwGetKey(GLFW_KEY_ESC) ) {
        glClear( GL_COLOR_BUFFER_BIT|GL_DEPTH_BUFFER_BIT );

        drawBackground(Config.top, Config.mid, Config.bot);

        setupCamera(camera);
        
        /*camera.update();	
	    camera.lookAt();*/

        /*GLfloat lpos[4] = { camera.GetEye().x, camera.GetEye().y, camera.GetEye().z, 1 };
        glLightfv(GL_LIGHT0, GL_POSITION, lpos);*/
       

        //drawPlaneIntersection(intersectionVectors);
        drawRays();


				
        //glMaterialfv(GL_FRONT_AND_BACK, GL_DIFFUSE, &objdiff.x);
        glMaterialfv(GL_FRONT_AND_BACK, GL_SPECULAR, &Config.objspec.x);
        glMaterialf(GL_FRONT_AND_BACK, GL_SHININESS, Config.shine);
        
        if (Config.showObj){
          if(Config.showOriginal){
            drawMesh(originalMesh, vec3f(camera.GetDir().x,camera.GetDir().y,camera.GetDir().z));
          }else{
            drawMesh(dividedMesh, vec3f(camera.GetDir().x,camera.GetDir().y,camera.GetDir().z));
          }
        } 
          
        
        /*
        Triangle *t = &mesh.faces[550];
        t->ca = vec4d(1.0f, 0.0f, 0.0f, 0.7f);
        t->cb = vec4d(1.0f, 0.0f, 0.0f, 0.7f);
        t->cc = vec4d(1.0f, 0.0f, 0.0f, 0.7f);

        t = &mesh.faces[608];
        t->ca = vec4d(1.0f, 0.0f, 0.0f, 0.7f);
        t->cb = vec4d(1.0f, 0.0f, 0.0f, 0.7f);
        t->cc = vec4d(1.0f, 0.0f, 0.0f, 0.7f);
         */
        // Draw tweak bars
        TwDraw();

        // Present frame buffer
        glfwSwapBuffers();
    }

    return 0;
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

void setupDepth3d(){
  
  if(dc3d == NULL)
    dc3d = new DepthComplexity3D(Config.fboWidth, Config.fboHeight, Config.discretSteps);
  
  dc3d->setComputeMaximumRays(Config.doGoodRays);
  dc3d->setComputeHistogram(Config.doHist);
  dc3d->setThreshold(Config.threshold);
}

void setupRDepth3d(){
  
  if(dc3dr == NULL)
    dc3dr = new RDepthComplexity3D(Config.fboWidth, Config.fboHeight, Config.discretSteps);  
  
  dc3dr->setComputeMaximumRays(Config.doGoodRays);
  dc3dr->setComputeHistogram(Config.doHist);
  dc3dr->setThreshold(Config.threshold);
}

void setupTweakBar(TwBar* bar){
  TwDefine(" GLOBAL ");

  //TwAddVarRW(bar, "showPlanes", TW_TYPE_BOOLCPP, &showPlanes, " label='show discret. planes' ");

  //TwAddVarRW(bar, "goodRays", TW_TYPE_BOOLCPP, &doGoodRays, " label='show more rays' ");

  TwAddVarRW(bar, "goodThreshold", TW_TYPE_UINT32, &Config.threshold, " label='intersection threshold' min=0 ");

  TwAddVarRW(bar, "discretSteps", TW_TYPE_UINT32, &Config.discretSteps, " label='discret. steps' min=2 ");
  TwAddVarRW(bar, "Use normal dc3d", TW_TYPE_BOOLCPP, &Config.normal_dc3d, NULL);
  TwAddButton(bar, "recompute", recompute, (void*)&originalMesh, " label='Recompute' ");

  TwAddVarRO(bar, "maxDepth", TW_TYPE_UINT32, &Config.maxDepth, " label='Max. depth' ");

  TwAddVarRW(bar, "top", TW_TYPE_COLOR3F, &Config.top.x, " group='background' ");
  TwAddVarRW(bar, "mid", TW_TYPE_COLOR3F, &Config.mid.x, " group='background' ");
  TwAddVarRW(bar, "bot", TW_TYPE_COLOR3F, &Config.bot.x, " group='background' ");

  TwAddVarRW(bar, "Diff", TW_TYPE_COLOR4F, &Config.objdiff.x, " group='Object' ");
  TwAddVarRW(bar, "Spec", TW_TYPE_COLOR4F, &Config.objspec.x, " group='Object' ");
  TwAddVarRW(bar, "Shin", TW_TYPE_FLOAT, &Config.shine, " group='Object' min='1' max='128' ");
  TwAddVarRW(bar, "Show", TW_TYPE_BOOLCPP, &Config.showObj, " group='Object' ");
  TwAddVarRW(bar, "Show Original", TW_TYPE_BOOLCPP, &Config.showOriginal, " group='Object' ");

  //TwAddVarRW(bar, "radius", TW_TYPE_FLOAT, &Config.sphereRadius, "  group='Sphere' label='radius' min=0.0 step=0.001 max=2.0");
  //TwAddVarRW(bar, "Color", TW_TYPE_COLOR4F, &Config.sphereColor.x, " group='Sphere' ");
  //TwAddVarRW(bar, "Show Ray", TW_TYPE_UINT32, &showRayIndex, " group='rays' min=0");
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

void drawMesh(const TriMesh& mesh, const vec3f& dir)
{
    
    glColorMaterial(GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE);
    glEnable(GL_COLOR_MATERIAL);
    
    //std::clog << "sorting...";
    if(sorted_faces.empty()) {
      sorted_faces = mesh.faces;
    }
    
    std::sort(sorted_faces.begin(), sorted_faces.end(), ByDist(dir));
    //std::clog << "done" << std::endl;
    
    //glEnable(GL_BLEND);
    //glBlendFunc(GL_ONE, GL_ONE_MINUS_SRC_ALPHA);
    //glDisable(GL_DEPTH_TEST);
    glEnable(GL_DEPTH_TEST);
    glEnable(GL_LIGHTING);
    glLightModeli(GL_LIGHT_MODEL_TWO_SIDE, 1);    

    glEnable(GL_VERTEX_ARRAY);

    glEnableClientState(GL_VERTEX_ARRAY);    
    glVertexPointer(3, GL_DOUBLE, 2*sizeof(vec3d)+sizeof(vec4d), &sorted_faces[0].a.x);
    
    glEnableClientState(GL_COLOR_ARRAY);
    glColorPointer(4, GL_DOUBLE, 2*sizeof(vec3d)+sizeof(vec4d), &sorted_faces[0].ca.x);

    glEnableClientState(GL_NORMAL_ARRAY);
    glNormalPointer(GL_DOUBLE, 2*sizeof(vec3d)+sizeof(vec4d), &sorted_faces[0].na.x);

    glDrawArrays(GL_TRIANGLES, 0, sorted_faces.size()*3);

    glDisableClientState(GL_VERTEX_ARRAY);
    glDisable(GL_VERTEX_ARRAY);
    glDisableClientState(GL_NORMAL_ARRAY);
    glDisable(GL_NORMAL_ARRAY);
    glDisableClientState(GL_COLOR_ARRAY);
    glDisable(GL_COLOR_ARRAY);

    glDisable(GL_BLEND);
    glEnable(GL_DEPTH_TEST);
    glDisable(GL_COLOR_MATERIAL);
    
}
//TODO test
TriMesh recompute(TriMesh& mesh, int maxDepth){
  
  tree = new MeshTree(mesh);
  
  TriMesh meshReturn = tree.recompute(maxDepth);
  
  
}

//TODO test
TriMesh MeshTree::recompute(int maxDepth){
  vec4d cuttingPlane;
  if(Config.normal_dc3d){
    cuttingPlane = computeLinesNormal(maxDepth);
  }else{
    cuttingPlane = computeLinesRandom(maxDepth);
  }
  
  if(cuttingPlane.x == 0.0 && cuttingPlane.y == 0.0 && 
     cuttingPlane.z == 0.0 ){
    putMeshColor();
    front = NULL;
    back = NULL;
    return mesh;    
  }
  
  TriMesh meshFront, meshBack;
  subDivideModel(cuttingPlane, &mesh, &meshFront, &meshBack);
  delete mesh;
  mesh = NULL;
  
  front = new MeshTree(meshFront);
  back = new MeshTree(meshBack);
  
  TriMesh frontMesh = front->recompute(maxDepth);
  TriMesh backMesh = back->recompute(maxDepth);
  
  TriMesh meshReturn;
  meshReturn.faces.insert(meshReturn.faces.end(),
                          frontMesh.faces.begin(),
                          frontMesh.faces.end());
  meshReturn.faces.insert(meshReturn.faces.end(),
                          backMesh.faces.begin(),
                          backMesh.faces.end());
  
  meshReturn.aabb = frontMesh.aabb;
  meshReturn.aabb.merge(backMesh.aabb.max);
  meshReturn.aabb.merge(backMesh.aabb.min);
  
  return meshReturn;
}

//TODO teste
void MeshTree::putMeshColor(){
  
  double red = static_cast<double>(rand())/static_cast<double>(RAND_MAX),
         blue = static_cast<double>(rand())/static_cast<double>(RAND_MAX),
         green = static_cast<double>(rand())/static_cast<double>(RAND_MAX);
  this->color = vec4d(red, blue, green, 1.0f);
  
  std::vector<Triangle>::iterator meshIterator = this->mesh.faces.begin();
  std::vector<Triangle>::iterator meshEnd = this->mesh.faces.end();
  
  for(; meshIterator != meshEnd; ++meshIterator){
    meshIterator->ca = this->color;
    meshIterator->cb = this->color;
    meshIterator->cc = this->color;
  }
}

//TODO test
vec4d MeshTree::computeLinesNormal(int const maxDepth){
  dc3d->setComputeGoodRays(Config.doGoodRays);
  dc3d->setThreshold(Config.threshold);
  dc3d->setComputeHistogram(Config.doHist);
  dc3d->setComputeMaximumRays(Config.maxRays);
  dc3d->setDiscretSteps(Config.discretSteps);
  dc3d->setMaximum(Config.maxDepth);
  dc3d->process(&mesh);
  
  if(dc3dr->maximum() <= maxDepth){
    return vec4d(0.0, 0.0, 0.0, 0.0);
  }
  
  return defineCuttingPlane(vector<Segment>(dc3d->maximumRays().begin(), 
                                            dc3d->maximumRays().end()),
                                            mesh.aabb);  
}

//TODO test
vec4d MeshTree::computeLinesRandom(int const maxDepth){
  dc3dr->setComputeGoodRays(Config.doGoodRays);
  dc3dr->setThreshold(Config.threshold);
  dc3dr->setComputeHistogram(Config.doHist);
  dc3dr->setComputeMaximumRays(Config.maxRays);
  dc3dr->setDiscretSteps(Config.discretSteps);
  dc3dr->setMaximum(Config.maxDepth);
  dc3dr->process(&mesh);
  
  if(dc3dr->maximum() <= maxDepth){
    return vec4d(0.0, 0.0, 0.0, 0.0);
  }
  
  return defineCuttingPlane(vector<Segment>(dc3dr->maximumRays().begin(), 
                                            dc3dr->maximumRays().end()),
                                            mesh.aabb);
  
}