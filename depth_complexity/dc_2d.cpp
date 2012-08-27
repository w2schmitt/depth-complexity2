#include "dc_2d.h"

#include <iostream>
#include <cassert>

// Internal function to check if the framebuffer was created correctly.
bool checkFramebufferStatus();

// Draw polygons in framebuffer.
//void drawQuad(const Point &p1, const Point &p2, const Point &p3, const Point &p4);
//void drawTriangle(const Point &p1, const Point &p2, const Point &p3);
void drawPolygon(const std::vector<Point> &polygon);

bool shrinkSegment(Point &p1, Point &p2, vec3d d1, vec3d d2);

#define WIDTH  512
#define HEIGHT 512
#define DEPTH  512
//#define BYTES_PER_TEXEL 3

DepthComplexity2D::DepthComplexity2D(const int fboWidth, const int fboHeight){
  _fboWidth = fboWidth;
  _fboHeight = fboHeight;
  _computeHistogram = false;
  _computeMaximumRays = false;
  _computeGoodRays = false;
  _maximum = 0;
  _threshold = 0;
  _shaderclearBuffer = 0;
  _shaderCountDC = 0;
  _status = false;  
  
  CImg casa;
  
  // initialize 3d texture (size x, size y, size z, channels, px values)
  // --> channel 1 -> DC of the ray
  // --> chanell 2 -> ray counter
  _tex3D = CImg<unsigned int>(WIDTH,HEIGHT,DEPTH,2,0);
  
  assert(initFBO());
  
  //set up shaders
  ShaderMgr shaderMgr; 
  _shaderCountDC = shaderMgr.createShaderProgram("shader/dc.vert", "shader/dc.frag");        
  _shaderclearBuffer = shaderMgr.createShaderProgram("shader/clear.vert", "shader/clear.frag");
  
  shaderMgr.linkShaderProgram(_shaderCountDC); 
  shaderMgr.linkShaderProgram(_shaderclearBuffer);
    
  shaderMgr.checkShaderStatus();
}


DepthComplexity2D::~DepthComplexity2D(){

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
   
  findDepthComplexity2D();
  
  if (_computeHistogram or _computeMaximumRays or _computeGoodRays)
    findMaximumRaysAndHistogram();
}


void DepthComplexity2D::findDepthComplexity2D() {
  glBindFramebufferEXT(GL_FRAMEBUFFER_EXT, _fboId);
  
  glPushAttrib(GL_VIEWPORT_BIT | GL_ENABLE_BIT);
    glViewport(0, 0, _fboWidth, _fboHeight);

    // desable all anti-aliasing to avoid precision error
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
		
    // set shader to clear the counter buffer ---------
    setShaderClearCounterBuffer();
    
    // Ensure that all texture writing is done
    glMemoryBarrierEXT(GL_FRAMEBUFFER_BARRIER_BIT_EXT);

    // Set shader to write DC in the counter buffer 
    setShaderCountDC();
    
    for (unsigned i=0; i<_segments->size(); ++i) {
        if (!_segments->at(i).active) continue;

        Segment Tv0 = computeDualSegmentFromPoint(_segments->at(i).a);
        Segment Tv1 = computeDualSegmentFromPoint(_segments->at(i).b);

	if (Tv0.a.y < Tv1.a.y) std::swap(Tv0, Tv1); 
            double t1,t2;
      
        if (segmentIntersection2D(Tv0,Tv1, &t1, &t2)){
            std::vector<Point> polygonA;
            std::vector<Point> polygonB;
            if (lineIntersection2D(Tv0, Tv1, &t1, &t2)){
                Point inter = Tv0.a*(1.0f-t1) + t1*Tv0.b;
              
                clipPolygon(Tv0.a, inter, Tv1.a, polygonA);
                clipPolygon(Tv0.b, inter, Tv1.b, polygonB);          
                
                drawPolygon(polygonA);
                drawPolygon(polygonB);
             }                 
        }
        else {					
								
            std::vector<Point> polygonA;
            clipPolygon(Tv0.a, Tv0.b, Tv1.b, Tv1.a, polygonA);
		
            if (polygonA.size() == 0) continue;
            
            drawPolygon(polygonA);
        }
    }
    
    // Ensure that all texture writing is done
    glMemoryBarrierEXT(GL_TEXTURE_UPDATE_BARRIER_BIT_EXT);
    // Disable Shaders
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

void DepthComplexity2D::clipPolygon(const Point &p1, const Point &p2, const Point &p3, const Point &p4, std::vector<Point> &polygon){
		
    Segment upSeg(Point(0,1), Point(1,1));
    Segment lowerSeg(Point(0,0), Point(1,0));
    double t1,t2,t3,t4;
       
    if ((p1.y < p4.y) || (p2.y < p3.y)) return;
    Point np1, np2, np3, np4;
    
    if (p1.y > 1.0f){
      
      lineIntersection2D(upSeg, Segment(p1, p2), &t1, &t2); 
      if (p4.y > 1.0f){
        lineIntersection2D(upSeg, Segment(p4, p3), &t3, &t4);    
        np1 = Point(t1, 1.0f); 
				np4 = Point(t3, 1.0f);

        bool ok = shrinkSegment( np1, np4, vec3d::left(), vec3d::right() );
        polygon.push_back(np1);
        if (ok)	polygon.push_back(np4);      
      }
      
      else {
        np1 = Point(t1, 1.0f); 
	np2 = Point(0.0f, 1.0f);
	np4 = p4;
        
        bool ok1 = shrinkSegment( np4, np2, vec3d::up()  , vec3d::zero() ),
             ok2 = shrinkSegment( np1, np2, vec3d::left(), vec3d::zero() );
						
        if (ok2)	polygon.push_back(np1);						
        polygon.push_back(np2);						
        if (ok1 && ok2) polygon.push_back(np4);        
						}
      
    }
    else {      
      np1 = p1; 
      np4 = p4;
      
      bool ok = shrinkSegment( np1, np4, vec3d::down(), vec3d::up() );      
      polygon.push_back(np1);
      if (ok) polygon.push_back(np4); 
    }
    
    
    if (p3.y < 0.0f){
      
      lineIntersection2D(lowerSeg, Segment(p4, p3), &t1, &t2); 
      if (p2.y < 0.0f){
        lineIntersection2D(lowerSeg, Segment(p1, p2), &t3, &t4);
        np3 = Point(t1, 0.0f);
        np2 = Point(t3, 0.0f);         
        bool ok = shrinkSegment( np2, np3, vec3d::left(), vec3d::right() );

        polygon.push_back(np3);
        if (ok) polygon.push_back(np2);		
      }
      else{
        np2 = p2;
        np4 = Point(1.0f, 0.0f);
        np3 = Point(t1, 0.0f); ;

        bool ok1 = shrinkSegment( np2, np4, vec3d::down()  , vec3d::zero() ),
             ok2 = shrinkSegment( np3, np4, vec3d::right(), vec3d::zero() );
        
        if (ok2) polygon.push_back(np3);						
        polygon.push_back(np4);
        if (ok1 && ok2) polygon.push_back(np2);
      }
      
    }
    else{      
      np2 = p2; 
      np3 = p3;

      bool ok = shrinkSegment( np2, np3, vec3d::down(), vec3d::up() );
      polygon.push_back(np3);
      if (ok) polygon.push_back(np2);
    }
}


void DepthComplexity2D::clipPolygon(const Point &p1, const Point &p2, const Point &p3, std::vector<Point> &polygon){
		
		Segment upSeg(Point(0,1), Point(1,1));
		Segment lowerSeg(Point(0,0), Point(1,0));
		double t1,t2,t3,t4;

		// discard polygon
		if ( (p1.y > 1.0f && p2.y > 1.0f && p3.y > 1.0f) or
				 (p1.y < 0.0f && p2.y < 0.0f && p3.y < 0.0f))
		{
				return;
		}
		Point np1, np2, np3;

		if (p1.y > 1.0f){ // =================================================
				
				lineIntersection2D(upSeg, Segment(p1, p2), &t1, &t2); 
				
				if (p3.y > 1.0f){
						lineIntersection2D(upSeg, Segment(p3, p2), &t3, &t4);
						np1 = Point(t1, 1.0f); 
						np3 = Point(t3, 1.0f);
						bool ok = shrinkSegment( np1, np3, vec3d::left(), vec3d::right() );

						polygon.push_back(np1);
						if (ok)	polygon.push_back(np3);

        } 
				else{
						np1 = Point(t1, 1.0f); 
						np2 = Point(0.0f, 1.0f);
						np3 = p3;
	
						bool ok1 = shrinkSegment( np3, np2, vec3d::up()  , vec3d::zero() ),
                 ok2 = shrinkSegment( np1, np2, vec3d::left(), vec3d::zero() );
						
						if (ok1)	polygon.push_back(np3);						
						polygon.push_back(np2);						
						if (ok1 && ok2) polygon.push_back(np1);
				}
		} // =================================================================
		else if (p1.y > 0.0f){

				if (p3.y > 0.0f){
						np1 = p1; 
						np3 = p3;
						bool ok;

						if (p1.y > p3.y)
								ok = shrinkSegment( np1, np3, vec3d::down(), vec3d::up() );
						else
								ok = shrinkSegment( np1, np3, vec3d::up(), vec3d::down() );

						polygon.push_back(np1);
						if (ok) polygon.push_back(np3);
				}
		} // =================================================================
		else {

				lineIntersection2D(lowerSeg, Segment(p1, p2), &t1, &t2); 

				if (p3.y > 0.0f){
						np1 = Point(t1, 0.0f); 
						np2 = Point(1.0f, 0.0f);
						np3 = p3;

						bool ok1 = shrinkSegment( np3, np2, vec3d::down()  , vec3d::zero() ),
											ok2 = shrinkSegment( np1, np2, vec3d::right(), vec3d::zero() );
						
						if (ok1) polygon.push_back(np3);						
						polygon.push_back(np2);
						if (ok1 && ok2) polygon.push_back(np1);
				}
				if (p3.y < 0.0f){
						lineIntersection2D(lowerSeg, Segment(p3, p2), &t3, &t4);
						np1 = Point(t3, 0.0f); 
						np3 = Point(t1, 0.0f);
						bool ok = shrinkSegment( np1, np3, vec3d::left(), vec3d::right() );
	
						polygon.push_back(np1);
						if (ok) polygon.push_back(np3);		
            }
		}
		Point middle = (p1+p3)*0.5;
		vec3d dir = middle - p2; dir.normalize();
		np2 = p2;
		shrinkSegment( np2, middle, dir, vec3d::zero() );

		polygon.push_back(np2);
}

bool shrinkSegment(Point &p1, Point &p2, vec3d d1, vec3d d2){
			
    double pixelSize = 1.0/512.0;
    double step = 0.5*pixelSize*sqrt(2);

    bool compare_x = p1.x < p2.x;
    bool compare_y = p1.y < p2.y;

    Point cp1, cp2;

    // Shrink segment
    cp1 = p1 + step*d1;
    cp2 = p2 + step*d2;

    // Check if it is not crossed
    if ( (compare_x != (cp1.x < cp2.x)) || (compare_y != (cp1.y < cp2.y)) ){
        Point middle = (p1 + p2)*0.5;
        p1 = p2 = middle;
        return false;
    }

    p1 = cp1;
    p2 = cp2;

    return true;
}

Segment DepthComplexity2D::computeDualSegmentFromPoint(const Point &p) {
    Segment seg;
    double t1, t2;

    // Finding left point
    bool r1 = lineIntersection3D(_to, Segment(_from.a, p), &t1, &t2);
                //printf("Point P(%f, %f) - at1 = %f\n",p.x, p.y, t1);
    assert(r1);
    seg.a = Point(0.0, t1);

    // Finding right point
    bool r2 = lineIntersection3D(_to, Segment(_from.b, p), &t1, &t2);
    //printf("Point P(%f, %f) - at1 = %f\n",p.x, p.y, t1);
    assert(r2);
    seg.b = Point(1.0, t1);

    return seg;
}

void drawPolygon(const std::vector<Point> &polygon){
    glBegin(GL_POLYGON);
        for (unsigned i=0; i< polygon.size(); ++i){
                glVertex2f(polygon[i].x,polygon[i].y);
        }
    glEnd();
}


/*
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
*/

unsigned int DepthComplexity2D::findMaxValueInCounterBuffer() {
  const int pixelNumber = _fboWidth * _fboHeight;
  unsigned int colorBuffer[pixelNumber];
  
  glBindTexture(GL_TEXTURE_2D, _counterBuffId);
  glGetTexImage( GL_TEXTURE_2D, 0 , GL_RED_INTEGER, GL_UNSIGNED_INT, colorBuffer ); 
  glBindTexture(GL_TEXTURE_2D, 0);
  
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
          double t1 = ((double)c+0.5)/(double)_fboWidth;
          double t2 = ((double)r+0.5)/(double)_fboHeight;
          seg.a = _from.a*(1.f-t1) + _from.b*t1;
          seg.b = _to.a*(1.f-t2) + _to.b*t2;          
          seg.sortPoints();
          
          updateTexture3D(seg, val);
          
          
          /*
          if (val == _maximum){
            //if (_maximumRays.size() < 1)
              _maximumRays.insert(seg);
          }
          else if (val >= _threshold){			   
            //if (_goodRays[val].size() < 1)
		_goodRays[val].insert(seg);            
          }
          */
        }
     }
  }
  glPopAttrib();
  glBindFramebufferEXT(GL_FRAMEBUFFER_EXT, 0);
}



const unsigned int colorTableSize=6;

float CTable[colorTableSize][4] = {
    {0.0, 0.0, 1.0, 0.85}, // Azul
    {1.0, 0.0, 1.0, 1.0}, // Roxo/Rosa
    {1.0, 0.5, 0.0, 1.0}, // Laranja
    {1.0, 1.0, 0.0, 1.0}, // Amarelo
    {0.0, 1.0, 0.0, 1.0}, // Verde
    {0.0, 1.0, 0.0, 1.0} // Verde
};

/* Linear Interpolation between color1 and color2*/
void lerpColor(float* newColor, const float *color1, const float *color2, float value){
    
    if (value<=0){
        //memcpy(newColor, color1, sizeof(float*)*3);
        newColor[0] = color1[0];
        newColor[1] = color1[1];
        newColor[2] = color1[2];
        newColor[3] = color1[3];
    }
    else if (value>=1){
        newColor[0] = color2[0];
        newColor[1] = color2[1];
        newColor[2] = color2[2];
        newColor[3] = color2[3];
        //std::cout << value << " - [" << color2[0] << " " << color2[1] << " " << color2[2] << "]\n";
        //memcpy(newColor, color2, sizeof(float*)*3);
    }
    else {
        newColor[0] = (1.0-value)*color1[0] + (value)*color2[0];
        newColor[1] = (1.0-value)*color1[1] + (value)*color2[1];
        newColor[2] = (1.0-value)*color1[2] + (value)*color2[2];
        newColor[3] = (1.0-value)*color1[3] + (value)*color2[3];
    }
}


void findColor(float *pxColor, float normalizedDC){
    
    float intensity;
    float tableIndex = normalizedDC*(float)(colorTableSize-1);
    
    int firstColor = (int)tableIndex;
    
    if (firstColor==colorTableSize-1){
        intensity=1.0;
        firstColor--;
    }
    else{
        intensity = tableIndex - (float)firstColor;
    }
    
    int secondColor = ((int)firstColor)+1;
    lerpColor(pxColor, CTable[firstColor], CTable[secondColor], intensity);
}


void DepthComplexity2D::drawPixel(vec3d pos, unsigned int color){
    
    unsigned int w = _tex3D.width();
    unsigned int h = _tex3D.height();
    unsigned int d = _tex3D.depth();
    
    if ( (pos.x > 0 && pos.x < w) &&
         (pos.y > 0 && pos.y < h) &&
         (pos.z > 0 && pos.z < d)){

         _tex3D(pos.x, pos.y, pos.z, 0) += color;
         _tex3D(pos.x, pos.y, pos.z, 1) += 1; // ray counter

    }

}



/*
 * line3d was dervied from DigitalLine.c published as "Digital Line Drawing"
 * by Paul Heckbert from "Graphics Gems", Academic Press, 1990
 * 
 * 3D modifications by Bob Pendleton. The original source code was in the public
 * domain, the author of the 3D version places his modifications in the
 * public domain as well.
 * 
 * line3d uses Bresenham's algorithm to generate the 3 dimensional points on a
 * line from (x1, y1, z1) to (x2, y2, z2)
 * 
 */

/* find maximum of a and b */
#define MAX(a,b) (((a)>(b))?(a):(b))

/* absolute value of a */
#define ABS(a) (((a)<0) ? -(a) : (a))

/* take sign of a, either -1, 0, or 1 */
#define ZSGN(a) (((a)<0) ? -1 : (a)>0 ? 1 : 0)


void DepthComplexity2D::bresenham_line_3D(vec3d &v1, vec3d &v2, unsigned int color){
    int x1, y1, x2, y2, z1, z2;

    x1 = v1.x; y1 = v1.y; z1 = v1.z;
    x2 = v2.x; y2 = v2.y; z2 = v2.z;
    
    int xd, yd, zd;
    int x, y, z;
    int ax, ay, az;
    int sx, sy, sz;
    int dx, dy, dz;

    dx = x2 - x1;
    dy = y2 - y1;
    dz = z2 - z1;

    ax = ABS(dx) << 1;
    ay = ABS(dy) << 1;
    az = ABS(dz) << 1;

    sx = ZSGN(dx);
    sy = ZSGN(dy);
    sz = ZSGN(dz);

    x = x1;
    y = y1;
    z = z1;

    if (ax >= MAX(ay, az))            /* x dominant */
    {
        yd = ay - (ax >> 1);
        zd = az - (ax >> 1);
        for (;;)
        {
            drawPixel(vec3d(x, y, z), color);
            if (x == x2)
            {
                return;
            }

            if (yd >= 0)
            {
                y += sy;
                yd -= ax;
            }

            if (zd >= 0)
            {
                z += sz;
                zd -= ax;
            }

            x += sx;
            yd += ay;
            zd += az;
        }
    }
    else if (ay >= MAX(ax, az))            /* y dominant */
    {
        xd = ax - (ay >> 1);
        zd = az - (ay >> 1);
        for (;;)
        {
            drawPixel(vec3d(x, y, z),color);
            if (y == y2)
            {
                return;
            }

            if (xd >= 0)
            {
                x += sx;
                xd -= ay;
            }

            if (zd >= 0)
            {
                z += sz;
                zd -= ay;
            }

            y += sy;
            xd += ax;
            zd += az;
        }
    }
    else if (az >= MAX(ax, ay))            /* z dominant */
    {
        xd = ax - (az >> 1);
        yd = ay - (az >> 1);
        for (;;)
        {
            drawPixel(vec3d(x, y, z),color);
            if (z == z2)
            {
                return;
            }

            if (xd >= 0)
            {
                x += sx;
                xd -= az;
            }

            if (yd >= 0)
            {
                y += sy;
                yd -= az;
            }

            z += sz;
            xd += ax;
            yd += ay;
        }
    }
}



GLuint createTexture3D(int width, int height, int depth, const float* texels){
    //allocate texture name
    GLuint id;
    glGenTextures(1, &id);
    glBindTexture(GL_TEXTURE_3D, id);
    
    glTexEnvf(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_MODULATE);
    glTexParameterf(GL_TEXTURE_3D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
    glTexParameterf(GL_TEXTURE_3D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
    
    glTexParameterf(GL_TEXTURE_3D, GL_TEXTURE_WRAP_S, GL_REPEAT);
    glTexParameterf(GL_TEXTURE_3D, GL_TEXTURE_WRAP_T, GL_REPEAT);
    glTexParameterf(GL_TEXTURE_3D, GL_TEXTURE_WRAP_R, GL_REPEAT);
    
    glTexImage3D(GL_TEXTURE_3D, 0, GL_RGBA, width, height, depth, 0, GL_RGBA, GL_FLOAT, texels);
    
    glBindTexture(GL_TEXTURE_3D, 0);
    
    return id;
    
}



#define LAYER(r, c)	(WIDTH * HEIGHT * r * c)
// 2->1 dimension mapping function
#define TEXEL2(s, t, c)	(c * (s * WIDTH + t))
// 3->1 dimension mapping function
#define TEXEL3(s, t, r, c)	(TEXEL2(s, t, c) + LAYER(r, c))


void DepthComplexity2D::cimg2Tex(unsigned int maxDC){
    // each layer
    unsigned int w = _tex3D.width();
    unsigned int h = _tex3D.height();
    unsigned int d = _tex3D.depth();
    
    unsigned int texsize = w*h*d;    
    float *interlaced_data = new float[4*texsize];
    
    for (unsigned int i=0; i<texsize*4; i++){
        interlaced_data[i] = 0.0f;
    }
    
    for (unsigned int r=0; r<d; r++){
        for (unsigned int s=0; s<w; s++){
            for (unsigned int t=0; t<h; t++){
                float pxColor[4] = {0.0, 0.0, 0.0, 0.0};
                float Nray = _tex3D(s,t,r,1);
                Nray = (Nray==0)? 1 : Nray; //if (Nray==0) Nray=1;
                findColor(pxColor, (_tex3D(s,t,r,0)/Nray)/(float)maxDC);
                
                interlaced_data[TEXEL3(t,s,r,4)+0] = pxColor[0];
                interlaced_data[TEXEL3(t,s,r,4)+1] = pxColor[1];
                interlaced_data[TEXEL3(t,s,r,4)+2] = pxColor[2];
                interlaced_data[TEXEL3(s,t,r,4)+3] = pxColor[3];
            }
        }
    }

    std::cout << "Creating 3D texture..." << std::endl;
    _texID = createTexture3D(_tex3D.width(), _tex3D.height(), _tex3D.depth(), interlaced_data);
    delete[] interlaced_data;

}


void DepthComplexity2D::updateTexture3D(Segment line, unsigned int dc){
     
    unsigned int w = _tex3D.width();
    unsigned int h = _tex3D.height();
    unsigned int d = _tex3D.depth();
    
    vec3d t = _aabb.min;
    vec3d size = _aabb.extents();     
    
    vec3d start = line.a;
    start -= t;
    start.x = floor( (start.x/size.x)*(w-1) + 0.5);
    start.y = floor( (start.y/size.y)*(h-1) + 0.5);
    start.z = floor( (start.z/size.z)*(d-1) + 0.5);
    
    vec3d end = line.b;
    end -= t;
    end.x = floor( (end.x/size.x)*(w-1) + 0.5);
    end.y = floor( (end.y/size.y)*(h-1) + 0.5);
    end.z = floor( (end.z/size.z)*(d-1) + 0.5);
    
    // rasterize ray into a 3D texture, using dc as color
    bresenham_line_3D(start, end, dc);  
}






//void DepthComplexity2D::computeIntersectionPoints(std::list<Point> &pts){
    
    
    
//}

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



void DepthComplexity2D::setShaderClearCounterBuffer(){
    glUseProgram(_shaderclearBuffer);
    
    // Pass MODELVIEW and PROJECTION matrices to Shader
    GLfloat model[16],projection[16];
    glGetFloatv(GL_MODELVIEW_MATRIX, model);
    glGetFloatv(GL_PROJECTION_MATRIX, projection);
    glProgramUniformMatrix4fv(_shaderclearBuffer, glGetUniformLocation(_shaderclearBuffer, "modelViewMat"), 1, GL_FALSE, model);
    glProgramUniformMatrix4fv(_shaderclearBuffer, glGetUniformLocation(_shaderclearBuffer, "projectionMat"), 1, GL_FALSE, projection);

    // Pass Counter Buffer Texture
    glProgramUniform1iEXT(_shaderclearBuffer, glGetUniformLocation(_shaderclearBuffer, "counterBuff"), 0);
      
    glBegin(GL_QUADS);
      glVertex2f(-1.0f, 2.0f);
      glVertex2f(2.0f, 2.0f);
      glVertex2f(2.0f, -1.0f);
      glVertex2f(-1.0f, -1.0f);
    glEnd();
}

void DepthComplexity2D::setShaderCountDC(){
    // Set shader to count DC
    glUseProgram(_shaderCountDC);   
	
    // Pass MODELVIEW and PROJECTION matrices to Shader
    GLfloat model[16],projection[16];
    glGetFloatv(GL_MODELVIEW_MATRIX, model);
    glGetFloatv(GL_PROJECTION_MATRIX, projection);
    glProgramUniformMatrix4fv(_shaderCountDC, glGetUniformLocation(_shaderCountDC, "modelViewMat"), 1, GL_FALSE, model);
    glProgramUniformMatrix4fv(_shaderCountDC, glGetUniformLocation(_shaderCountDC, "projectionMat"), 1, GL_FALSE, projection);

    // Pass counter buff texture
    glProgramUniform1iEXT(_shaderCountDC, glGetUniformLocation(_shaderCountDC, "counterBuff"), 0);
}
