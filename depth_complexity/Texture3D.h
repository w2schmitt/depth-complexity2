/* 
 * File:   Texture3D.h
 * Author: wagner
 *
 * Created on 14 de Setembro de 2012, 10:53
 */
#include <iostream>
#include <GL/glew.h>
#include <cstdlib>
#include <stdio.h>
#include <iomanip>
#include "util.h"
#include "CImg.h"
#include "ShaderMgr.h"


using namespace cimg_library;


#ifndef TEXTURE3D_H
#define	TEXTURE3D_H


//#define BYTES_PER_TEXEL 3

class Texture3D {
    
private:
  CImg<unsigned int>                                    _tex3D;
  vec3d                                                 _texSize;
  BoundingBox                                           _aabb;
  GLuint                                                _texID;
  
  GLuint                                                _shaderTex3d;
  bool                                                  _shaderReady;
  unsigned int                                          *_glTex3d;
  
  
  GLuint createOpenglTexture3D(int width, int height, int depth, const unsigned int* texels);
  GLuint createOpenglTexture3D2(int width, int height, int depth, const float* texels);
  
public:
    int _w,_h,_d;
    Texture3D(){}
    void CreateTexture3D(unsigned int _w, unsigned int _h, unsigned int _d, unsigned int _ch, float _init);
    Texture3D(const Texture3D& orig);
    virtual ~Texture3D();
    //GLuint textureID;
    
    void setTex3dSize(vec3d size) { this->_texSize = size;}
    void updateTexture3D(Segment s, unsigned int value);
    void cimg2Tex(unsigned int maxDC);
    void buildGLTexture();
    
    GLuint texture3DId() const       { return _texID; }
    void setMeshBoundingbox(BoundingBox aabb) { this->_aabb = aabb; }
    
    void drawPixel(vec3d pos, unsigned int color);
    void bresenham_line_3D(vec3d &v1, vec3d &v2, unsigned int color);
    
    //shader methods
    void setShaderTex3d();
    
    
    
};

#endif	/* TEXTURE3D_H */

