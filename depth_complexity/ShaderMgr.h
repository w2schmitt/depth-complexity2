#ifndef SHADER_MGR_H
#define SHADER_MGR_H

#define STRING_BUFFER 2048

#include <GL/glew.h>
#ifdef __APPLE__
#include <GLUT/glut.h>
#else
#include <GL/glut.h>
#include <algorithm>
#endif

#include <iostream>
#include <fstream>

class ShaderMgr{
public:
  ShaderMgr():_shaderStatus(true){};
  GLuint createShaderProgram(const char* vFilename, const char* pFilename);  
  void linkShaderProgram(GLuint shaderProgramId);
  
  void checkShaderStatus(){ 
    if (_shaderStatus)
      std::clog << "[OK] Shaders Ok." << std::endl;
    else
      std::clog << "[ERROR] Verify Shaders" << std::endl;
  }
   
private:
  bool _shaderStatus;
  GLuint createShader(const char* filename, GLuint shaderType);
  char *readShaderFile(const char *filename);
  // Internal function to check if shader was compiled correctly.
  bool checkShadersStatus(GLuint shaderID, GLenum statusType, const char* msg);
  
};

#endif

