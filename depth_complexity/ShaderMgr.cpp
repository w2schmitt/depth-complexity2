
#include "ShaderMgr.h"


GLuint ShaderMgr::createShaderProgram(const char* vFilename, const char* pFilename){
	
  GLuint shaderProgram = glCreateProgram();
  
  GLuint vShaderId = createShader(vFilename, GL_VERTEX_SHADER);
  GLuint pShaderId = createShader(pFilename, GL_FRAGMENT_SHADER);
  
  glAttachShader(shaderProgram, vShaderId);
  glAttachShader(shaderProgram, pShaderId);
  
  return shaderProgram;
}


GLuint ShaderMgr::createShader(const char* filename, GLuint shaderType){
  
  GLuint shaderID = glCreateShader(shaderType);
  
  // Read shader files
	char *shaderData = readShaderFile(filename);
	
	// Pointer to shader data
	const char *ptVData = shaderData;
  
  glShaderSource(shaderID,  1, &ptVData, NULL);
  glCompileShader(shaderID);
  
  if (shaderType==GL_VERTEX_SHADER)
    checkShadersStatus(shaderID, GL_COMPILE_STATUS, "Vertex Shader");
  else if (shaderType==GL_FRAGMENT_SHADER)
    checkShadersStatus(shaderID, GL_COMPILE_STATUS, "Fragment Shader");
    
	free(shaderData);
  
  return shaderID;
}


void ShaderMgr::linkShaderProgram(GLuint shaderProgramId){
  glLinkProgram(shaderProgramId);
  checkShadersStatus(shaderProgramId, GL_LINK_STATUS, "Link Error" );

  glValidateProgram(shaderProgramId);
  checkShadersStatus(shaderProgramId, GL_VALIDATE_STATUS, "Validation Error" );  
}


char* ShaderMgr::readShaderFile(const char *filename){
	
	FILE *fp;
	char *content = NULL;

	int count=0;

	if (filename != NULL) {
		fp = fopen(filename,"rt");
		if (fp != NULL) {      
			  fseek(fp, 0, SEEK_END);
			  count = ftell(fp);
			  rewind(fp);

				if (count > 0) {
					content = (char *)malloc(sizeof(char) * (count+1));
					count = fread(content,sizeof(char),count,fp);
					content[count] = '\0';
			}
			fclose(fp);
		} else {
			std::cerr << "  --> [ERROR] Cannot load Vertex or Fragment Shader File! " << std::endl;
		}
	}
	
	return content;
}


bool ShaderMgr::checkShadersStatus(GLuint shaderID, GLenum statusType, const char* msg){
	
	GLint ok=0;
  
  if (statusType==GL_COMPILE_STATUS)
    glGetShaderiv(shaderID, statusType, &ok);
  else
    glGetProgramiv(shaderID, statusType, &ok);
	
	if (!ok){
		int ilength; char stringBuffer[STRING_BUFFER];
		glGetShaderInfoLog(shaderID, STRING_BUFFER, &ilength, stringBuffer);
		std::cerr<<"  --> [ERROR] Shader Error : (" << msg << ") " << stringBuffer<< std::endl; 
		return (_shaderStatus=false);    
	}
  
	return true;
}

