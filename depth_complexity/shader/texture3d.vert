#version 330
 
uniform mat4 projectionMat;
uniform mat4 modelViewMat;

//in vec3 vertexPos;
//in vec3 texCoord;
//out vec3 texCoordOut;

out Data {
    vec3 normal;
    vec4 eye;
    smooth vec3 texCoord;
} DataOut;

attribute vec3 gl_MultiTexCoord0;

attribute vec4 position;   // local space
attribute vec3 gl_Normal;     // local space
//in vec3 texCoord;

void main(void) {
    gl_Position = projectionMat * modelViewMat * vec4(position.xyz, 1.0);
    DataOut.texCoord = gl_MultiTexCoord0;
}

