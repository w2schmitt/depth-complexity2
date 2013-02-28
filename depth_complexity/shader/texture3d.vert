#version 130
 
uniform mat4 projectionMat;
uniform mat4 modelViewMat;

in vec3 vertexPos;
in vec3 texCoord;
out vec3 texCoordOut;
 
void main(void) {
    gl_Position = projectionMat * modelViewMat * vec4(vertexPos.xyz, 1.0);
    texCoordOut = texCoord;
}

