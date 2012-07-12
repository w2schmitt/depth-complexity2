#version 400
 
uniform mat4 projectionMat;
uniform mat4 modelViewMat;

in vec3 vertexPos;
 
void main(void) {
    gl_Position = projectionMat * modelViewMat * vec4(vertexPos.xyz, 1.0);
}

