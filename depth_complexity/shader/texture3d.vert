#version 330

// opengl Matrices 
uniform mat4  gl_ModelViewProjectionMatrix;
uniform mat3  gl_NormalMatrix;


// Inputs
attribute vec3 gl_MultiTexCoord0;
attribute vec4 gl_Vertex;      // local space
attribute vec3 gl_Normal;     // local space
attribute vec4 gl_Color;


// Outputs
out vec3 texCoord;
out vec3 normal;
//out vec4 color;


//in vec3 texCoord;

void main(void) {
    gl_Position = gl_ModelViewProjectionMatrix * vec4(gl_Vertex.xyz, 1.0);

    normal = normalize(gl_NormalMatrix * gl_Normal);
    texCoord = gl_MultiTexCoord0;

    

    

    //color = gl_FrontColor;
}

