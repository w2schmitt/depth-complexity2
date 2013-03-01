#version 330

//-- enable opengl extentions
#extension GL_EXT_gpu_shader4 : enable    //Include support for this extension, which defines usampler2D
//#extension GL_NV_gpu_shader5 : enable
//#extension GL_EXT_shader_image_load_store : enable

//-- fragment counter texture
//coherent uniform layout(size1x32) uimage3D tex3d;
uniform sampler3D tex3d;

in Data {
    vec3 normal;
    vec4 eye;
    smooth vec3 texCoord;
} DataIn;
 
//in vec3 texCoordOut;
out vec4 colorOut;
in vec3 normal;

void main(void){
    //uint texel = texture(tex3d, ivec3(texCoordOut)).x; //imageLoad(tex3d, ivec3(texCoordOut)).x;

    vec4 texColor = texture(tex3d, DataIn.texCoord);
    colorOut = texColor; //vec4(DataIn.texCoord, 1.0); //

}
