#version 130

//-- enable opengl extentions
#extension GL_EXT_gpu_shader4 : enable    //Include support for this extension, which defines usampler2D
//#extension GL_NV_gpu_shader5 : enable
//#extension GL_EXT_shader_image_load_store : enable

//-- fragment counter texture
//coherent uniform layout(size1x32) uimage3D tex3d;
uniform sampler3D tex3d;

in vec3 texCoordOut;
out vec4 colorOut;

void main(void){
    //uint texel = texture(tex3d, ivec3(texCoordOut)).x; //imageLoad(tex3d, ivec3(texCoordOut)).x;

    vec4 texel = texture(tex3d, texCoordOut);
    colorOut = texel;

    //if (texel < 2U)
    //    colorOut = vec4(0.0, 0.0, 0.0, 1.0);
    //else if (texel > 12U)
    //   colorOut = vec4(1.0, 0.0, 0.0, 1.0);
    //else 
    //    colorOut = vec4(0.0, 0.0, 1.0, 1.0);
}
