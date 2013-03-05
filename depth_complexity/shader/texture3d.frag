#version 330

//-- enable opengl extentions
#extension GL_EXT_gpu_shader4 : enable    //Include support for this extension, which defines usampler2D
//#extension GL_NV_gpu_shader5 : enable
//#extension GL_EXT_shader_image_load_store : enable

//-- fragment counter texture
//coherent uniform layout(size1x32) uimage3D tex3d;
uniform sampler3D tex3d;

//in Data {
//    vec3 normal;
//    vec4 eye;
//    smooth vec3 texCoord;
//} DataIn;
 
//in vec3 texCoordOut;

out vec4 colorOut;

in vec3 texCoord;
in vec3 normal;
//in vec4 color;

void main(void){
    vec4 color;

    // lighting
    vec3 lightDir = vec3(0,0,1);
    vec4 lightDiffuse = vec4(1,1,1,1);
    float NdotL;
    
    // material
    vec4 triDiffuseColor = vec4(0.45f, 0.45f, 0.45f, 0.85f);
    vec4 triAmbientColor = vec4(0.1f, 0.1f, 0.1f, 1.0f);
    vec4 diffuse;

    //uint texel = texture(tex3d, ivec3(texCoordOut)).x; //imageLoad(tex3d, ivec3(texCoordOut)).x;
    vec4 texColor = texture(tex3d, texCoord.yxz);

    
    lightDir = normalize(lightDir);
    NdotL = max(dot(normal, lightDir), 0.0);

    diffuse = texColor * lightDiffuse;
    color = max(NdotL * diffuse, triAmbientColor);
    color.a = 0.75f;

    colorOut = color;

    //colorOut = color; //vec4(DataIn.texCoord, 1.0); //

}
