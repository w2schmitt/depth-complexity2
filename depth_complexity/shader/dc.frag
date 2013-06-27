#version 420

//-- enable opengl extentions
#extension GL_NV_gpu_shader5 : enable
#extension GL_EXT_shader_image_load_store : enable

//-- Whole number pixel offsets
layout(pixel_center_integer) in vec4 gl_FragCoord;

//-- fragment counter texture
coherent uniform layout(size1x32) uimage2D counterBuff;
coherent uniform layout(size2x32) image2D thicknessBuff;

//-- screen resolution
uniform ivec2 resolution;

//-- znear and far values
uniform vec2 zvalues;

float zLinearization(float zvalue){
    return (2.0 * zvalues.x) / (zvalues.y + zvalues.x - zvalue * (zvalues.y-zvalues.x));
}

//void ComputeThickness(z){

//}

void main(void){
   //-- get coords
   ivec2 coords = ivec2(gl_FragCoord.xy);
   float zval = gl_FragCoord.z;

   if(coords.x>=0 && coords.y>=0 && coords.x<resolution.x && coords.y<resolution.y ){
     int abidx = (int)imageAtomicIncWrap( counterBuff, coords, 65535 );
     float zval_linear = zLinearization(zval);
     imageStore(thicknessBuff, coords, vec4(zval_linear, zval_linear, 0.1,0.1));
   }
  
   //-- Discard fragment so nothing is writen to the framebuffer
   discard;
}


