#version 400

//-- enable opengl extentions
#extension GL_NV_gpu_shader5 : enable
#extension GL_EXT_shader_image_load_store : enable

//-- Whole number pixel offsets
layout(pixel_center_integer) in vec4 gl_FragCoord;

//-- fragment counter texture
coherent uniform layout(size1x32) uimage2D counterBuff;
uniform ivec2 resolution;

void main(void){
   //-- get coords
   ivec2 coords = ivec2(gl_FragCoord.xy);

   if(coords.x>=0 && coords.y>=0 && coords.x<resolution.x && coords.y<resolution.y ){
     //int abidx = (int)imageAtomicIncWrap( counterBuff, coords, 65535 );
     imageStore(counterBuff, coords, uvec4(1));
     //imageAtomicAdd(counterBuff, coords, 1);
   }
  
   //-- Discard fragment so nothing is writen to the framebuffer
   discard;
}
