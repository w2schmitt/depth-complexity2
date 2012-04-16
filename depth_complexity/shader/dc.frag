#version 400

//-- enable opengl extentions
#extension GL_NV_gpu_shader5 : enable
#extension GL_EXT_shader_image_load_store : enable
#extension GL_EXT_bindable_uniform : enable
#extension GL_NV_shader_buffer_load : enable
#extension GL_NV_shader_buffer_store : enable

//-- Whole number pixel offsets
layout(pixel_center_integer) in vec4 gl_FragCoord;

//-- fragment counter texture
coherent uniform layout(size1x32) uimage2D counterBuff;

//-- Output fragment color
out vec4 outFragColor;

void main(void){
   //-- get coords
   ivec2 coords = ivec2(gl_FragCoord.xy);
   //memoryBarrier();
   //imageStore(counterBuff, coords, ivec4(20000));
   //if (imageLoad(counterBuff,coords) == uvec4(32550).r){
   //if(coords.x>=0 && coords.y>=0 && coords.x<512 && coords.y<512 ){
     int abidx = (int)imageAtomicIncWrap( counterBuff, coords, 25000 );
     //outFragColor = vec4(0.5f); //vec4(imageLoad(counterBuff, coords));
   //}
   //else{
   //  outFragColor = vec4(0.0f);
   //}
   discard;
   //}
   //else
   //   discard;
   //outFragColor = vec4(150);
   //-- Check if we are inside the framebuffer
   
     
    // imageStore(counterBuff, coords, uvec4(1,1,1,1));
		//-- increment counter
		//int abidx = (int)imageAtomicIncWrap( counterBuff, coords, 2000 );
    //imageStore(counterBuff, coords, uvec4(0));
    //int abidx = (int)imageAtomicAdd( counterBuff, coords, 1 );
  //}
   
   //-- Discard fragment so nothing is writen to the framebuffer
   //discard;
}
