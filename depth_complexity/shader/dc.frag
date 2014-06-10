#version 420

//-- enable opengl extentions
#extension GL_NV_gpu_shader5 : enable
#extension GL_EXT_shader_image_load_store : enable

//-- Whole number pixel offsets
layout(pixel_center_integer) in vec4 gl_FragCoord;

//-- fragment counter texture
coherent uniform layout(size1x32) uimage2D counterBuff;
coherent uniform layout(size2x32) image2D thicknessBuff;
coherent uniform layout(binding = 2, r32f) image2DArray arrayBuff;

//-- screen resolution
uniform ivec2 resolution;
uniform int maxDepth;

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
     // record thickness
     float zval_linear = zLinearization(zval);
     for (int i=0; i<64; i++) {
       vec4 stored_frag = imageLoad(arrayBuff, ivec3(coords, i));
       if (zval_linear < stored_frag.x) {
         // move all
         for (int j=(64-2); j>=i; j--) {
         	vec4 toMoveFrag =  imageLoad(arrayBuff, ivec3(coords, j));
         	if (toMoveFrag.x >= 998) continue;
        	imageStore(arrayBuff, ivec3(coords, j+1), vec4(toMoveFrag.x, 0.0, 0.0, 0.0));
         }	
         imageStore(arrayBuff, ivec3(coords, i), vec4(zval_linear, 0.0, 0.0, 0.0));
       }
   	 }

   	 // record depth complexity
     int abidx = (int)imageAtomicIncWrap( counterBuff, coords, 65535 );
     

     // old thickness
     //float zval_linear = zLinearization(zval);
     //vec4 stored_frag = imageLoad(thicknessBuff, coords);
     //if (zval < stored_frag.x){
     //   imageStore(thicknessBuff, coords, vec4(zval, stored_frag.y, 0.0,0.0));
     //} else if (zval < stored_frag.y){
     //   imageStore(thicknessBuff, coords, vec4(stored_frag.x, zval, 0.0,0.0));
     //}
   }
 
   //-- Discard fragment so nothing is writen to the framebuffer
   discard;
}


