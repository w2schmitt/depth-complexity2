//FRAGMENT SHADER
//#version 130

//Include support for this extension, which defines usampler2D
//#extension GL_EXT_gpu_shader4 : enable

void main(){
	
	//vec4(0.15f, 0.0f, 0.0f, 0.0f);
	gl_FragData[0] = gl_Color;
}
