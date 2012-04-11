//VERTEX SHADER
//#version 130

void main(){
	gl_FrontColor = gl_Color;
	gl_Position = ftransform();
    //TexCoord0 = InTexCoord0;
}
