#version 330

//-- enable opengl extentions
#extension GL_EXT_gpu_shader4 : enable    //Include support for this extension, which defines usampler2D


//-- fragment counter texture
uniform usampler3D tex3d;
uniform uint maxDC;
uniform uint minValue;
uniform uint maxValue;

//-- output
out vec4 colorOut;

//-- inputs
in vec3 texCoord;
in vec3 normal;


//-- contants
vec4 ColorTable[6] = vec4[6]( 
    vec4(0.0, 0.0, 1.0, 0.85),
    vec4(1.0, 0.0, 1.0, 1.0),
    vec4(1.0, 0.5, 0.0, 1.0),
    vec4(1.0, 1.0, 0.0, 1.0),
    vec4(0.0, 1.0, 0.0, 1.0),
    vec4(0.0, 1.0, 0.0, 1.0) 
);

vec4 lerpColor(vec4 color1, vec4 color2, float value){
     if (value<=0.0){
        return vec4(color1);
    }
    else if (value>=1.0){
        return vec4(color2);
    }
    else {
        return vec4((1.0-value)*color1 + (value)*color2);
    }
}

vec4 findColor(float normalizedDC){
    
    float intensity;
    float tableIndex = normalizedDC*float((ColorTable.length()-1));
    
    int firstColor = int(tableIndex);
    
    if (firstColor==ColorTable.length()-1){
        intensity=1.0;
        firstColor--;
    }
    else{
        intensity = tableIndex - float(firstColor);
    }
    
    int secondColor = int(firstColor)+1;
    return lerpColor(ColorTable[firstColor], ColorTable[secondColor], intensity);
    
}

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
    
    //vec4 texColor;
    uint texel = texture(tex3d, texCoord).x; //imageLoad(tex3d, ivec3(texCoordOut)).x;
    //vec4 texColor = texture(tex3d, texCoord).x;

    float normalizedDC;
    if (texel < minValue || texel > maxValue){
        normalizedDC = 0.0;
    } else {    
        normalizedDC = float(texel)/float(maxDC);
    }
     
    vec4 texColor = findColor(normalizedDC);
    
    lightDir = normalize(lightDir);
    NdotL = max(dot(normal, lightDir), 0.0);

    diffuse = texColor * lightDiffuse;
    color = max(NdotL * diffuse, triAmbientColor);
    //color.a = 0.75f;

    colorOut = color;

    //colorOut = color; //vec4(DataIn.texCoord, 1.0); //

}
