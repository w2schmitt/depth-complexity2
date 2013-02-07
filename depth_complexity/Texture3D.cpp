/* 
 * File:   Texture3D.cpp
 * Author: wagner
 * 
 * Created on 14 de Setembro de 2012, 10:53
 */

#include "Texture3D.h"

#define WIDTH  512
#define HEIGHT 512
#define DEPTH  512


void Texture3D::CreateTexture3D(unsigned int _w, unsigned int _h, unsigned int _d, unsigned int _ch, float _init) {
    _tex3D = CImg<unsigned int>(WIDTH, HEIGHT,DEPTH,_ch,0);
}

Texture3D::Texture3D(const Texture3D& orig) {
    
}

Texture3D::~Texture3D() {
    
}


const unsigned int colorTableSize=6;

float CTable[colorTableSize][4] = {
    {0.0, 0.0, 1.0, 0.85}, // Azul
    {1.0, 0.0, 1.0, 1.0}, // Roxo/Rosa
    {1.0, 0.5, 0.0, 1.0}, // Laranja
    {1.0, 1.0, 0.0, 1.0}, // Amarelo
    {0.0, 1.0, 0.0, 1.0}, // Verde
    {0.0, 1.0, 0.0, 1.0} // Verde
};

/* Linear Interpolation between color1 and color2*/
void lerpColor(float* newColor, const float *color1, const float *color2, float value){
    
    if (value<=0){
        //memcpy(newColor, color1, sizeof(float*)*3);
        newColor[0] = color1[0];
        newColor[1] = color1[1];
        newColor[2] = color1[2];
        newColor[3] = color1[3];
    }
    else if (value>=1){
        newColor[0] = color2[0];
        newColor[1] = color2[1];
        newColor[2] = color2[2];
        newColor[3] = color2[3];
        //std::cout << value << " - [" << color2[0] << " " << color2[1] << " " << color2[2] << "]\n";
        //memcpy(newColor, color2, sizeof(float*)*3);
    }
    else {
        newColor[0] = (1.0-value)*color1[0] + (value)*color2[0];
        newColor[1] = (1.0-value)*color1[1] + (value)*color2[1];
        newColor[2] = (1.0-value)*color1[2] + (value)*color2[2];
        newColor[3] = (1.0-value)*color1[3] + (value)*color2[3];
    }
}


void findColor(float *pxColor, float normalizedDC){
    
    float intensity;
    float tableIndex = normalizedDC*(float)(colorTableSize-1);
    
    int firstColor = (int)tableIndex;
    
    if (firstColor==colorTableSize-1){
        intensity=1.0;
        firstColor--;
    }
    else{
        intensity = tableIndex - (float)firstColor;
    }
    
    int secondColor = ((int)firstColor)+1;
    lerpColor(pxColor, CTable[firstColor], CTable[secondColor], intensity);
}

GLuint Texture3D::createOpenglTexture3D(int width, int height, int depth, const float* texels){
    //allocate texture name
    GLuint id;
    glGenTextures(1, &id);
    glBindTexture(GL_TEXTURE_3D, id);
    
    glTexEnvf(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_MODULATE);
    glTexParameterf(GL_TEXTURE_3D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
    glTexParameterf(GL_TEXTURE_3D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
    
    glTexParameterf(GL_TEXTURE_3D, GL_TEXTURE_WRAP_S, GL_REPEAT);
    glTexParameterf(GL_TEXTURE_3D, GL_TEXTURE_WRAP_T, GL_REPEAT);
    glTexParameterf(GL_TEXTURE_3D, GL_TEXTURE_WRAP_R, GL_REPEAT);
    
    glTexImage3D(GL_TEXTURE_3D, 0, GL_RGBA, width, height, depth, 0, GL_RGBA, GL_FLOAT, texels);
    
    glBindTexture(GL_TEXTURE_3D, 0);
    
    return id;
    
}




void Texture3D::drawPixel(vec3d pos, unsigned int color){
    
    unsigned int w = _tex3D.width();
    unsigned int h = _tex3D.height();
    unsigned int d = _tex3D.depth();
    
    if ( (pos.x > 0 && pos.x < w) &&
         (pos.y > 0 && pos.y < h) &&
         (pos.z > 0 && pos.z < d)){

        if (color > _tex3D(pos.x, pos.y, pos.z, 0)){
            _tex3D(pos.x, pos.y, pos.z, 0) = color;
            //_tex3D(pos.x, pos.y, pos.z, 1) += 1; // ray counter
        }

    }

}


#define LAYER(r, c)	(WIDTH * HEIGHT * r * c)
// 2->1 dimension mapping function
#define TEXEL2(s, t, c)	(c * (s * WIDTH + t))
// 3->1 dimension mapping function
#define TEXEL3(s, t, r, c)	(TEXEL2(s, t, c) + LAYER(r, c))


void Texture3D::cimg2Tex(unsigned int maxDC){
    // each layer
    unsigned int w = _tex3D.width();
    unsigned int h = _tex3D.height();
    unsigned int d = _tex3D.depth();
    
    //std::cout << w << " " << h << " " << d << std::endl;
    
    unsigned int texsize = w*h*d;    
    float *interlaced_data = new float[4*texsize];
    
    for (unsigned int i=0; i<texsize*4; i++){
        interlaced_data[i] = 0.0f;
    }
    
    
    for (unsigned int r=0; r<d; r++){
        for (unsigned int s=0; s<w; s++){
            for (unsigned int t=0; t<h; t++){
                float pxColor[4] = {0.0, 0.0, 0.0, 0.0};
                float Nray = _tex3D(s,t,r,1);
                Nray = (Nray==0)? 1 : Nray;
                findColor(pxColor, (_tex3D(s,t,r,0)/*/nray*/)/(float)maxDC);
                
                interlaced_data[TEXEL3(t,s,r,4)+0] = pxColor[0];
                interlaced_data[TEXEL3(t,s,r,4)+1] = pxColor[1];
                interlaced_data[TEXEL3(t,s,r,4)+2] = pxColor[2];
                interlaced_data[TEXEL3(t,s,r,4)+3] = pxColor[3];
            }
        }
    }
    
    
    std::cout << "Creating 3D texture..." << std::endl;
    _texID = createOpenglTexture3D(_tex3D.width(), _tex3D.height(), _tex3D.depth(), interlaced_data);
    if (_texID!=0) { std::cout << "[OK] 3D texture succefful" << std::endl;}
    delete[] interlaced_data;

}



void Texture3D::updateTexture3D(Segment line, unsigned int dc){
     
    unsigned int w = _tex3D.width();
    unsigned int h = _tex3D.height();
    unsigned int d = _tex3D.depth();
    
    vec3d t = _aabb.min;
    vec3d size = _aabb.extents();     
    
    vec3d start = line.a;
    start -= t;
    start.x = floor( (start.x/size.x)*(w-1) + 0.5);
    start.y = floor( (start.y/size.y)*(h-1) + 0.5);
    start.z = floor( (start.z/size.z)*(d-1) + 0.5);
    
    vec3d end = line.b;
    end -= t;
    end.x = floor( (end.x/size.x)*(w-1) + 0.5);
    end.y = floor( (end.y/size.y)*(h-1) + 0.5);
    end.z = floor( (end.z/size.z)*(d-1) + 0.5);
    
    // rasterize ray into a 3D texture, using dc as color
    bresenham_line_3D(start, end, dc);  
}





/* find maximum of a and b */
#define MAX(a,b) (((a)>(b))?(a):(b))
/* absolute value of a */
#define ABS(a) (((a)<0) ? -(a) : (a))
/* take sign of a, either -1, 0, or 1 */
#define ZSGN(a) (((a)<0) ? -1 : (a)>0 ? 1 : 0)

/*
 * line3d was dervied from DigitalLine.c published as "Digital Line Drawing"
 * by Paul Heckbert from "Graphics Gems", Academic Press, 1990
 * 
 * 3D modifications by Bob Pendleton. The original source code was in the public
 * domain, the author of the 3D version places his modifications in the
 * public domain as well.
 * 
 * line3d uses Bresenham's algorithm to generate the 3 dimensional points on a
 * line from (x1, y1, z1) to (x2, y2, z2)
 * 
 */

void Texture3D::bresenham_line_3D(vec3d &v1, vec3d &v2, unsigned int color){
    int x1, y1, x2, y2, z1, z2;

    x1 = v1.x; y1 = v1.y; z1 = v1.z;
    x2 = v2.x; y2 = v2.y; z2 = v2.z;
    
    int xd, yd, zd;
    int x, y, z;
    int ax, ay, az;
    int sx, sy, sz;
    int dx, dy, dz;

    dx = x2 - x1;
    dy = y2 - y1;
    dz = z2 - z1;

    ax = ABS(dx) << 1;
    ay = ABS(dy) << 1;
    az = ABS(dz) << 1;

    sx = ZSGN(dx);
    sy = ZSGN(dy);
    sz = ZSGN(dz);

    x = x1;
    y = y1;
    z = z1;

    if (ax >= MAX(ay, az))            /* x dominant */
    {
        yd = ay - (ax >> 1);
        zd = az - (ax >> 1);
        for (;;)
        {
            drawPixel(vec3d(x, y, z), color);
            if (x == x2)
            {
                return;
            }

            if (yd >= 0)
            {
                y += sy;
                yd -= ax;
            }

            if (zd >= 0)
            {
                z += sz;
                zd -= ax;
            }

            x += sx;
            yd += ay;
            zd += az;
        }
    }
    else if (ay >= MAX(ax, az))            /* y dominant */
    {
        xd = ax - (ay >> 1);
        zd = az - (ay >> 1);
        for (;;)
        {
            drawPixel(vec3d(x, y, z),color);
            if (y == y2)
            {
                return;
            }

            if (xd >= 0)
            {
                x += sx;
                xd -= ay;
            }

            if (zd >= 0)
            {
                z += sz;
                zd -= ay;
            }

            y += sy;
            xd += ax;
            zd += az;
        }
    }
    else if (az >= MAX(ax, ay))            /* z dominant */
    {
        xd = ax - (az >> 1);
        yd = ay - (az >> 1);
        for (;;)
        {
            drawPixel(vec3d(x, y, z),color);
            if (z == z2)
            {
                return;
            }

            if (xd >= 0)
            {
                x += sx;
                xd -= az;
            }

            if (yd >= 0)
            {
                y += sy;
                yd -= az;
            }

            z += sz;
            xd += ax;
            yd += ay;
        }
    }
}

