/* 
 * File:   HistMatching.cpp
 * Author: wagner
 *
 * Created on 19 de Maio de 2013, 14:42
 */

#include <cstdlib>
#include <iostream>
#include <fstream>
#include <dirent.h>
#include <string.h>
#include <vector>
#include <algorithm>
#include <math.h>
#include <numeric>
#include <float.h>

using namespace std;

// STRUCTURES ------------------------------>
struct Hist {
    int index;
    std::string name;
    std::vector<int> dc;
    std::vector<float> rays;
};

struct Comp {
    std::string firstName;
    std::string secondName;
    double bc;
};

// GLOBALS --------------------------------->
std::vector<std::string> filelist;
std::vector<Hist> histlist;
std::string fileInput;
std::string output_dir;
std::vector<Comp> compList;                         // bhattacharya coef.
std::vector< std::vector<float> > distanceMatrix;
bool useWeights = false;


// FUNCTIONS DECLARATION ------------------->
std::vector<std::string>        searchFiles(std::string folder);
void                            openFiles(std::string filepath);
bool                            compFunction (Hist h);
bool                            sortFunction  (Comp c1, Comp c2);
void                            matchShape(std::string name);
double                          computeBC(Hist h1, Hist h2);
void                            createMatrixDistance();
void                            printdataInJs();
bool                            compIndexFunction (Hist h1, Hist h2);
bool                            IsNumber(double x);
bool                            IsFiniteNumber(double x);
Hist                            findHist(int index);

// DEFINITIONS ----------------------------->
int main(int argc, char** argv) {
    //"thickness//"

    openFiles (argv[1]);    
    std::sort (histlist.begin(), histlist.end(), compIndexFunction);  
    
    output_dir = argv[2];
    
    //allocate matrix size
    unsigned int NoModels = histlist.size();
    distanceMatrix.resize(NoModels);
    for (unsigned int i=0; i<NoModels; i++){
        
        distanceMatrix[i].resize(NoModels,0);
    }
    
    std::cout << "creating matrix\n" << std::endl;
    createMatrixDistance();
    
    std::cout << "printing data\n" << std::endl;
    printdataInJs();
    //FILE * pFile;
    //pFile = fopen ("distanceMatrix.bin", "wb");
    
    //fwrite (distanceMatrix.data(), sizeof(float), distanceMatrix.size(), pFile);
    //fclose (pFile);
    
    return 0;
}

bool IsNumber(double x) 
{
    // This looks like it should always be true, 
    // but it's false if x is a NaN.
    return (x == x); 
}
bool IsFiniteNumber(double x) 
{
    return (x <= DBL_MAX && x >= -DBL_MAX); 
} 

void printdataInJs(){
    for (unsigned int i=0; i<distanceMatrix.size(); i++){
        char filename[100];
        sprintf(filename,"%s//m%d.js", output_dir.c_str(), i);
        
        //create new js file
        FILE * pFile;        
        pFile = fopen (filename,"w");
        
        char comma=' ';
        fprintf(pFile, "m%d = [", i);
        for (unsigned int j=0; j<distanceMatrix[i].size(); j++){
                fprintf(pFile, "%c%f", comma, distanceMatrix[i][j]);
                comma=',';
        }
        fprintf(pFile, "];");
        fclose (pFile);
    }
}

bool compIndexFunction (Hist h1, Hist h2){
    return h1.index < h2.index;
}

Hist findHist(int index){
    for (unsigned int i=0; i < histlist.size(); i++){
        if (histlist[i].index == index)
            return histlist[i];
    }
    
    return histlist[0];
}

bool compFunction (Hist h){
    if (h.name.compare(fileInput) == 0){
        return true;
    }
    return false;
}

bool sortFunction (Comp c1, Comp c2){
    return (c1.bc > c2.bc);
}

void createMatrixDistance(){
    //distanceMatrix.resize(histlist.size()*histlist.size(), 0);
    
    for (unsigned int i=0; i < histlist.size(); i++){
        std::cout << "computing... " << i << std::endl; 
        Hist h1 = findHist(i);
        for (unsigned int j=0; j < histlist.size(); j++){
            //std::cout << histlist.size()  << " - " << j << std::endl;
            Hist h2 = findHist(j);
            //std::cout << h2.index << std::endl;
            float bc = computeBC(h1,h2);
            //std::cout << bc << std::endl;
            float similarity = -log(bc);
       
            if (IsNumber(similarity) && IsFiniteNumber(similarity)){
                similarity = (similarity<1.E-7)? 0.00f : similarity;
            } else {
                //std::clog << "Warning: histogram " << i << " and " << j << " --> " << similarity << std::endl;
                similarity = 999.0;

            }
                 //std::cout << bc << " - " << similarity << std::endl;
            //std::cin.get();
            
            //if (i==1119 && j==1134){
            //    std::clog << similarity << std::endl;
            //}
            distanceMatrix[i][j] = similarity;
        }        
    }
}

void matchShape(std::string name){
    
    // check if name is registred
    std::vector<Hist>::iterator it = std::find_if(histlist.begin(), histlist.end(), compFunction);
    
    for (unsigned int i=0; i<histlist.size(); i++){
        Comp c;
        c.firstName = it->name;
        c.secondName = histlist[i].name;
        c.bc = computeBC(*it, histlist[i] );
        compList.push_back(c);
    }
}

double computeBC(Hist h1, Hist h2){
    //std::cout << "pica2\n";
    double bc = 0.0;
    unsigned int minSize;
    //unsigned long long int sumCoef1 = 0, sumCoef2 = 0;
    //unsigned int maxSize;
    
    
    //maxSize = (h1.dc.size() > h2.dc.size())? h1.dc.size() : h2.dc.size();
    minSize = (h1.rays.size() > h2.rays.size())? h2.rays.size() : h1.rays.size();
    //std::cout << minSize << std::endl;
    //for (unsigned int i=0; i<h1.rays.size(); i++){
    //    std::cout << h1.rays[i] << " ";
    //}
    //std::cout << std::endl;

    //std::cout <<  << std::endl;
    
    //std::vector<double> coef1;//(h1.rays);
    //std::vector<double> coef2;//(h2.rays);
    
    // normalize coefs
    //sumCoef1 = std::accumulate(h1.rays.begin(), h1.rays.end(), 0);
    //sumCoef2 = std::accumulate(h2.rays.begin(), h2.rays.end(), 0);
    
    // coef must have the same domain...
    //coef1.resize(maxSize, 0);
    //coef2.resize(maxSize, 0);
    
    // and normalized
    //for (unsigned int i=0; i<h1.rays.size(); i++){
    //    coef1[i] = h1.rays[i]/(double)sumCoef1;
   // }
    //for (unsigned int i=0; i<h2.rays.size(); i++){
    //    coef2[i] = h2.rays[i]/(double)sumCoef2;
    //}
    
    // find bhat... coeficient
    for (unsigned int i=0; i<minSize; i++){
        bc += (sqrt(h1.rays[i] * h2.rays[i])); //*(useWeights? i*10: 1));        
    }

    //std::cout << bc << std::endl;
    //std::cin.get();
    
    return bc;
}



/*
void openFiles(std::string filepath){
    
    filelist = searchFiles(filepath);
    
    char label1[50], label2[50];    
    unsigned long long int Nrays;
    int dc;
    
    for (unsigned int i=0; i<filelist.size(); i++){
         std::filebuf fb; 
         if (fb.open ((filepath+filelist[i]).c_str() ,std::ios::in)){
            std::istream in(&fb);
            
            // new histogram
            Hist h;
            char temp;
            h.name = filelist[i];
            sscanf(filelist[i].c_str(), "%c%d", &temp, &h.index );           
           
            //std::clog << "Reading: " << filelist[i] << std::endl;
            in >> label1 >> label2;            
            while (in){            
                in >> dc >> Nrays;
                if (!in) break;
                h.dc.push_back(dc);
                h.rays.push_back(Nrays);
                //std::cout << dc << " " << Nrays << std::endl;
            }     
            histlist.push_back(h);
        } else {
             std::cerr << "Cannot open " << filepath+filelist[i] << std::endl;
        }
    }
    
    std::clog << "loaded " << histlist.size() << " histograms" << std::endl; 
}*/

    void openFiles(std::string filepath){
    
    filelist = searchFiles(filepath);
    
    char var_name[50], label2[50];    
    //unsigned long long int Nrays;
    //int dc,thickness;
    //unsigned long int value;
    int dc = 0;
    float frequency = 0.0f;
    int v=0;
    bool success;
    
    for (unsigned int i=0; i<filelist.size(); i++){
        dc = 0;
        FILE* pFile = fopen( (filepath+filelist[i]).c_str(), "r+" );
        if (pFile){
            success = true;     
            // new histogram
            Hist h;
            char temp;
            h.name = filelist[i];
            v = sscanf(filelist[i].c_str(), "%c%d", &temp, &h.index );           

            v = fscanf(pFile, "%s = [ ", var_name);
            char c = ',';
            while ( c==','){            
                v = fscanf(pFile,"%f%c", &frequency, &c);
                h.dc.push_back(dc++);
                h.rays.push_back(frequency);
                //std::cout << h.name << frequency << std::endl;
                //std::cin.get();
                //std::pair<int,int> p = std::make_pair(thickness, dc);
                //h->values.insert( std::pair< std::pair<int,int>, unsigned long> (p,value));
                if (v<2){
                    success=false;
                    break;
                }
            }   

            if (success){  
                //printf("loaded %d hist: %s\n", i,var_name);
                //std::cout << "testando" << std::endl;
                histlist.push_back(h);
                //std::cout << "testando2" << std::endl;
            } else {
                printf("error loading %s\n", var_name);
            }
            fclose(pFile);
        } else {
             std::cerr << "Cannot open " << filepath+filelist[i] << std::endl;
        }
    }
    
    std::clog << "loaded " << histlist.size() << " histograms" << std::endl; 
}

std::vector<std::string> searchFiles(std::string folder){
    DIR *dir;
    struct dirent *ent;
    std::vector<std::string> f;
    
    if ((dir = opendir (folder.c_str())) != NULL) {
        /* print all the files and directories within directory */
        while ((ent = readdir (dir)) != NULL) {
            if ( (strcmp(ent->d_name,".")!=0) && (strcmp(ent->d_name,"..")!=0) ){
                //std::cout << ent->d_name << std::endl;
                f.push_back(std::string(ent->d_name));
            }
    }
    closedir (dir);
    } else {
        /* could not open directory */
        std::cerr << "Cannot open Directory\n";
    }
    
    return f;
}




