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

using namespace std;

// STRUCTURES ------------------------------>
struct Hist {
    int index;
    std::string name;
    std::vector<int> dc;
    std::vector<unsigned long long int> rays;
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
std::vector<Comp> compList;                         // bhattacharya coef.
std::vector<float> distanceMatrix;


// FUNCTIONS DECLARATION ------------------->
std::vector<std::string>        searchFiles(std::string folder);
void                            openFiles(std::string filepath);
bool                            compFunction (Hist h);
bool                            sortFunction  (Comp c1, Comp c2);
void                            matchShape(std::string name);
double                          computeBC(Hist h1, Hist h2);
void                            createMatrixDistance();
bool                            compIndexFunction (Hist h1, Hist h2);
Hist                            findHist(int index);

// DEFINITIONS ----------------------------->
int main(int argc, char** argv) {


    openFiles ("Histograms//");    
    std::sort (histlist.begin(), histlist.end(), compIndexFunction);  
    
    createMatrixDistance();
    
    FILE * pFile;
    pFile = fopen ("distanceMatrix.bin", "wb");
    
    fwrite (distanceMatrix.data(), sizeof(float), distanceMatrix.size(), pFile);
    fclose (pFile);
    
    return 0;
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
    distanceMatrix.resize(histlist.size()*histlist.size(), 0);
    
    for (unsigned int i=0; i < histlist.size(); i++){
        Hist h1 = findHist(i);
        for (unsigned int j=0; j < histlist.size(); j++){
            //std::cout << histlist.size()  << " - " << j << std::endl;
            Hist h2 = findHist(j);
            //std::cout << h2.index << std::endl;
            float bc = computeBC(h1,h2);
            float similarity = -log(bc);
            similarity = (similarity<1.E-7)? 0.00f : similarity;
            
            if (i==1119 && j==1134){
                std::clog << similarity << std::endl;
            }
            distanceMatrix[i*histlist.size() + j] = similarity;
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
    unsigned long long int sumCoef1 = 0, sumCoef2 = 0;
    unsigned int maxSize;
    
    
    maxSize = (h1.dc.size() > h2.dc.size())? h1.dc.size() : h2.dc.size();
    
    
    std::vector<double> coef1;//(h1.rays);
    std::vector<double> coef2;//(h2.rays);
    
    // normalize coefs
    sumCoef1 = std::accumulate(h1.rays.begin(), h1.rays.end(), 0);
    sumCoef2 = std::accumulate(h2.rays.begin(), h2.rays.end(), 0);
    
    // coef must have the same domain...
    coef1.resize(maxSize, 0);
    coef2.resize(maxSize, 0);
    
    // and normalized
    for (unsigned int i=0; i<h1.rays.size(); i++){
        coef1[i] = h1.rays[i]/(double)sumCoef1;
    }
    for (unsigned int i=0; i<h2.rays.size(); i++){
        coef2[i] = h2.rays[i]/(double)sumCoef2;
    }
    
    // find bhat... coeficient
    for (unsigned int i=0; i<maxSize; i++){
        bc += sqrt(coef1[i] * coef2[i]);        
    }
    
    return bc;
}




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



