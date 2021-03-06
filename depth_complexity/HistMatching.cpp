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
    std::string name;
    std::vector<int> dc;
    std::vector<unsigned long long int> rays;
};

struct Comp {
    std::string firstName;
    std::string secondName;
    double bc;
    float similarity;
};


// GLOBALS --------------------------------->
std::vector<std::string> filelist;
std::vector<Hist> histlist;
std::string fileInput;
std::vector<Comp> compList;                         // bhattacharya coef.


// FUNCTIONS DECLARATION ------------------->
std::vector<std::string>        searchFiles(std::string folder);
void                            openFiles(std::string filepath);
bool                            compFunction (Hist h);
bool                            sortFunction  (Comp c1, Comp c2);
void                            matchShape(std::string name);
double                          computeBC(Hist h1, Hist h2);

// DEFINITIONS ----------------------------->
int main(int argc, char** argv) {

    // input model
    if (argc==2){
        fileInput = argv[1];         
    } else {
        std::cerr << "--->      ./hist_matching <inputFile>" << std::endl;
        return -1;
    }
    std::cout << "FILE INPUT: " << fileInput << std::endl << std::endl;  
    // read histogram files
    openFiles("Histograms//");    
    matchShape(fileInput);
    
    std::sort (compList.begin(), compList.end(), sortFunction);
    
    // print similarity coefficient
    for (unsigned int i=0; i<compList.size(); i++){
        std::cout << compList[i].secondName << " : " << compList[i].similarity << std::endl;
    }
    
    return 0;
}

bool compFunction (Hist h){
    if (h.name.compare(fileInput) == 0){
        return true;
    }
    return false;
}

bool sortFunction (Comp c1, Comp c2){
    return (c1.similarity < c2.similarity);
}

void matchShape(std::string name){
    // check if name is registred    
    std::vector<Hist>::iterator it = std::find_if(histlist.begin(), histlist.end(), compFunction);
    //std::clog << histlist[10].name << std::endl;
    
    if (it != histlist.end()){
        for (unsigned int i=0; i<histlist.size(); i++){
            Comp c;
            c.firstName = it->name;
            c.secondName = histlist[i].name;
            c.bc = computeBC(*it, histlist[i] );
            c.similarity = -log(c.bc);
            c.similarity = (c.similarity<1E-7)? 0 : c.similarity;
            
            compList.push_back(c);
        }
    } else {
        std::cerr << "Hist not found." << std::endl;
    }
   
}

double computeBC(Hist h1, Hist h2){
    double bc = 0.0;
    unsigned long long int sumCoef1 = 0, sumCoef2 = 0;
    unsigned int maxSize;
    
    maxSize = (h1.dc.size() > h2.dc.size())? h1.dc.size() : h2.dc.size();
    
    std::vector<double> coef1;//(h1.rays);
    std::vector<double> coef2;//(h2.rays);
    
    // normalize coefs
    sumCoef1 = std::accumulate(h1.rays.begin(), h1.rays.end(), 0);
    sumCoef2 = std::accumulate(h2.rays.begin(), h2.rays.end(), 0);
    
    //std::cout << sumCoef1 << std::endl;
    
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
    
    //std::cout << bc << std::endl;
    
    
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
            h.name = filelist[i];
           
            //std::cout << "Reading: " << h.name << std::endl;
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



