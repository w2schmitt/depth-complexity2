
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <dirent.h>
#include <string.h>
#include <vector>
#include <algorithm>
#include <math.h>
#include <numeric>

// STRUCTURES ------------------------------>
struct Hist {
    int index;
    std::string name;
    std::vector<int> dc;
    std::vector<unsigned long long int> rays;
    std::vector<double> normalized_rays;
};

// GLOBALS --------------------------------->
std::vector<std::string> filelist;
std::vector<Hist> histlist;
unsigned int i_max=0, max=0;

// DECLARATIONS ---------------------------->
std::vector<std::string>        searchFiles(std::string folder);
void                            openFiles(std::string filepath);
bool                            compIndexFunction (Hist h1, Hist h2);
void                            getMaximumSize(std::vector<Hist> h);
void                            printdata();
void                            normalize();

// DEFINITIONS ----------------------------->
int main(int argc, char** argv) {
    openFiles ("Histograms//");    
    std::sort (histlist.begin(), histlist.end(), compIndexFunction); 
    
    
    getMaximumSize(histlist);
    normalize();
    printdata();
    
    return 0;
}

void getMaximumSize(std::vector<Hist> h){
    unsigned int tempmax=0;
    
    //std::cout << "hsize: " << h.size() << std::endl;
    for (unsigned int i=0; i<h.size(); i++){
        tempmax = h[i].dc.size();
        if (tempmax > max){
            max = tempmax;
            i_max = i;
        }
    } 
}

void printdata(){
    //std::string str="";
    //std::cout << "max: " << max << std::endl;
    for (unsigned int i=0; i<max; i++){
        for (unsigned int j=0; j<=histlist.size(); j++){
            if (i==0){
                if (j==0){
                    printf("Intersections");
                } else 
                    printf(",m%d", j-1);
            } else {
                if (j==0)
                    printf("%d", i-1);
                else                   
                    if (i <= histlist[j-1].normalized_rays.size())
                        printf(",%f", histlist[j-1].normalized_rays[i-1]);
                    else
                        printf(",0");
            }
        }
        printf("\n");
    
        
    }
    
}

void normalize(){
    unsigned long long int sumCoef1 = 0;
    
    for (unsigned int i=0; i<histlist.size(); i++){
        sumCoef1 = 0;
        sumCoef1 = std::accumulate(histlist[i].rays.begin(), histlist[i].rays.end(), 0);
        
        histlist[i].normalized_rays.resize(max,0);
        
        for (unsigned int j=0; j<histlist[i].rays.size(); j++){
            if (sumCoef1!=0)
                histlist[i].normalized_rays[j] = histlist[i].rays[j]/(double)sumCoef1;
        }
    }
}


bool compIndexFunction (Hist h1, Hist h2){
    return h1.index < h2.index;
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

