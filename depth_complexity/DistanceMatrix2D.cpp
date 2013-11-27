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
#include <utility> 
#include <algorithm>
#include <math.h>
#include <numeric>
#include <float.h>
#include <map>

using namespace std;

enum TYPE {
    EMD,
    BASIC,
    PROPORTION,
    PROPORTION_2
};

TYPE algorithm = PROPORTION_2;


// STRUCTURES ------------------------------>
struct Hist {
    int index;
    std::string name;
    std::vector<int> dc;
    std::vector<unsigned long long int> rays;
};

struct Hist2D {
    int index;
    std::string name;
    std::map< std::pair<int,int>, unsigned long int> values; 
    vector<unsigned long> sumByDc;
    unsigned long totalSum;
};

struct Comp {
    std::string firstName;
    std::string secondName;
    double bc;
};

// GLOBALS --------------------------------->
std::vector<std::string> filelist;
std::vector<Hist> histlist;
std::vector<Hist2D*> histlist2D;
std::string fileInput;
std::string output_dir;
std::vector<Comp> compList;                         // bhattacharya coef.
std::vector< std::vector<float> > distanceMatrix;
int numOfFiles = 0;
bool useWeights = false;



//std::vector<double> similarity;


// FUNCTIONS DECLARATION ------------------->
std::vector<std::string>        searchFiles(std::string folder);
void                            openFiles(std::string filepath);
bool                            compFunction (Hist h);
bool                            sortFunction  (Comp c1, Comp c2);
void                            matchShape(std::string name);
double                          computeBC(Hist h1, Hist h2);
double                          computeSimilarity2D(Hist2D* h1, Hist2D* h2);
double                          computeSimilarity2DByDC(Hist2D* h1,Hist2D* h2);
double                          proportionalSimilarity(Hist2D *h1, Hist2D *h2);
void                            createMatrixDistance();
bool                            compIndexFunction (Hist2D *h1, Hist2D *h2);
bool                            IsNumber(double x);
bool                            IsFiniteNumber(double x);
Hist2D*                         findHist(int index);
void                            saveSimilarityToJs(int file, std::vector<double> similarityList);
double                          proportionalSimilarity2(Hist2D *h1, Hist2D *h2);


// DEFINITIONS ----------------------------->
int main(int argc, char** argv) {
    
    numOfFiles = atoi(argv[1]);
    histlist2D.resize(numOfFiles, NULL); 

    std::cout << "Loading Histograms..." << std::endl;
    openFiles (argv[2]);

    //std::cout << "Sorting Histograms..." << std::endl;
    //std::sort (histlist2D.begin(), histlist2D.end(), compIndexFunction);  
    

    //std::cout << "Calculating distance..." << std::endl;
    output_dir = argv[3];

    createMatrixDistance();

    
    
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
/*
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
*/
bool compIndexFunction (Hist2D* h1, Hist2D* h2){
    return h1->index < h2->index;
}

Hist2D* findHist(int index){
    for (unsigned int i=0; i < histlist2D.size(); i++){
        if (histlist2D[i]->index == index){
            return histlist2D[i];
        }
    }

    return NULL;
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
    std::vector<double> similarityList(numOfFiles);

    for (unsigned int i=0; i < numOfFiles; i++){        
        //similarityList.clear();

        Hist2D *h1 = histlist2D[i];
        if (!h1){
            similarityList.assign(numOfFiles, 9.0);
            saveSimilarityToJs(i, similarityList);
            continue;
        }
        std::clog << "similarity for [m" << i << ".js]" << std::endl; 

        for (unsigned int j=0; j < numOfFiles; j++){
            Hist2D *h2 = histlist2D[j];
            double similarity = 99999.0;
            if (h2){
                //

                if (algorithm == PROPORTION){
                    similarity = proportionalSimilarity(h1,h2);
                } else if (algorithm == BASIC){
                    similarity = computeSimilarity2D(h1,h2);
                    //similarity = computeSimilarity2DByDC(h1,h2);
                } else if (algorithm == PROPORTION_2){
                    similarity = proportionalSimilarity2(h1,h2);
                }
                

                if (IsNumber(similarity) && IsFiniteNumber(similarity)){
                    similarity = (similarity<1.E-7)? 0.00f : similarity;
                } else {
                    similarity = 999.0;
                }
            }
            similarityList[j] = similarity;
        }  

        saveSimilarityToJs(i, similarityList);    
    }
}



void saveSimilarityToJs(int file, std::vector<double> similarityList){
    
    char filename[100];
    sprintf(filename,"%s//m%d.js", output_dir.c_str(), file);

    FILE * pFile;        
    pFile = fopen (filename,"w");

    char comma=' ';
    fprintf(pFile, "m%d = [", file); 
    for (unsigned int j=0; j<similarityList.size(); j++){
        fprintf(pFile, "%c%f", comma, similarityList[j]);
        comma=',';
    }

    fprintf(pFile, "];");
    fclose (pFile);
}

unsigned int selectedDC = 1;


int addMaps(int total, const std::map< std::pair<int,int>, unsigned long int >::value_type &data){
    return total + data.second;
}

int addMapsByDC(int total, const std::map< std::pair<int,int>, unsigned long int >::value_type &data){
    if (data.first.first == selectedDC)
        return total + data.second;
    return total;
}

double proportionalSimilarity2(Hist2D *h1, Hist2D *h2){
    double similarity = 0.0;
    unsigned long long int sumCoef1 = 0, sumCoef2 = 0;

    //sumCoef1 = std::accumulate(h1->values.begin(), h1->values.end(), 0, addMaps);
    //sumCoef2 = std::accumulate(h2->values.begin(), h2->values.end(), 0, addMaps);
    
    std::map< std::pair<int,int>, unsigned long int >::iterator it1=h1->values.begin();
    std::map< std::pair<int,int>, unsigned long int >::iterator it2=h2->values.begin();

    while(true){
        double v1 = (double)it1->second/(double)h1->totalSum;
        double v2 = (double)it2->second/(double)h2->totalSum;

        if (it1->first == it2->first){
            similarity += ((v1 > v2)? (v2/v1)*(v1-v2) : (v1/v2)*(v2-v1));
            ++it1;
            ++it2;
        } else if (it1->first < it2->first){
            similarity += 1;
            ++it1;
        } else {
            similarity += 1;
            ++it2;
        }

        if (it1 == h1->values.end()){
            while (it2 != h2->values.end()){
                similarity += 1;
                ++it2;
            }
            break;
        }
        if (it2 == h2->values.end()){
            while (it1 != h1->values.end()){
                similarity += 1;
                ++it1;
            }
            break;
        }
    }

    return (similarity);
}

double proportionalSimilarity(Hist2D *h1, Hist2D *h2){
    double similarity = 0.0;
    unsigned long long int sumCoef1 = 0, sumCoef2 = 0;

    //sumCoef1 = std::accumulate(h1->values.begin(), h1->values.end(), 0, addMaps);
    //sumCoef2 = std::accumulate(h2->values.begin(), h2->values.end(), 0, addMaps);
    
    std::map< std::pair<int,int>, unsigned long int >::iterator it1=h1->values.begin();
    std::map< std::pair<int,int>, unsigned long int >::iterator it2=h2->values.begin();

    while(true){
        double v1 = (double)it1->second/(double)h1->totalSum;
        double v2 = (double)it2->second/(double)h2->totalSum;

        if (it1->first == it2->first){
            similarity += 1.0f - ((v1 > v2)? v2/v1 : v1/v2);
            ++it1;
            ++it2;
        } else if (it1->first < it2->first){
            similarity += 1.0f;
            ++it1;
        } else {
            similarity += 1.0f;
            ++it2;
        }

        if (it1 == h1->values.end()){
            while (it2 != h2->values.end()){
                similarity += 1.0f;
                ++it2;
            }
            break;
        }
        if (it2 == h2->values.end()){
            while (it1 != h1->values.end()){
                similarity += 1.0f;
                ++it1;
            }
            break;
        }
    }

    return (similarity);
}

double computeSimilarity2DByDC(Hist2D *h1, Hist2D *h2){
    //double similarity = 0.0;
    unsigned long long int sumCoef1 = 0, sumCoef2 = 0;
    std::map<int,long long unsigned> sumByDC;
    std::vector<double> similarity;

    //sumCoef1 = std::accumulate(h1->values.begin(), h1->values.end(), 0, addMapsByDC);
    //sumCoef2 = std::accumulate(h2->values.begin(), h2->values.end(), 0, addMapsByDC);

    std::map< std::pair<int,int>, unsigned long int >::iterator it1=h1->values.begin();
    std::map< std::pair<int,int>, unsigned long int >::iterator it2=h2->values.begin();

    while(true){

        if (it1->first == it2->first){
            selectedDC = it1->first.second;
            if (similarity.size() < selectedDC+1) similarity.resize(selectedDC+1,0.0f);
            similarity[selectedDC] += fabs(it1->second/(double)(h1->sumByDc[selectedDC]) - it2->second/(double)(h2->sumByDc[selectedDC]));
            ++it1;
            ++it2;
        } else if (it1->first < it2->first){
            selectedDC = it1->first.second;
            if (similarity.size() < selectedDC+1) similarity.resize(selectedDC+1,0.0f);
            //sumCoef1 = std::accumulate(h1->values.begin(), h1->values.end(), 0, addMapsByDC);
            similarity[selectedDC] += (it1->second/(double)(h1->sumByDc[selectedDC]));
            ++it1;
        } else {
            selectedDC = it2->first.second;
            if (similarity.size() < selectedDC+1) similarity.resize(selectedDC+1,0.0f);
            //sumCoef2 = std::accumulate(h2->values.begin(), h2->values.end(), 0, addMapsByDC);
            similarity[selectedDC] += (it2->second/(double)h2->sumByDc[selectedDC]);
            ++it2;
        }

        if (it1 == h1->values.end()){
            while (it2 != h2->values.end()){
                selectedDC = it2->first.second;
                if (similarity.size() < selectedDC+1) similarity.resize(selectedDC+1,0.0f);
                //sumCoef2 = std::accumulate(h2->values.begin(), h2->values.end(), 0, addMapsByDC);
                similarity[selectedDC] += it2->second/(double)(h2->sumByDc[selectedDC]);
                ++it2;
            }
            break;
        }
        if (it2 == h2->values.end()){
            while (it1 != h1->values.end()){
                selectedDC = it1->first.second;
                if (similarity.size() < selectedDC+1) similarity.resize(selectedDC+1,0.0f);
                //sumCoef1 = std::accumulate(h1->values.begin(), h1->values.end(), 0, addMapsByDC);
                similarity[selectedDC] += it1->second/(double)(h1->sumByDc[selectedDC]);
                ++it1;
            }
            break;
        }
    }

    //for (unsigned int i=0; i< similarity.size(); i++){
    //    cout << similarity[i] << endl;
    //}
    //std::cout <<  similarity[2] << endl;

    return (std::accumulate(similarity.begin(), similarity.end(), 0.0));
}

double computeSimilarity2D(Hist2D *h1, Hist2D *h2){
    double similarity = 0.0;
    unsigned long long int sumCoef1 = 0, sumCoef2 = 0;

    sumCoef1 = h1->totalSum; //std::accumulate(h1->values.begin(), h1->values.end(), 0, addMaps);
    sumCoef2 = h2->totalSum; //std::accumulate(h2->values.begin(), h2->values.end(), 0, addMaps);

    std::map< std::pair<int,int>, unsigned long int >::iterator it1=h1->values.begin();
    std::map< std::pair<int,int>, unsigned long int >::iterator it2=h2->values.begin();

    while(true){
        if (it1->first == it2->first){
            similarity += fabs(it1->second/(double)sumCoef1 - it2->second/(double)sumCoef2);
            ++it1;
            ++it2;
        } else if (it1->first < it2->first){
            similarity += (it1->second/(double)sumCoef1);
            ++it1;
        } else {
            similarity += (it2->second/(double)sumCoef2);
            ++it2;
        }

        if (it1 == h1->values.end()){
            while (it2 != h2->values.end()){
                similarity += it2->second/(double)sumCoef2;
                ++it2;
            }
            break;
        }
        if (it2 == h2->values.end()){
            while (it1 != h1->values.end()){
                similarity += it1->second/(double)sumCoef1;
                ++it1;
            }
            break;
        }
    }

    return (similarity);
}


void openFiles(std::string filepath){
    
    filelist = searchFiles(filepath);
    
    char var_name[50], label2[50];    
    unsigned long long int Nrays;
    int dc,thickness;
    unsigned long int value;
    int v=0;
    bool success;
    
    for (unsigned int i=0; i<filelist.size(); i++){
        //std::cout << "file:" << filepath+filelist[i] << std::endl;
        FILE* pFile = fopen( (filepath+filelist[i]).c_str(), "r+" );
        if (pFile){
            success = true;     
            // new histogram
            Hist2D* h = new Hist2D();
            char temp;
            h->name = filelist[i];
            v = sscanf(filelist[i].c_str(), "%c%d", &temp, &h->index );           
          
            v = fscanf(pFile, "%s = [ ", var_name);
            char c = ',';
            while ( c==','){            
                v = fscanf(pFile,"[%d,%d,%lu]%c", &thickness, &dc, &value, &c);
                std::pair<int,int> p = std::make_pair(thickness, dc);
                h->values.insert( std::pair< std::pair<int,int>, unsigned long> (p,value));
                //std::cout << value << std::endl;
                if (h->sumByDc.size() < dc+1) h->sumByDc.resize(dc+1,0);
                h->sumByDc[dc] += value;
                if (v<4){
                    success=false;
                    break;
                }
            }  
            h->totalSum = std::accumulate(h->sumByDc.begin(), h->sumByDc.end(), 0); 
            if (success){  
                //printf("loaded %d hist: %s\n", i,var_name);
                
                histlist2D[h->index] = h;

            } else {
                //printf("error loading %s\n", var_name);
            }
            fclose(pFile);
        } else {
             std::cerr << "Cannot open " << filepath+filelist[i] << std::endl;
        }
    }
    
    std::clog << "loaded " << histlist2D.size() << " histograms" << std::endl; 
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




