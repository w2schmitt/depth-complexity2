
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <dirent.h>
#include <string.h>
#include <vector>
#include <algorithm>
#include <math.h>
#include <numeric>
#include <map>


void parseFile(std::istream& out, std::map<std::string, std::vector<int>* > &classes);
void writeClaFormat(std::map<std::string, std::vector<int>* > &classes);

int main(int argc, char** argv){
    std::map<std::string, std::vector<int>* > classes;

    for (unsigned i=1; i<argc; i++){
        //std::cout << argv[i] << std::endl;
        std::ifstream out(argv[i]);
        if (out.is_open()) {
            parseFile(out, classes);
            out.close();
        }        
    }   
    writeClaFormat(classes);
}

void writeClaFormat(std::map<std::string, std::vector<int>* > &classes){
  std::map<std::string, std::vector<int>* >::iterator it = classes.begin();
  std::map<std::string, std::vector<int>* >::iterator itend = classes.end();
  
  char comma=' ';
  printf("var m_class = {");
  for (;it!=itend; ++it){
    std::string className = it->first;
    std::vector<int>* modelsList = it->second;

    char comma2=' ';
    printf("%c\"%s\":[",comma, className.c_str());
    for (unsigned i=0; i<modelsList->size(); i++){
        printf("%c%d", comma2, modelsList->at(i) );
        comma2=',';
    }
    comma=',';
    printf("]");
  }  
  printf("};");  
}

void parseFile(std::istream& out, std::map<std::string, std::vector<int>* > &classes){
    std::string temp1;
    int temp2, NoClasses, NoModels, NoPerClass, currentModel;
    std::string class1,class2;    
    
    out >> temp1 >> temp2;
    out >> NoClasses >> NoModels;
    
    for (int i=0; i<NoClasses; i++){        
        out >> class1 >> class2 >> NoPerClass;
        std::string key = (class2.compare("0")==0? class1 : class2);

        if (NoPerClass==0) continue;

        std::vector<int>* list;
        std::map<std::string, std::vector<int>* >::iterator it = classes.find(key);
        if (it != classes.end() ){
            list = it->second;
        } else {
            list = new std::vector<int>();
            classes.insert( std::pair<std::string, std::vector<int>* >(key, list));
        }
        
        for (int j=0; j<NoPerClass; j++){
            out >> currentModel;
            list->push_back(currentModel);
        }
    }
}
