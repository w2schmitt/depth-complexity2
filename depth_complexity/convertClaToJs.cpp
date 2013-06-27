
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <dirent.h>
#include <string.h>
#include <vector>
#include <algorithm>
#include <math.h>
#include <numeric>


void parseFile();


int main(int argc, char** argv){
        
    parseFile();
}

void parseFile(){
    
    std::string temp1;
    int temp2, NoClasses, NoModels, NoPerClass, currentModel;
    std::string class1,class2;
    
    printf("var m_class = {");
    std::cin >> temp1 >> temp2;
    std::cin >> NoClasses >> NoModels;
    char comma2 = ' ';
    
    for (int i=0; i<NoClasses; i++){        
        std::cin >> class1 >> class2 >> NoPerClass;
        if (NoPerClass==0) continue;
        printf("%c\"%s-%s\":[",comma2, class1.c_str(), class2.c_str());
        char comma=' ';
        for (int j=0; j<NoPerClass; j++){
            std::cin >> currentModel;
            printf("%c%d",comma, currentModel);
            comma=',';
        }
        comma2=',';
        printf("]");
    }
    printf("};");    
}
