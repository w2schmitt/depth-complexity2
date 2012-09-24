/* 
 * File:   bps_offline.cpp
 * Author: sombra
 *
 * Created on 6 de Agosto de 2012, 16:03
 */

#include <cstdlib>
#include "bps_model.h"
#include "flags.h"

std::vector<Segment> loadOFFLines(std::istream& in);
void writeOFFModel(const TriMesh& model, std::ostream& out);
std::string getExtension(const std::string& filename);

/*
 * 
 */
int main(int argc, char** argv) {
  cmd_usage("Program that divides a model, based in the lines generated by depthcomplexity3d_offline");
  
  const char *filename = cmd_option("-f", "", "Model in the OBJ or OFF format");
  const char *modelRays = cmd_option("-r", "", "OFF File with maximum intersections lines");
  
  const char *modelA = cmd_option("-fa", "", "OFF file with the first part of the model");
  const char *modelB = cmd_option("-fb", "", "OFF file with the other part of the model");
  
  try{
    if(strcmp(filename, "")==0)
      throw "No model file specified";
    
    //carrega modelo a ser dividido
    tic();
    std::ifstream file(filename);
    std::string ext = getExtension(filename);
    TriMesh mesh;
    if (ext == "off" || ext == "OFF")
      mesh = loadOFFMesh(file);
    else if (ext == "obj" || ext == "OBJ")
      mesh = loadOBJMesh(file);
    else
      throw "Unknown file type for model!";
    file.close();
    toc("Loading Mesh");
    
    //Carrega linhas com o maximo de interseções
    tic();
    std::ifstream line_file(modelRays);
    std::string ext_line = getExtension(modelRays);
    std::vector<Segment> lines;
    if (ext_line == "off" || ext_line == "OFF")
      lines = loadOFFLines(line_file);
    else
      throw "Unknown file type for lines!";
    file.close();
    toc("Loading Lines");
    
    TriMesh part_a, part_b;
    
    /*vec3d planeNormal = lines.at(0).a - lines.at(0).b,
          model_center = ((mesh.aabb.max + mesh.aabb.min)/2.0);
    planeNormal.normalize();
    cout << model_center << std::endl;*/
    vec4d plane = defineCuttingPlane(lines, mesh.aabb);//makePlane(planeNormal, model_center);
    subDivideModel(plane, mesh, part_a, part_b);
    
    //write part_a of the model
    if (strcmp(modelA, "")!=0) {
      std::string extOff = getExtension(modelA);
      if (extOff == "off" || extOff == "OFF") {
        std::ofstream fileRays(modelA);
        writeOFFModel(part_a,fileRays);
        fileRays.close();
      } else throw "Model A file should be *.off!";
    }
    
    //write part_a of the model
    if (strcmp(modelB, "")!=0) {
      std::string extOff = getExtension(modelB);
      if (extOff == "off" || extOff == "OFF") {
        std::ofstream fileRays(modelB);
        writeOFFModel(part_b,fileRays);
        fileRays.close();
      } else throw "Model B file should be *.off!";
    }
    
  } catch (const char* msg)  {
    std::cerr << "Failed: " << msg << std::endl;
    return 1;
  } catch (std::string msg) {
    std::cerr << "Failed: " << msg << std::endl;
    return 1;
  }
  
  return 0;
}


std::vector<Segment> loadOFFLines(std::istream& in){
  
  std::clog << "Loading OFF file" << std::endl;

  std::string head;
  in >> head;
  if (head != "OFF")
    throw "Does not start with OFF!";

  int nverts, nfaces, nedges;
  
  if (!(in >> nverts >> nfaces >> nedges))
    throw "Could not read number of vertices, faces, edges";

  std::vector<Segment> lines;
  lines.reserve(nfaces);
  
  std::vector<vec3d> vertices(nverts);
  for (int i=0; i<nverts; ++i) {
    in >> vertices[i].x >> vertices[i].y >> vertices[i].z;
  }
  
  for (int i=0; i<nfaces; ++i) {
    int sz, a, b;
    Segment s;
    in >> sz >> a >> b;
    s.a = vertices.at(a);
    s.b = vertices.at(b);
    lines.push_back(s);
  }
  
  std::clog << "Loaded " << lines.size() << " lines" << std::endl;
  
  return lines;
}

void writeOFFModel(const TriMesh& model, std::ostream& out){
  out << "OFF" << std::endl;
  int m_size= model.faces.size();
  
  out << m_size*3 << " " << m_size << " " << m_size*3 << std::endl;
  
  std::vector<Triangle>::const_iterator ite = model.faces.begin();
  std::vector<Triangle>::const_iterator end = model.faces.end();
  for(;ite != end; ++ite){
    out << ite->a.x << " " << ite->a.y << " " << ite->a.z << std::endl;
    out << ite->b.x << " " << ite->b.y << " " << ite->b.z << std::endl;
    out << ite->c.x << " " << ite->c.y << " " << ite->c.z << std::endl;
  }
  
  for(int i= 0; i < m_size; i++){
    out << "3 " << 3*i << " " << 3*i+1 << " " << 3*i+2 << std::endl;
  }
}

std::string getExtension(const std::string& filename)
{
    std::string::size_type dotpos = filename.rfind(".");
    if (dotpos != std::string::npos)
        return filename.substr(dotpos+1);
    return "";
}