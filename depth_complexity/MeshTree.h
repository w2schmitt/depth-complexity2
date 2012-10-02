/* 
 * File:   MeshTree.h
 * Author: sombra
 *
 * Created on 20 de Setembro de 2012, 11:18
 */

#ifndef MESHTREE_H
#define	MESHTREE_H

#include <iostream>
#include <algorithm>
#include <istream>
#include <cmath>
#include <string>
#include <sstream>
#include <cstdlib>
#include <vector>
#include <assert.h>
#include <fstream>
#include <map>
#include <list>
#include <set>
#include <boost/lexical_cast.hpp>

#include "util.h"
#include "dc_3d.h"
#include "dc_3d_random.h"
#include "bps_model.h"

typedef struct{
  int nodeCont;
  std::vector<string> rankList;
  std::string edgeList,
              nodeList; 
}WriteGraphAux;

class MeshTree {
public:
  MeshTree();
  MeshTree(bool useNormalDepth);
  virtual ~MeshTree();
  
  void MakeTree(TriMesh &mesh, int minDepth);
  void MakeCompleteTree(TriMesh &mesh, int maxHeight);
  string WriteTreeGraph(std::ostream &file);
  string TreeGraph(int height, WriteGraphAux& graph);
  
private:
  int depthComplexity;
  MeshTree *front, *back;
  TriMesh* mesh;
  static DepthComplexity3D *dc_3d;
  static RDepthComplexity3D *dcr_3d;
  static bool normalDepth;
};

#endif	/* MESHTREE_H */

