
#include <stdexcept>

#include "MeshTree.h"

DepthComplexity3D *MeshTree::dc_3d = NULL;
RDepthComplexity3D *MeshTree::dcr_3d = NULL;
bool MeshTree::normalDepth = false;

MeshTree::MeshTree(){
  depthComplexity = 0;
  front = NULL;
  back = NULL;
  
  if(MeshTree::dc_3d != NULL || MeshTree::dcr_3d != NULL){
    normalDepth = false;
    MeshTree::dcr_3d = new RDepthComplexity3D(512, 512, 10);
    MeshTree::dcr_3d->setComputeGoodRays(false);
    MeshTree::dcr_3d->setComputeHistogram(false);
    MeshTree::dcr_3d->setComputeMaximumRays(true);
    MeshTree::dc_3d = NULL;
  }
}

MeshTree::MeshTree(bool useNormalDepth){
  
  MeshTree::normalDepth = useNormalDepth;
  
  if(MeshTree::normalDepth){
    MeshTree::dc_3d = new DepthComplexity3D(512, 512, 10);
    MeshTree::dc_3d->setComputeGoodRays(false);
    MeshTree::dc_3d->setComputeHistogram(false);
    MeshTree::dc_3d->setComputeMaximumRays(true);
    MeshTree::dcr_3d = NULL;
  }else{
    MeshTree::dcr_3d = new RDepthComplexity3D(512, 512, 10);
    MeshTree::dcr_3d->setComputeGoodRays(false);
    MeshTree::dcr_3d->setComputeHistogram(false);
    MeshTree::dcr_3d->setComputeMaximumRays(true);
    MeshTree::dc_3d = NULL;
  }
  
  depthComplexity = 0;
  front = NULL;
  back = NULL;
}

MeshTree::~MeshTree(){
  delete front;
  delete back;
  delete MeshTree::dc_3d;
  delete MeshTree::dcr_3d;
  delete mesh;
}

void MeshTree::MakeTree(TriMesh &mesh, int minDepth){
  std::vector<Segment> segList;
  if(MeshTree::normalDepth){
    MeshTree::dc_3d->process(mesh);
    segList.assign(MeshTree::dc_3d->maximumRays().begin(), MeshTree::dc_3d->maximumRays().end());
    depthComplexity = MeshTree::dc_3d->maximum();
    
  }else{
    MeshTree::dcr_3d->process(mesh);
    segList.assign(MeshTree::dcr_3d->maximumRays().begin(), MeshTree::dcr_3d->maximumRays().end());
    depthComplexity = MeshTree::dcr_3d->maximum();
  }
  
  if(depthComplexity < minDepth){
    this->mesh = new TriMesh;
    this->mesh->faces.assign(mesh.faces.begin(), mesh.faces.end());
    this->mesh->aabb = mesh.aabb;
    return;
  }  
  
  vec4d cutPlane = defineCuttingPlane(segList, mesh.aabb);
  segList.clear();
  TriMesh meshFront, meshBack;
  
  subDivideModel(cutPlane, mesh, meshFront, meshBack);
  
  cout << meshFront.faces.size() << " " << cutPlane << endl;
  front = new MeshTree();
  front->MakeTree(meshFront, minDepth);
  
  back = new MeshTree();
  back->MakeTree(meshBack, minDepth);
  
  mesh.faces.clear();
}

void MeshTree::MakeCompleteTree(TriMesh &mesh, int maxHeight){
  
  if(maxHeight == 0){
    this->mesh = new TriMesh;
    this->mesh->faces.assign(mesh.faces.begin(), mesh.faces.end());
    this->mesh->aabb = mesh.aabb;
    return;
  }
  
  std::vector<Segment> segList;
  if(MeshTree::normalDepth){
    MeshTree::dc_3d->process(mesh);
    segList.assign(MeshTree::dc_3d->maximumRays().begin(), MeshTree::dc_3d->maximumRays().end());
    depthComplexity = MeshTree::dc_3d->maximum();
    
  }else{
    MeshTree::dcr_3d->process(mesh);
    segList.assign(MeshTree::dcr_3d->maximumRays().begin(), MeshTree::dcr_3d->maximumRays().end());
    depthComplexity = MeshTree::dcr_3d->maximum();
  }
  
  
  vec4d cutPlane = defineCuttingPlane(segList, mesh.aabb);
  segList.clear();
  TriMesh meshFront, meshBack;
  
  subDivideModel(cutPlane, mesh, meshFront, meshBack);
  
  cout << meshFront.faces.size() << " " << cutPlane << endl;
  front = new MeshTree();
  front->MakeCompleteTree(meshFront, maxHeight-1);
  
  back = new MeshTree();
  back->MakeCompleteTree(meshBack, maxHeight-1);
  mesh.faces.clear();
  
}

string MeshTree::WriteTreeGraph(std::ostream &file){
  WriteGraphAux graph;
  graph.edgeList = "";
  graph.nodeCont = 0;
  graph.nodeList = "";
  graph.rankList.push_back("");
  
  TreeGraph(0, graph);
  
  file << "digraph G{\n";
  file << graph.nodeList;
  
  std::vector<string>::iterator ite;
  for(ite=graph.rankList.begin(); ite != graph.rankList.end(); ite++){
    file << "{ rank = same;" << *ite << "};\n";
  }
  
  file << graph.edgeList << "}";
  
  graph.nodeList.clear();
  graph.edgeList.clear();
  graph.rankList.clear();
  
  return "done";
}

string MeshTree::TreeGraph(int height, WriteGraphAux& graph){
  string nodeName = "node"+boost::lexical_cast<string>(graph.nodeCont++);
  graph.nodeList += nodeName+"[label = \""+boost::lexical_cast<string>(depthComplexity)+"\"];\n";
  
  try{
    graph.rankList.at(height) += "\""+nodeName+"\";";
  }catch(out_of_range){
    graph.rankList.push_back("\""+nodeName+"\";");
  }
  
  string nFront, nBack;
  
  if(front != NULL){
    nFront = front->TreeGraph(height+1, graph);
    graph.edgeList += "\""+nodeName+"\":sw -> \""+nFront+"\":n;\n";
  }
  if(back != NULL){
    nBack = back->TreeGraph(height+1, graph);
    graph.edgeList += "\""+nodeName+"\":se -> \""+nBack+"\":n;\n";
  }
  return nodeName;
}
