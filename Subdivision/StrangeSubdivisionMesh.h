#ifndef _strange_dubdivmesh_
#define _strange_dubdivmesh_

#include "AdaptiveLoopSubdivisionMesh.h"

class StrangeSubdivisionMesh : public AdaptiveLoopSubdivisionMesh
{
public:
  virtual void Subdivide() {
    // ....
    AdaptiveLoopSubdivisionMesh::Subdivide();
  }

protected:
  bool Subdividable(unsigned int fi){
    // Every 4th face is not subdividable - kinda strange!
    // Do something more interesting...
    return (fi % 4);
  }

};

#endif
