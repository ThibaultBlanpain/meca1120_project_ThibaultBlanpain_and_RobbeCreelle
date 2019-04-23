
import numpy as np

# -------------------------------------------------------------------------

def readMesh(fileName) :
  with open(fileName,"r") as f :
    nNode = int(f.readline().split()[3])
    xyz   = np.array(list(list(float(w) for w in f.readline().split()[2:]) for i in range(nNode)))
    nElem = int(f.readline().split()[3])
    elem  = np.array(list(list(int(w)   for w in f.readline().split()[2:]) for i in range(nElem)))
  X = xyz[:,0]
  Y = xyz[:,1]
  H = xyz[:,2]
  return [nNode,X,Y,H,nElem,elem]

theMeshFile = "PacificTriangleTiny.txt"
[nNode,X,Y,H,nElem,elem] = readMesh(theMeshFile)

# -------------------------------------------------------------------------

def readResult(fileBaseName,iter,nElem) :
  fileName = fileBaseName % iter
  with open(fileName,"r") as f :
    nSize = int(f.readline().split()[3])
    if (nElem != nSize) :
      print(" ==== Error : incoherent sizes : %d != %d" % (nElem,nSize))
    E = np.array(list(list(float(w) for w in f.readline().split()[2:5]) for i in range(nElem)))
    print(" === iteration %6d : reading %s ===" % (iter,fileName))
  return E

# -------------------------------------------------------------------------

def writeResult(fileBaseName,iter,E) :
  fileName = fileBaseName % iter
  nElem = E.shape[0]
  with open(fileName,"w") as f :
    f.write("Number of elements %d\n" % nElem)
    for i in range(nElem):
      f.write("%6d : %14.7e %14.7e %14.7e\n" % (i,*E[i,:]))
    print(" === iteration %6d : writing %s ===" % (iter,fileName))

# -------------------------------------------------------------------------

def initialConditionOkada(x,y) :
  R = 6371220;
  x3d = 4*R*R*x / (4*R*R + x*x + y*y);
  y3d = 4*R*R*y / (4*R*R + x*x + y*y);
  z3d = R*(4*R*R - x*x - y*y) / (4*R*R + x*x + y*y);
  lat = np.arcsin(z3d/R)*180/np.pi;
  lon = np.arctan2(y3d,x3d)*180/np.pi;
  lonMin = 142;
  lonMax = 143.75;
  latMin = 35.9;
  latMax = 39.5;
  olon = (lonMin+lonMax)/2;
  olat = (latMin+latMax)/2;
  angle = -12.95*np.pi/180;
  lon2 = olon + (lon-olon)*np.cos(angle) + (lat-olat)*np.sin(angle);
  lat2 = olat - (lon-olon)*np.sin(angle) + (lat-olat)*np.cos(angle);
  return np.all([lon2 <= lonMax,lon2 >= lonMin,lat2 >= latMin,lat2 <= latMax],axis=0).astype(int)

# -------------------------------------------------------------------------
def interpollation2D(U,xsi,eta):
    return U[0]*xsi+U[1]*eta+U[2]*(np.ones(3)-xsi-eta)
def interpollation1D(U,xsi):
    return 0.5*(U[0]*(1-xsi)+U[1]*(1+xsi))

#classe pour the Mesh pour eviter la repetition de readMesh(theMeshFile)
class Mesh(object):
    def __init__(self,theMeshFile):
        [self.nNode,self.X,self.Y,self.H,self.nElem,self.elem]=readMesh(theMeshFile)
theMesh=Mesh(theMeshFile)

def bathymetrie(theMesh,ielem,xsi,eta):# renvoie une interpollation de la bathymetrie au point xsi eta
    Nodes=theMesh.elem[ielem]
    return interpollation2D(theMesh.H[Nodes],xsi,eta)

def computeShapeTriangle(theMesh,Element,jaco) :#Element pas vectoriel malheureusement
  dphidxsi = np.array([ 1.0, 0.0,-1.0])
  dphideta = np.array([ 0.0, 1.0,-1.0])
  nodes = theMesh.elem[Element]
  x = theMesh.X[nodes]
  y = theMesh.Y[nodes]
  dxdxsi = x @ dphidxsi
  dxdeta = x @ dphideta
  dydxsi = y @ dphidxsi
  dydeta = y @ dphideta

  dphidx = (dphidxsi * dydeta - dphideta * dydxsi) / jaco[Element]
  dphidy = (dphideta * dxdxsi - dphidxsi * dxdeta) / jaco[Element]
  return [dphidx,dphidy]
def computeJacobian(theMesh):
    X=theMesh.X
    Y=theMesh.Y
    elem=theMesh.elem
    return abs((X[elem[:,2]]-X[elem[:,0]])*(Y[elem[:,1]]-Y[elem[:,0]])-(X[elem[:,1]]-X[elem[:,0]])*(Y[elem[:,2]]-Y[elem[:,0]]))


# strategie pour creer des edges: 1) repertorier les aretes 2) enlever les doublons
#3) classer Boundary et non Boundary
class Edges(object):
  def __init__(self,theMesh):
    self.mesh=theMesh
    self.nEdges = theMesh.nElem * 3
    tab= -1*np.ones((self.nEdges,4),dtype='int')
    for i in range(theMesh.nElem):
        tab[3*i:3*i+4,2]=i
        Tripoints=theMesh.elem[i]
        tab[3*i,0:2]=Tripoints[0:2]
        tab[3*i+1,0:2]=Tripoints[1:]
        tab[3*i+2,0:2]=Tripoints[2],Tripoints[0]
    self.edges=sorted(tab,key = lambda item: (min(item[0],item[1]),max(item[0],item[1])),reverse=True)
    decal=0
    for i in range(len(self.edges)-1):
        if self.edges[i][0]==self.edges[i+1][1] and self.edges[i+1][0]==self.edges[i][1]:
            decal+=1
            self.edges[i][3]=self.edges[i+1][2]
            self.edges[i+1]=[-2,-2,-2,-2]
    self.edges=sorted(self.edges,key = lambda item: item[3],reverse=True)
    n=self.nEdges
    self.nEdges=n-decal
    self.nBoundary=n-2*decal
    self.edges=self.edges[0:self.nEdges]
    self.edges=sorted(self.edges,key = lambda item: item[3])
  def printf(self):

      print("Number of edges %d" % self.nEdges)
      print("Number of boundary edges %d" % self.nBoundary)
      for i in range(self.nEdges):
        print("%6d : %4d %4d : %4d %4d" % (i,*self.edges[i]))
theEdges=Edges(theMesh)

def computeShapeEdge(theEdges,iEdge):#iEdge non vectoriel malheureusement

  nodes = theEdges.edges[iEdge][0:2]
  x = theEdges.mesh.X[nodes]
  y = theEdges.mesh.Y[nodes]
  dx = x[1] - x[0]
  dy = y[1] - y[0]
  jac = np.sqrt(dx*dx+dy*dy)#longueur du segment (attention au 0.5)
  nx =  dy / jac
  ny = -dx / jac
  return[nx,ny,jac]


def compute(theMeshFile,theResultFiles,U,V,E,dt,nIter,nSave):

  [nNode,X,Y,H,nElem,elem] = readMesh(theMeshFile)
  theMesh=Mesh(theMeshFile)
  theEdges=Edges(theMesh)


  Counter=0
  Ainverse = np.array([[18.0,-6.0,-6.0],[-6.0,18.0,-6.0],[-6.0,-6.0,18.0]])
  jaco=computeJacobian(theMesh)
  #initialisation des matrices B
  B_U=np.zeros((len(U),3))
  B_V=B_U
  B_E=B_U

  while Counter!= nIter:

      #integrales de surfaces
      for iElem in range(nElem):
          continue
      # integrales des edges interieures
      for iEdge in range(theEdges.nBoundary,theEdges.nEdges):
          continue
      for iEdgeB in range(theEdges.nBoundary):
          continue
      #inversion des matrices
      for iElem in range(nElem):
         
          B_U[iElem] +=Ainverse@B_U[iElem]/jaco[iElem]# nx3= 3x3 @ nx1
          B_V[iElem] +=Ainverse@B_V[iElem]/jaco[iElem]
          B_E[iElem] += Ainverse@B_E[iElem]/jaco[iElem]
      # euler explicite
      U+= dt*B_U
      V+= dt*B_V
      E+= dt*B_E
      Counter+=1
      #sauvegarde
      #if (Counter % nSave == 0):
          #writeResult(theResultFiles,Counter,E)

  return [U,V,E],Counter


U = np.zeros([nElem,3])
V = np.zeros([nElem,3])
E=np.zeros([nElem,3])
theResultFiles = "eta-%06d.txt"
print(compute(theMeshFile,theResultFiles,U,V,E,1,8,1))
