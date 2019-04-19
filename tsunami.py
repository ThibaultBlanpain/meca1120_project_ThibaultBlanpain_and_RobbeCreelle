
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
<<<<<<< HEAD
theMeshFile = "/Users/thibaultblanpain/Documents/GitHub/meca1120_project_ThibaultBlanpain_and_RobbeCreelle/PacificTriangleTiny.txt"
=======
theMeshFile = "PacificTriangleTiny.txt"
>>>>>>> 131f326a102d405f2f0723af5980b9a35d45380e
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
def bathymetrie(theMeshFile,ielem,xsi,eta):# renvoie une interpollation de la bathymetrie au point xsi eta
    theMesh=readMesh(theMeshFile)
    Nodes=elem[ielem]
    return interpollation2D(H[Nodes],xsi,eta)

def computeShapeTriangle(theMesh,theElement) :
  dphidxsi = np.array([ 1.0, 0.0,-1.0])
  dphideta = np.array([ 0.0, 1.0,-1.0])
  nodes = theMesh.elem[theElement]
  x = theMesh.X[nodes]
  y = theMesh.Y[nodes]
  dxdxsi = x @ dphidxsi
  dxdeta = x @ dphideta
  dydxsi = y @ dphidxsi
  dydeta = y @ dphideta
  jac = abs(dxdxsi*dydeta - dxdeta*dydxsi)
  dphidx = (dphidxsi * dydeta - dphideta * dydxsi) / jac
  dphidy = (dphideta * dxdxsi - dphidxsi * dxdeta) / jac
  return [dphidx,dphidy,jac]

def computeShapeEdge(theEdges,iEdge):
  nodes = theEdges.edges[iEdge][0:2]
  x = theEdges.mesh.X[nodes]
  y = theEdges.mesh.Y[nodes]
  dx = x[1] - x[0]
  dy = y[1] - y[0]
  jac = np.sqrt(dx*dx+dy*dy)
  nx =  dy / jac
  ny = -dx / jac
  return[nx,ny,jac]


def compute(theMeshFile,theResultFiles,U,V,E,dt,nIter,nSave):

  #calcul des matrices
  Counter=0
  while Counter!= nIter:
      #calcul des matrices
      #inversion et euler explicite
      MatU=np.zeros((len(U),3))
      MatE=MatU
      MatV=MatU
      U = U+dt*MatU #euler explicite
      V = V+dt*MatV
      E = E+dt*MatE
      Counter+=1
      #sauvegarde
      #if (Counter % nSave == 0):
          #writeResult(theResultFiles,Counter,E)

  return [U,V,E]


U = np.zeros([nElem,3])
V = np.zeros([nElem,3])
E=np.zeros([nElem,3])
theResultFiles = "eta-%06d.txt"
print(compute(theMeshFile,theResultFiles,U,V,E,1,8,1))
