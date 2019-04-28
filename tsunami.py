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
    self.mapEdgeLeft = np.zeros((self.nEdges,2),dtype=np.int)
    self.mapEdgeRight = np.zeros((self.nEdges,2),dtype=np.int)
    for iEdge in range(self.nBoundary,self.nEdges):
      myEdge = self.edges[iEdge]
      elementLeft  = myEdge[2]
      elementRight = myEdge[3]
      nodesLeft    = self.mesh.elem[elementLeft]
      nodesRight   = self.mesh.elem[elementRight]
      self.mapEdgeLeft[iEdge,:]  = [3*elementLeft  + np.nonzero(nodesLeft  == myEdge[j])[0][0] for j in range(2)]
      self.mapEdgeRight[iEdge,:] = [3*elementRight + np.nonzero(nodesRight == myEdge[j])[0][0] for j in range(2)]
    for iEdgeB in range(self.nBoundary):
      myEdge = self.edges[iEdgeB]
      elementLeft  = myEdge[2]
      nodesLeft    = self.mesh.elem[elementLeft]
      self.mapEdgeLeft[iEdgeB,:]  = [3*elementLeft  + np.nonzero(nodesLeft  == myEdge[j])[0][0] for j in range(2)]

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
  R = 6371220
  omega= 2*np.pi/86400
  gamma=10**-7
  g=9.81

  [nNode,X,Y,H,nElem,elem] = readMesh(theMeshFile)
  theMesh=Mesh(theMeshFile)
  theEdges=Edges(theMesh)

  Counter=0
  Ainverse = np.array([[18.0,-6.0,-6.0],[-6.0,18.0,-6.0],[-6.0,-6.0,18.0]])
  jaco=computeJacobian(theMesh)
  #interation rule 2D
  xsi=np.array([0.166666666666667,0.666666666666667,0.166666666666667])
  eta=np.array([0.166666666666667,0.166666666666667,0.666666666666667])
  weight=np.array([0.166666666666667,0.166666666666667,0.166666666666667])
  #integration 1D
  xsi1D=np.array([-0.5773502691896257, 0.5773502691896257])
  weight1D=np.array([1.0,1.0])
  #fonctions de formes
  phi= np.array([1-xsi-eta,xsi,eta])
  phi1D= np.asarray([1.0-xsi1D,1.0+xsi1D])/ 2

  while Counter!= nIter:
      #initialisation des matrices B
      B_U=np.zeros((len(U),3))
      B_V=np.zeros((len(V),3))
      B_E=np.zeros((len(E),3))
      #integrales de surfaces
      for iElem in range(nElem):
          h=bathymetrie(theMesh,iElem,xsi,eta)
          u=interpollation2D(U[iElem],xsi,eta)
          v=interpollation2D(V[iElem],xsi,eta)
          x=interpollation2D(theMesh.X[elem[iElem]],xsi,eta)
          y=interpollation2D(theMesh.Y[elem[iElem]],xsi,eta)
          [dphidx,dphidy] =computeShapeTriangle(theMesh,iElem,jaco)
          bigR= (4*R**2+x**2+y**2)/(4*R**2)

          B_E[iElem] += sum( (np.outer(u*h*bigR,dphidx) + np.outer(v*h*bigR,dphidy)) * weight *jaco[iElem]) #1
          B_E[iElem] += ((h*(x*u+y*v)/(R**2))*weight)@ phi/jaco[iElem] #2

          f=2*omega*np.sin(y/bigR)

          B_U[iElem] += ((f*v-gamma*u)*weight)@ phi /jaco[iElem]#3
          n_elev=interpollation2D(E[iElem],xsi,eta)# n_elev=H-h

          B_U[iElem] += ((g*x*n_elev)*weight)@phi /jaco[iElem]*1/(2*R**2)#4
          B_U[iElem] += sum (np.outer(g*n_elev*bigR,dphidx)*weight*jaco[iElem])#5

          B_V[iElem] += ((-1*f*v-gamma*u)*weight)@ phi *jaco[iElem]#idem 6 que 3
          B_V[iElem] += ((g*y*n_elev)*weight)@phi /jaco[iElem]*1/(2*R**2)#7 meme que 4
          B_V[iElem] += sum (np.outer(g*n_elev*bigR,dphidy)*weight*jaco[iElem])#8 idem que 5
      # integrales des edges interieures
      B_U=np.ravel(B_U)
      B_V=np.ravel(B_V)
      B_E=np.ravel(B_E)
      for iEdge in range(theEdges.nBoundary,theEdges.nEdges):
          mapLeft=theEdges.mapEdgeLeft[iEdge]
          mapRight=theEdges.mapEdgeRight[iEdge]
          [nx,ny,jac]=computeShapeEdge(theEdges,iEdge)
          nx,ny= -1*nx,-1*ny # car c est une normal entrante
          x=interpollation1D(theMesh.X[theEdges.edges[iEdge][0:2]],xsi1D)
          y=interpollation1D(theMesh.Y[theEdges.edges[iEdge][0:2]],xsi1D)
          h=interpollation1D(theMesh.H[theEdges.edges[iEdge][0:2]],xsi1D)
          bigR=(4*R**2+x**2+y**2)/(4*R**2)
          #calcul de u star et de e_star
          u_l=np.array([U[(mapLeft[0]/3).astype(int)][mapLeft[0]%3],U[(mapLeft[1]/3).astype(int)][mapLeft[1]%3]])
          u_r=np.array([U[(mapRight[0]/3).astype(int)][mapRight[0]%3],U[(mapRight[1]/3).astype(int)][mapRight[1]%3]])
          e_l=np.array([E[(mapLeft[0]/3).astype(int)][mapLeft[0]%3],E[(mapLeft[1]/3).astype(int)][mapLeft[1]%3]])
          e_r=np.array([E[(mapRight[0]/3).astype(int)][mapRight[0]%3],E[(mapRight[1]/3).astype(int)][mapRight[1]%3]])
          u_star=0.5*(u_l+u_r)+np.sqrt(g/h)*0.5*(e_l-e_r)
          e_star=0.5*(e_l+e_r)+np.sqrt(h/g)*0.5*(u_l-u_r)
          #remplissage de B
          B_E[mapLeft] += (weight1D*u_star*h*bigR)@phi1D*0.5*jac
          B_E[mapRight] -= (weight1D*u_star*h*bigR)@phi1D*0.5*jac #9

          B_U[mapLeft] += nx*g*(weight1D*e_star)@phi1D*0.5*jac
          B_U[mapRight] -= nx*g*(weight1D*e_star)@phi1D*0.5*jac #10

          B_V[mapLeft] += ny*g*(weight1D*e_star)@phi1D*0.5*jac
          B_V[mapRight] -= ny*g*(weight1D*e_star)@phi1D*0.5*jac #11 idem que 10

      # integrales des edges frontieres

      for iEdgeB in range(theEdges.nBoundary):
           mapLeft=theEdges.mapEdgeLeft[iEdgeB]
           [nx,ny,jac]=computeShapeEdge(theEdges,iEdgeB)
           nx,ny= -1*nx,-1*ny # car c est une normal entrante
           x=interpollation1D(theMesh.X[theEdges.edges[iEdgeB][0:2]],xsi1D)
           y=interpollation1D(theMesh.Y[theEdges.edges[iEdgeB][0:2]],xsi1D)
           h=interpollation1D(theMesh.H[theEdges.edges[iEdgeB][0:2]],xsi1D)
           bigR=(4*R**2+x**2+y**2)/(4*R**2)
           #calcul de e_star (u_star=0)
           u_l=np.array([U[(mapLeft[0]/3).astype(int)][mapLeft[0]%3],U[(mapLeft[1]/3).astype(int)][mapLeft[1]%3]])
           e_l=np.array([E[(mapLeft[0]/3).astype(int)][mapLeft[0]%3],E[(mapLeft[1]/3).astype(int)][mapLeft[1]%3]])
           e_star= e_l+np.sqrt(h/g)*u_l
           #remplissage des matrice B #12=0 car u_star=0
           B_U[mapLeft] += nx*g*(weight1D*e_star)@phi1D*0.5*jac#13
           B_V[mapLeft] += ny*g*(weight1D*e_star)@phi1D*0.5*jac#14
      #remise en nx3
      B_U=np.reshape(B_U,(nElem,3))
      B_V=np.reshape(B_V,(nElem,3))
      B_E=np.reshape(B_E,(nElem,3))
      #inversion des matrices
      for iElem in range(nElem):
          B_U[iElem] =Ainverse@B_U[iElem]/jaco[iElem]# nx3= 3x3 @ nx1
          B_V[iElem] =Ainverse@B_V[iElem]/jaco[iElem]
          B_E[iElem] = Ainverse@B_E[iElem]/jaco[iElem]

      # euler explicite
      U+= dt*B_U
      V+= dt*B_V
      E+= dt*B_E
      Counter+=1
      #sauvegarde
      #if (Counter % nSave == 0):
          #writeResult(theResultFiles,Counter,E)
  return [U,V,E]


U = np.ones([nElem,3])
V = np.zeros([nElem,3])
E=np.zeros([nElem,3])
E[3:100,:]=1
theResultFiles = "eta-%06d.txt"
print(compute(theMeshFile,theResultFiles,U,V,E,1,1,20))
#impositions des conditions aux frontieres1: E gauche = E droite  et vitesse a gauche = -1*vitesse a droite(pour U,V)
