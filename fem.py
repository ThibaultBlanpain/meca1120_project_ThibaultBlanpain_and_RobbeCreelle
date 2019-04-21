# -------------------------------------------------------------------------
#
# PYTHON for FEM DUMMIES 18-19
# FEM : Small library for teaching purpose
#
#  Copyright (C) 2019 UCL-EPL : Vincent Legat
#  All rights reserved.
#
#  Vincent Legat
#
# -------------------------------------------------------------------------
#


import numpy as np

# -------------------------------------------------------------------------

_gaussTri3Xsi    = np.array([0.166666666666667,0.666666666666667,0.166666666666667])
_gaussTri3Eta    = np.array([0.166666666666667,0.166666666666667,0.666666666666667])
_gaussTri3Weight = np.array([0.166666666666667,0.166666666666667,0.166666666666667])

_gaussEdg2Xsi    = np.array([-0.5773502691896257, 0.5773502691896257])
_gaussEdg2Weight = np.array([1.0,1.0])

class IntegrationRule(object):

  def __init__(self,elementType,n):
    if (elementType == "Triangle" and n == 3) :
      self.name = "Gauss with 3 points"
      self.n = n
      self.xsi    = _gaussTri3Xsi
      self.eta    = _gaussTri3Eta
      self.weight = _gaussTri3Weight
    elif (elementType == "Edge" and n == 2) :
      self.name = "Gauss with 2 points"
      self.n = n
      self.xsi    = _gaussEdg2Xsi
      self.weight = _gaussEdg2Weight
    else :
      self.name = "Unknown rule"
      self.n = 0

  def printf(self):
    print(" Integration rule : %s " % self.name)
    print(" Number of nodes = %d " % self.n)
    print(" xsi     = ",self.xsi)
    print(" eta     = ",self.eta)
    print(" weights = ",self.weight)

# -------------------------------------------------------------------------

class LinearSystem(object):

  def __init__(self,n):
    self.A = np.zeros((n,n))
    self.B = np.zeros(n)
    self.n = n

  def printf(self):
    for i in range(self.n):
      for j in range(self.n):
        if (self.A[i,j] == 0) :
          print("         ", end='')
        else :
          print(" %+.1e" % self.A[i,j],end='')
      print(" :  %+.1e" % self.B[i]);

  def printff(self):
    for i in range(self.n):
      for j in range(self.n):
        if (self.A[i,j] == 0) :
          print("     ", end='')
        else :
          print(" %+.1f" % self.A[i,j],end='')
      print(" :  %+.1f" % self.B[i]);



  def reset(self):
    self.A = np.zeros((self.n,self.n))
    self.B = np.zeros(self.n)


# ---- Eliminate in place the system

  def eliminate(self):
    A = self.A; B = self.B; n = self.n

    # Gauss elimination

    for k in range(n):
      if (abs(A[k,k]) <= 1e-16 ) :
        print("Pivot index %d  " % k)
        print("Pivot value %e  " % A[k,k])
        raise ValueError("Cannot eliminate with such a pivot")
      for i in range(k+1,n):
        factor = A[i,k] / A[k,k]
        for j in range(k+1,n):
          A[i,j] = A[i,j] - A[k,j] * factor
        B[i] = B[i] - B[k] * factor

    # Backward substitution

    for i in range(n-1,-1,-1):
      factor = 0;
      for j in range(i+1,n):
        factor += A[i,j] * B[j];
      B[i] = ( B[i] - factor)/A[i,i];

# ---- Constrain a node to a given value

  def constrain(self,myNode,myValue):
    A = self.A; B = self.B; n = self.n
    for i in range(n):
      B[i] -= myValue * A[i,myNode]
      A[i,myNode] = 0
    for i in range(n):
      A[myNode,i] = 0
    A[myNode,myNode] = 1;
    B[myNode] = myValue;

# -------------------------------------------------------------------------

class Mesh(object):

  def __init__(self,filename):
    with open(filename,"r") as f :
      self.nNode = int(f.readline().split()[3])
      self.xy    = np.array(list(list(float(w)
        for w in f.readline().split()[2:]) for i in range(self.nNode)))
      self.nElem = int(f.readline().split()[3])
      self.elem  = np.array(list(list(int(w)
        for w in f.readline().split()[2:]) for i in range(self.nElem)))
      self.X     = self.xy[:,0]
      self.Y     = self.xy[:,1]

  def printf(self):
    print("Number of nodes %d" % self.nNode)
    for i in range(self.nNode):
      print("%6d : %14.7e %14.7e" % (i,*self.xy[i,:]))
    print("Number of triangles %d" % self.nElem)
    for i in range(self.nElem):
      print("%6d : %6d %6d %6d" % (i,*self.elem[i,:]))

  def write(self,filename):
    with open(filename,"w") as f :
      f.write("Number of nodes %d\n" % self.nNode)
      for i in range(self.nNode):
        f.write("%6d : %14.7e %14.7e\n" % (i,*self.xy[i,:]))
      f.write("Number of triangles %d\n" % self.nElem)
      for i in range(self.nElem):
        f.write("%6d : %6d %6d %6d\n" % (i,*self.elem[i,:]))

  def writeShift(self,filename):
    with open(filename,"w") as f :
      f.write("Number of nodes %d\n" % self.nNode)
      self.xy = self.xy - 1.0
      for i in range(self.nNode):
        f.write("%6d : %14.7e %14.7e\n" % (i,*self.xy[i,:]))
      f.write("Number of triangles %d\n" % self.nElem)
      for i in range(self.nElem):
        f.write("%6d : %6d %6d %6d\n" % (i,*self.elem[i,:]))



# -------------------------------------------------------------------------

class Edges(object):

  def __init__(self,mesh):
    self.mesh = mesh
    self.nEdges = mesh.nElem * 3
    self.nBoundary = 0
    self.edges = [[0 for i in range(4)] for i in range(self.nEdges)]
    for i in range (mesh.nElem) :
      for j in range(3) :
        id = i*3 + j
        self.edges[id][0] = mesh.elem[i][j]
        self.edges[id][1] = mesh.elem[i][(j+1)%3]
        self.edges[id][2] = i
        self.edges[id][3] = -1
    self.edges.sort(key = lambda item : -(min(item[0:2])*self.nEdges)-max(item[0:2]))
    index = 0
    for i in range(self.nEdges) :
      if (self.edges[i][0:2] != self.edges[i-1][1::-1]) :
         self.edges[index] = self.edges[i]
         index += 1
      else :
         self.edges[index-1][3] = self.edges[i][2]
    del self.edges[index:]
    self.edges.sort(key = lambda item : item[3])
    self.nBoundary = 2*index - self.nEdges
    self.nEdges = index

  def printf(self):
    print("Number of edges %d" % self.nEdges)
    print("Number of boundary edges %d" % self.nBoundary)
    for i in range(self.nEdges):
      print("%6d : %4d %4d : %4d %4d" % (i,*self.edges[i]))


# -------------------------------------------------------------------------

class Grid(object):

  def __init__(self,filename):
    with open(filename,"r") as f :
      self.nx     = int(f.readline().split()[1])
      self.ny     = int(f.readline().split()[1])
      self.Ox     = float(f.readline().split()[1])
      self.Oy     = float(f.readline().split()[1])
      self.dx     = float(f.readline().split()[1])
      self.nodata = int(f.readline().split()[1])
      self.elem   = np.array(list(list(float(w)
        for w in f.readline().split()[:]) for i in range(self.ny)))
      self.elem   = np.flipud(self.elem)

  def printf(self):
    print("nx = %d" % self.nx)
    print("ny = %d" % self.ny)
    print("Ox = %.7e" % self.Ox)
    print("Oy = %.7e" % self.Oy)
    print("dx = %.7e" % self.dx)
    print("Size of the table : ",np.shape(self.elem))

# -------------------------------------------------------------------------



class Advection(object):

  def __init__(self, filename):
    self.mesh = Mesh(filename)
    self.edges = Edges(self.mesh)
    self.rule = IntegrationRule("Triangle", 3)
    self.ruleEdge = IntegrationRule("Edge",2)

    size = 3 * self.mesh.nElem
    self.size = size
    self.C = np.zeros(size)
    self.U = np.zeros(size)
    self.V = np.zeros(size)
    self.F = np.zeros(size)

    self.mapEdgeLeft = np.zeros((self.edges.nEdges,2),dtype=np.int)
    self.mapEdgeRight = np.zeros((self.edges.nEdges,2),dtype=np.int)
    for iEdge in range(self.edges.nBoundary,self.edges.nEdges):
      myEdge = self.edges.edges[iEdge]
      elementLeft  = myEdge[2]
      elementRight = myEdge[3]
      nodesLeft    = self.mesh.elem[elementLeft]
      nodesRight   = self.mesh.elem[elementRight]
      self.mapEdgeLeft[iEdge,:]  = [3*elementLeft  + np.nonzero(nodesLeft  == myEdge[j])[0][0] for j in range(2)]
      self.mapEdgeRight[iEdge,:] = [3*elementRight + np.nonzero(nodesRight == myEdge[j])[0][0] for j in range(2)]

    self.mapTriangle = np.zeros((self.mesh.nElem,3),dtype=np.int)
    for iElem in range(self.mesh.nElem) :
      self.mapTriangle[iElem,:] = [3*iElem+j for j in range(3)]


  def initialize(self):
    theMesh  = self.mesh
    for iElem in range(theMesh.nElem) :
      map = [3*iElem+j for j in range(3)]
      nodes = theMesh.elem[iElem]
      x = theMesh.X[nodes]
      y = theMesh.Y[nodes]
      c = np.exp(-(x-0.8)*(x-0.8)/0.01) * np.exp(-(y-0.5)*(y-0.5)/0.01)
      u,v,eta = self.stommel(x,y)
      self.C[map] = c
      self.U[map] = u
      self.V[map] = v


  def stommel(self,x,y):

#
#   Solution analytique de Stommel dans un carre [0,1]x[0,1]
#   Modelisation de l'elevation de l'ocean Atlantique dans un carre adimensionnel
#   Ce modele que l'on attribue generalement au grand oceanographe Henry M.
#   Stommel (1920-1992), est considere comme le premier modele qualitativement correct du Gulf Stream
#
    tau0 = 0.1
    L = 1e6
    gamm = 1e-6
    rho = 1000
    delta = 1
    g = 9.81
    h = 1000
    f0 = 1e-4
    beta = 0.5e-11

    Y = y - 0.5
    epsilon = gamm / (L * beta)

    Z1 = (-1 + np.sqrt(1 + (2 * np.pi * delta * epsilon) * (2 * np.pi * delta * epsilon))) / (2 * epsilon)
    Z2 = (-1 - np.sqrt(1 + (2 * np.pi * delta * epsilon) * (2 * np.pi * delta * epsilon))) / (2 * epsilon)

    D = ((np.exp(Z2) - 1) * Z1 + (1 - np.exp(Z1)) * Z2) / (np.exp(Z1) - np.exp(Z2))

    f1 = np.pi / D * (1 + ((np.exp(Z2) - 1) * np.exp(x * Z1) + (1 - np.exp(Z1)) * np.exp(x * Z2)) / (np.exp(Z1) - np.exp(Z2)))
    f2 = 1 / D* (((np.exp(Z2) - 1) * Z1 * np.exp(x * Z1) + (1 - np.exp(Z1)) * Z2 * np.exp(x * Z2)) / (np.exp(Z1) - np.exp(Z2)))

    eta = (D * tau0 * f0 * L / (np.pi * gamm * rho * delta * g * h) *
          ( - gamm / (f0 * delta * np.pi) * f2 * np.sin(np.pi * Y)
          + 1 / np.pi * f1 * (np.cos(np.pi * Y) * (1 + beta * Y)
          - beta / np.pi * np.sin(np.pi * Y) ) ))

    u = D * tau0 / (np.pi * gamm * rho * h) * f1 * np.sin(np.pi * Y)
    v = D * tau0 / (np.pi * gamm * rho * delta * h) * f2 * np.cos(np.pi * Y)

    return u,v,eta

# -------------------------------------------------------------------------

from timeit import default_timer as timer

def tic(message = ''):
  global startTime
  startTime = timer()

def toc(message = ''):
  global startTime
  stopTime = timer()
  if message:
    message = ' (' + message + ')' ;
  print("Elapsed time is %.6f seconds %s" % ((stopTime - startTime),message) )
  elapsedTime = stopTime - startTime;
  startTime = timer()
  return elapsedTime

# -------------------------------------------------------------------------
