#! /usr/bin/env python
# -*- coding: UTF-8 -*-

from clases import *

LOGp = ['Report from %s.py'%(__name__)]
WARNINGp = ['WARNINGs from %s.py'%(__name__)]


def PZB_grid(vecs,nk=10):
   """
     Returns a list of points sampling the PBZ. There will be nk^len(vecs) points
   """
   ndim = len(vecs)
   if ndim == 1:
     coefs = vecs[0]
     points = []
     for c in coefs:
        points.append(c*vecs[0])
   elif ndim == 2:
     coefsX = np.linspace(0,1,nk)
     coefsY = np.linspace(0,1,nk)
     points = []
     for c1 in coefsX:
        for c2 in coefsY:
           points.append(c1*vecs[0]+c2*vecs[1])
   elif ndim == 3:
     coefsX = np.linspace(0,1,nk)
     coefsY = np.linspace(0,1,nk)
     coefsZ = np.linspace(0,1,nk)
     points = []
     for c1 in coefsX: # nk veces        #
        for c2 in coefsY: # nk veces     # ndim veces
           for c3 in coefsZ: # nk veces  #
              points.append(c1*vecs[0]+c2*vecs[1]+c3*vecs[2])
   return points

def dens2band(H):
   M = []
   for j in range(H.shape[1]):
      aux = []
      for im in range(H.shape[0]-j):
         aux.append(H[im,im+j])
      for _ in range(H.shape[0]-j,H.shape[1]):
         aux = [0] + aux
      M.append(aux)
   M.reverse()
   M = np.array(M)
   return M


def find_neig(at,base,bonds,n=1):
   firsts = []
   for bo in bonds:
      mat_neigh = bo[1]
      cont = 0
      for elem in mat_neigh[at]:
         if elem == n:
           firsts.append((base[cont].place,bo[0]))
         cont += 1
   return firsts


def randorotator(AT,cell,latt_vec,rotate=True,rand=True):
   if rotate:
     if rand:
       v = np.array([random.random(),random.random(),random.random()])
       Q = mapeo(random.random(),0,360)
       print 'Rotado %s deg alrededor del eje'%(Q),v
       cellp = []
       for vec in cell:
          cellp.append(rotation(vec,Q,v))
       latt_vecp = []
       for vec in latt_vec:
          latt_vecp.append(rotation(vec,Q,v))
     else:
       v = np.array([0,0,1])
       Q = 90
       print 'Rotado %s deg alrededor del eje'%(Q),v
       cellp = []
       for vec in cell:
          cellp.append(rotation(vec,Q,v))
       latt_vecp = []
       for vec in latt_vec:
          latt_vecp.append(rotation(vec,Q,v))
   else:
     v=np.array([0,0,0])
     Q=0
     cellp = []
     for vec in cell:
        cellp.append(vec)
     latt_vecp = []
     for vec in latt_vec:
        latt_vecp.append(vec)
   return AT,cellp,latt_vecp


def reciprocal(latt_vec):
   """
     Should be valid for 0D, 1D and 2D
     CAREFUL!!! NOT TESTED!!!
   """
   if len(latt_vec) == 1:
     a1 = latt_vec[0]
     # Define more direct vectors to use the same formula
     rand = random.random
     a2 = np.cross( a1,np.array([rand(),rand(),rand()]) )
     a2 = a2*np.linalg.norm(a1)/np.linalg.norm(a2)
     a3 = np.cross(a1,a2)
     a3 = a3*np.linalg.norm(a1)/np.linalg.norm(a3)
   elif len(latt_vec) == 2:
     a1 = latt_vec[0]
     a2 = latt_vec[1]
     a3 = np.cross(a1,a2)
     a3 = a3*np.linalg.norm(a1)/np.linalg.norm(a3)
   elif len(latt_vec) == 3:
     a1 = latt_vec[0]
     a2 = latt_vec[1]
     a3 = latt_vec[2]

   b1 = 2*np.pi*(np.cross(a2,a3)/(np.dot(a1,np.cross(a2,a3))))
   b2 = 2*np.pi*(np.cross(a3,a1)/(np.dot(a1,np.cross(a2,a3))))
   b3 = 2*np.pi*(np.cross(a1,a2)/(np.dot(a1,np.cross(a2,a3))))

   vecs = [b1,b2,b3]
   recip_vec = []
   for i in range(len(latt_vec)): # returns only the necessary
      recip_vec.append(vecs[i])   # reciprocal vectors
   return recip_vec


# ===================================================
#    Para recorrer todas las posibles celdas
#    en busqueda de 2os, 3os,... vecinos
#    de momento recorre todas las combinaciones
#    de vectores {-1,0,1} para cada lattice_vector
#               [Es N dimensional!!!]
# ===================================================
def vec_neig(latt_vec,n=1):
   """ Generates a N-dimensional grid (N = 1,2,3) and returns all the
   combinations of the vectors 
   """
   li = [range(-n,n+1) for _ in range(len(latt_vec))]
   if len(latt_vec) == 0:
      vecs = [np.array([0,0,0])]
   elif len(latt_vec) == 1:
      vecs = []
      for i in li[0]:
         vecs.append(i*latt_vec[0])
   elif len(latt_vec) == 2:
      mxy = np.meshgrid(*li)
      matx = mxy[0]
      vecs = []
      for i in range(matx.shape[0]):
         for j in range(matx.shape[1]):
            aux = 0.0
            for iv in range(len(mxy)):
               mat = mxy[iv]
               aux += mat[i,j]*latt_vec[iv]
            vecs.append(aux)
   elif len(latt_vec) == 3:
      mxy = np.meshgrid(*li)
      matx = mxy[0]
      vecs = []
      for i in range(matx.shape[0]):
         for j in range(matx.shape[1]):
            for k in range(matx.shape[2]):
               aux = 0.0
               for iv in range(len(mxy)):
                  mat = mxy[iv]
                  aux += mat[i,j,k]*latt_vec[iv]
               vecs.append(aux)
   return vecs

def vecinlist(vec,lista):
   """
     Checks if a vector vec is in a list lista
   """
   isin = False
   for v in lista:
      eq = vec - v
      if np.allclose(eq,np.zeros(eq.shape,dtype=float)):
        isin = True
        break
      else:
        isin = False
   return isin


def sort(lista,eps=eps):
   if isinstance(lista,list):
     lista.sort() #sorts in ascending order
     x = range(0, len(lista), 1) #gets range
     x.reverse() #reverses
     for k in x:
        if abs(lista[k] - lista[k-1]) < eps: #if the list value -1 is equal to the next,
          del(lista[k-1])     #remove it
     return lista
   else:
     orde = np.sort(lista)
     res = []
     for item in orde:
        if len(res) == 0 or item > res[-1] + eps:
          res.append(item)
     return res


def trip2coo(lista):
   """
     Returns a Coo Matrix from a given list (lista) containing
   3 lists, X, Y, data such that  A[ X[i],Y[i] ] = data[i]
   """
   X, Y, data = [], [], []
   for trip in lista:
      X.append(trip[0])
      Y.append(trip[1])
      data.append(trip[2])
   return scipy.sparse.coo_matrix((data,(X,Y)))


def list2mat(lista):
   """
     Returns a 2D matrix form a given list
   * Uses the type of the first element
   """
   mat = np.matrix(np.zeros((len(lista),lista[0].shape[1]),dtype = type(lista[0][0][0,0])))
   for i in range(len(lista)):
      fila = lista[i]
      for j in range(fila.shape[1]):
         mat[i,j] = fila[0,j]
   return mat


def symmetrize(A):
   """
     This subroutine returns the complete hemitian matrix of a given one.
   """
   zeros_diag=np.matrix(np.zeros((len(A),len(A)),dtype=complex))
   diag=np.array([])
   for i in range(len(A.diagonal())):
      diag=np.append(diag, A.diagonal()[i])
   for i in range(len(A)):
      zeros_diag[i,i]=diag[i]
#   print A + A.transpose().conjugate() - zeros_diag
#   print 'hola'
   return A + A.transpose().conjugate() - zeros_diag


def mobiusator(q,h):
   x=(1.+(h/2.)*np.sin(q/2.))*np.cos(q)
   y=(1.+(h/2.)*np.cos(q/2.))*np.sin(q)
   z=(h/2.)*np.cos(q/2.)
   return x,y,z


def tobm(m):
  """TRansforms in an upper band matrix"""
  from numpy import array
  n=len(m)
  bm=array([[0.0j for i in range(n)] for j in range(n)])
  for i in range(n):
    for j in range(i,n):
      bm[n-j+i-1][j]=m[i,j]
  return array(bm)


def read_params(archivo):
   f = open(sys.argv[1], "r")
   params=[]
   while True:
        linea = f.readline()
        if not linea: break
        linea = linea.lstrip()
        if linea[0] != '#':
          params.append(linea.rstrip('\n').split("=")[0].lstrip().rstrip())
          params.append(linea.rstrip('\n').split("=")[1].lstrip().rstrip())
   f.close()
   return params



"""~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ K-PATHS ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"""
def recorrido(lista,nkx):
   RE = []
   for ipunto in range(len(lista)-1):
      P1 = lista[ipunto]
      P2 = lista[ipunto+1]
      x=np.linspace(P1[0],P2[0],nkx+1)
      y=np.linspace(P1[1],P2[1],nkx+1)
      z=np.linspace(P1[2],P2[2],nkx+1)
      for kx,ky,kz in zip(x[1:],y[1:],z[1:]):
         RE.append((kx,ky,kz))
   return RE # [RECORRIDO_X,RECORRIDO_Y,RECORRIDO_Z]


def Kp_near(nkx):
    '''
      This functions returns 2 lists containing all the kx and ky in which
    we want to calculate the bands.
    '''
    RECORRIDO_X=[]
    RECORRIDO_Y=[]
    x=np.linspace(0,0,nkx)                       #
    y=np.linspace(2.*np.pi/3.,4.*np.pi/3.,nkx)   #  Gamma/2 --> K
    for kx,ky in zip(x,y):                       #
       RECORRIDO_X.append(kx)                    #
       RECORRIDO_Y.append(ky)                    #

    x=np.linspace(0,np.pi/np.sqrt(3.),nkx)     #
    y=np.linspace(4.*np.pi/3.,np.pi,nkx)    #  K --> K'
    for kx,ky in zip(x[1:len(x)],y[1:len(y)]): #
       RECORRIDO_X.append(kx)                  #
       RECORRIDO_Y.append(ky)                  #

    return RECORRIDO_X,RECORRIDO_Y
"""~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"""



"""~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Visualization ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"""
def structure(latt_vec, nvec, base,name='structure.xyz'):
   """
     XYZ file of all the cells taken into account
   """
   vecs_cells = vec_neig(latt_vec,nvec)
   archivo=open(name,'w')
   archivo.write(str(len(vecs_cells)*len(base))+'\n\n')
   for r in vecs_cells:
      for E in base:
         vec = E.position + r
         archivo.write(E.atom.name+'  '+str(vec[0])+'  '+str(vec[1])+'  '+str(vec[2])+'\n')
   archivo.close()



def show_base(base):
   """
     Quick and easy way to see all the properties of every element of the basis
   """
   for elemento in base:
      print 'Element',elemento.place,':'
      print '  Atom:',elemento.atom.name, ' Sz=', elemento.spin, ' Position:', elemento.position
      print '    Orbitals:'
      for orbital in elemento.atom.orbitals:
         print '    ',orbital

def data2file1(name,lista):
   archivo=open(name,'w')
   for fila in lista:
      for col in fila:
         archivo.write(str(col)+'  ')
      archivo.write('\n')
   archivo.close()
#   return name

def data2file(name,MAT):
   """
     MAT should have been constructed like this:
       MAT = np.matrix([X,Y,Z,...])
   """
   archivo=open(name,'w')
   MAT = MAT.transpose()
   for fila in MAT:
      archivo.write(str(fila).replace('[','').replace(']','')+'\n')
   archivo.close()
#   return name


def data2file2(name,lista):
   """
     This function returns the name in which the band structure is stored
    it will be a *.dat file
   """
   X,Y,Z = lista
   archivo=open(name,'w')
   for x,y,z in zip(X,Y,Z):
      archivo.write(str(x)+'   '+str(y)+'   '+str(z)+'\n')
#      archivo.write('\n')
   archivo.close()
#   return name


"""~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"""


def getname(orbital):
   if all(orbital[0][0] == np.array([0, 0, 1, 1])) or all(orbital[0][0] == np.array([0, 0, 1, -1])):
     name = 's'
   elif all(orbital[0][0] == np.array([ 1, -1,  1,  1])) or all(orbital[0][0] == np.array([ 1, -1,  1,  -1])):
     if orbital[0][1] == 1./np.sqrt(2.):
       name = 'px'
     else:
       name = 'py'
   elif all(orbital[0][0] == np.array([ 1, 0,  1,  1])) or all(orbital[0][0] == np.array([ 1, 0,  1,  -1])):
     name = 'pz'
   else:
     name='error'
   return name


def getindex(base,place):
   index = 0
   for elemento in base[0:place]:
      index = index + orbitalesnum.get(elemento.atom.name)
   return index


def getlayer(elemento):
   layer=listalayer.get(elemento)


def antidiagonal(dim,elem):
   lista = [[None for i in range(dim)] for j in range(dim)]
   for i,j in zip(range(dim),range(dim-1,-1,-1)):
      lista[i][j] = elem
   antidiag = scipy.sparse.bmat(lista)
   return antidiag



#def time(label,time_old,f):
#   """
#     Prints the time elapsed since time_old. and returns the time at which we
#   measured the time for future measurements
#     pz = 0
#   """
#   time_new = time.clock()
#   f.write(str(' %s:' % (label))+'  '+str(time_new-time_old)+'\n')
#   time_old=time_new
#   return time_old


def rotation(v,Q,u=np.array([0,0,1]),deg=True):
   """
     Rotates a vector v by an angle Q around the axis u
     Q in degrees by default
   """
   u = u/np.linalg.norm(u)
   ux = u[0]
   uy = u[1]
   uz = u[2]
   #print '==>',ux,uy,uz
   if deg :
     q = np.radians(Q)
   else:
     q = Q
   cos1 = 1.-np.cos(q)
   R = np.array(
   [[np.cos(q)+(ux**2)*cos1 , ux*uy*cos1-uz*np.sin(q), ux*uz*cos1+uy*np.sin(q)],
    [uy*ux*cos1+uz*np.sin(q), np.cos(q)+(uy**2)*cos1 , uy*uz*cos1-ux*np.sin(q)],
    [uz*ux*cos1-uy*np.sin(q), uz*uy*cos1+ux*np.sin(q), np.cos(q)+(uz**2)*cos1]])
   vec = np.dot(R, v)
   return vec


def find_closest(r,base):
   dists = []
   for e in base:
      dists.append((np.linalg.norm(e.position-r),e.place))
   return min(dists)[1]


def mapeo(x,a,b,x0=0,x1=1):
   """
     recta que mapea el intervalo [x0,x1] --> [a,b]
     by default [0,1] --> [a,b]
     returns X such that:
                 a<x<b ---->  x0<X<x1
   """
   m = (float(a)-float(b))/(float(x0)-float(x1))
   n = a-m*x0
   return m*x+n


def gap(H):
   """
     Returns the gap between the two middle eigenvalues of H
   """
   eig,eigvecs = np.linalg.eigh(H)
   gap = eig[len(eig)/2] - eig[len(eig)/2-1]
   return gap


def norm(v,a=1):
   return v*a/np.linalg.norm(v)


""" ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ OBSOLETOS ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ """
def read_params(archivo):
   f = open(sys.argv[1], "r")
   params=[]
   while True:
        linea = f.readline()
        if not linea: break
        linea = linea.lstrip()
        if linea[0] != '#':
          params.append(linea.rstrip('\n').split("=")[0].lstrip().rstrip())
          params.append(linea.rstrip('\n').split("=")[1].lstrip().rstrip())
   f.close()
   return params
