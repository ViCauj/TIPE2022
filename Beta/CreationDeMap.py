import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from matplotlib.cm import get_cmap

#fig = plt.figure()

plt.set_cmap('inferno')

def affiche(X):
    plt.axis("off")
    plt.imshow(X)
    plt.show()

n = 100
iDebutSortie = n//2 - n//20
iFinSorite = n//2 + n//20

####################################################################################################
#########################################Création Image#############################################
####################################################################################################

def ObstacleTriangle(X, i0, j0, l, col):
    for i in range(l):
        for j in range(i, l-i):
            X[j + i0, i + j0] = col

def ObstacleRectangle(X, i0, j0, i1, j1, col):
    n = len(X)
    if i1 > n or j1 > n:
        print("OUT OF RANGE")
        return -1
    for I in range(i0, i1):
        for J in range(j0, j1):
            X[I,J] = col
            
def ObstacleRond(X, x, y, r, col):
    n = len(X)
    for k in [x + i for i in range(-r, r + 1)]:
        for l in [y + i for i in range(-r, r + 1)]:
            if (k - x)**2 + (l - y)**2 <= r**2:
                if (0 <= k < n) and ( 0<= l < n) :
                    X[k,l] = col

Map = np.zeros((n,n))
for i in range(len(Map)):
    for j in range(len(Map[i])):
        if i in [0, n-1] or j in [0, n-1]:
            Map[i, j] = 0
        else : 
            Map[i, j] = 10*n
            
ObstacleRectangle(Map, iDebutSortie, 0, iFinSorite + 1, 1, 10*n) #sortie

#ObstacleRectangle(Map, n//2, n//2, n//2, n//2, 10*n) #rush au centre

ObstacleRond(Map, n//2, n//2, 45, 0)
ObstacleRond(Map, n//2, n//2, 44, 10*n)
ObstacleRectangle(Map, n//2 -2, n//2 + 44, n//2+2, n//2+46, 10*n)

ObstacleRond(Map, n//2, n//2, 35, 0)
ObstacleRond(Map, n//2, n//2, 34, 10*n)
ObstacleRectangle(Map, n//2 -2, n//2 - 36, n//2+2, n//2-33, 10*n)

ObstacleRond(Map, n//2, n//2, 25, 0)
ObstacleRond(Map, n//2, n//2, 24, 10*n)
ObstacleRectangle(Map, n//2 -26, n//2 -2, n//2-23, n//2+2, 10*n)

ObstacleRond(Map, n//2, n//2, 15, 0)
ObstacleRond(Map, n//2, n//2, 14, 10*n)
ObstacleRectangle(Map, n//2 +13, n//2 -2, n//2+16, n//2+2, 10*n)

#ObstacleTriangle(Map, n//2-n//5, n//10, 2*n//5, 0)
#ObstacleRond(Map, n//2, 15, 8, 0)
#ObstacleRectangle(Map, 10, 10, 80, 30, 0)
#ObstacleRectangle(Map, 50, 45, 90, 85, 0)

plt.axis("off")
plt.imshow(Map)
plt.savefig("Map/Map.png",dpi=1000)
plt.show()
np.savetxt('Map/Map.txt', Map) 

####################################################################################################
####################################################################################################
####################################################################################################

####################################################################################################
########################################Carte des distances#########################################
####################################################################################################
            
MapDistance = np.zeros((n, n)) #carte des distances
for i in range(n):
    for j in range(n):
        if Map[i, j] == 0: #si c'est un obstacle
            MapDistance[i, j] = np.inf #innatteignable
            
for i in range(iDebutSortie, iFinSorite + 1):
    MapDistance[i, 0] = 1
    
def voisin(x):
    i, j = x[0], x[1]
    # A = [i, i-1, i + 1]
    # B = [j, j-1, j + 1]
    # X = []
    # for k in A:
    #     for l in B:
    #         X.append([k, l])
    X = [[i-1, j], [i+1, j], [i, j+1], [i, j-1]] #pas de deplacement diag
    res = []
    for x in X:
        t = False
        for k in x:
            if k>=0 and k<n:
                t = True
        if t and (MapDistance[x[0], x[1]] == 0):
            res.append(x)
    return res   

#ims = []
#im = plt.imshow(MapDistance, animated=True)
#ims.append([im])
    
def dist(L, k):
    if len(L) == 0:
        return 0
    else :
        R = []
        k += 1
        for x in L:
            MapDistance[x[0], x[1]] = k
            vois = voisin(x)
            for y in vois:
                if not y in R:
                    R.append(y)
        #im = plt.imshow(MapDistance, animated=True)
        #ims.append([im])
        return dist(R, k)
    
#ani = animation.ArtistAnimation(fig, ims, interval=0, blit=True)
#plt.show()

dist([[i, 0] for i in range(iDebutSortie, iFinSorite + 1)], 0)
#dist([[n//2,n//2]], 0) #rush au centre

####################################################################################################
####################################################################################################
####################################################################################################


for i in range(n):
    for j in range(n):
        if MapDistance[i,j] == np.inf:
            Map[i, j] = 0
        else:
            Map[i, j] = 6*n - MapDistance[i, j]
            
np.savetxt('Map/mapDistance.txt', MapDistance)   

plt.axis("off")
plt.imshow(Map)
plt.savefig("Map/MapGradient.png",dpi=1000) 
plt.show()
