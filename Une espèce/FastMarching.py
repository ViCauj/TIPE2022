import numpy as np
import matplotlib.pyplot as plt

plt.set_cmap('hot')
#rem mieux si N impaire et r impaire

r = 31 #int(input("rayon d'une personne (doit divisser " + str(N) +") : "))
N = 33*r
while N%r != 0:
    r = int(input("rayon d'une personne (doit divisser " + str(N) +") : "))
    
n = N//r

f = open("Map/config.txt", "w")
f.write(str(N) + "\n" + str(r) + "\n" + str(n))
f.close()


MapPassage = []
for i in range(n):
    RES = []
    for j in range(n):
        res = []
        for a in range(r*i, r*(i + 1)):
            for b in range(r*j, r*(j + 1)):
                res.append([a, b])
        RES.append(res)
    MapPassage.append(RES)
MapPassage = np.array(MapPassage)

np.savetxt("Map\MapPassage.txt", MapPassage.flatten())
f = open("Map\MapPassageShape.txt", "w")
f.write(str(MapPassage.shape))
f.close()

print("fin setup coord correspondantes")

PetiteMap = np.ones((n, n), int)
GrandeMap = np.ones((N, N), int)

#
def ObstacleRectangle(X, i0, j0, i1, j1, couleur):
    for I in range(i0, i1):
        for J in range(j0, j1):
            X[I,J] = couleur
        
def ObstacleRond(X, x, y, r, couleur):
    for k in [x + i for i in range(-r, r + 1)]:
        for l in [y + i for i in range(-r, r + 1)]:
            if (k - x)**2 + (l - y)**2 <= r**2:
                X[k,l] = couleur
#
#bordures
#bordures GrandeMap
ObstacleRectangle(GrandeMap, 0, 0, N, r, 0)
ObstacleRectangle(GrandeMap, 0, 0, r, N, 0)
ObstacleRectangle(GrandeMap, 0, N-r, N, N, 0)
ObstacleRectangle(GrandeMap, N-r, 0, N, N, 0)

#bordures PetiteMap
ObstacleRectangle(PetiteMap, 0, 0, n, 1, 0)
ObstacleRectangle(PetiteMap, 0, 0, 1, n, 0)
ObstacleRectangle(PetiteMap, 0, n-1, n, n, 0)
ObstacleRectangle(PetiteMap, n-1, 0, n, n, 0)
#

DebutGSortie = N//2 - N//20
FinGSortie = N//2 + N//20

DebutPSortie = n//2 - n//20
FinPSortie = n//2 + n//20

iCercle, jCercle, rCercle = n//2, 1 + n//5, n//10

#Obstacles
#Obstacles GrandeMap
ObstacleRectangle(GrandeMap, DebutGSortie, 0, FinGSortie + 1, r, 1) #sortie
ObstacleRond(GrandeMap, r*iCercle + int(0.5*r), r*jCercle + int(0.5*r), (r-1)*rCercle, 0) #-1 pour ne pas dépasser du vrai cercle

#Obstacles PetiteMap
ObstacleRectangle(PetiteMap, DebutPSortie, 0, FinPSortie + 1, 1, 1) #sortie
ObstacleRond(PetiteMap, n//2, 1 + n//5, n//10, 0)

#

plt.axis("off")
plt.imshow(GrandeMap)
plt.savefig("Map/Map.png",dpi=1000) 

np.savetxt("Map/GrandeMap.txt", GrandeMap)
np.savetxt("Map/PetiteMap.txt", PetiteMap)

# plt.axis("off")
# plt.imshow(GrandeMap)
# plt.savefig("Map/Map.png",dpi=1000) 

# np.savetxt("Map/GrandeMap.txt", GrandeMap)
# np.savetxt("Map/PetiteMap.txt", PetiteMap)
