import numpy as np
import matplotlib.pyplot as plt

plt.set_cmap('hot')

#rem mieux si N impaire et r impaire

r = 31 
l = 100*r
L = 6*r
    
pL = L//r
pl = l//r

f = open("Map/config.txt", "w")
f.write(str(L) + "\n" + str(l) + "\n" + str(r) + "\n" + str(pL) + "\n" + str(pl))
f.close()


MapPassage = []
for i in range(pL):
    RES = []
    for j in range(pl):
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

PetiteMap = np.ones((pL, pl), int)
GrandeMap = np.ones((L, l), int)

plt.imshow(GrandeMap)

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
ObstacleRectangle(GrandeMap, 0, 0, L, r, 0)
ObstacleRectangle(GrandeMap, 0, 0, r, l, 0)
ObstacleRectangle(GrandeMap, 0, l-r, L, l, 0)
ObstacleRectangle(GrandeMap, L-r, 0, L, l, 0)

#bordures PetiteMap
ObstacleRectangle(PetiteMap, 0, 0, pL, 1, 0)
ObstacleRectangle(PetiteMap, 0, 0, 1, pl, 0)
ObstacleRectangle(PetiteMap, 0, pl-1, pL, pl, 0)
ObstacleRectangle(PetiteMap, pL-1, 0, pL, pl, 0)
#


DebutGSortie = L//2 - r
FinGSortie = L//2 + r

DebutPSortie = pL//2 - 1
FinPSortie = pL//2 + 1

#Obstacles
#Obstacles GrandeMap
ObstacleRectangle(GrandeMap, DebutGSortie, 0, FinGSortie, r, 1) #sortie
ObstacleRectangle(GrandeMap, DebutGSortie, l-r, FinGSortie, l, 1)
#ObstacleRond(GrandeMap, N//2, r + N//5, N//10, 0)

#Obstacles PetiteMap
ObstacleRectangle(PetiteMap, DebutPSortie, 0, FinPSortie, 1, 1) #sortie
ObstacleRectangle(PetiteMap, DebutPSortie, pl-1, FinPSortie, pl, 1)
#ObstacleRond(PetiteMap, n//2, 1 + n//5, n//10, 0)

#

plt.axis("off")
plt.imshow(GrandeMap)
plt.savefig("Map/Map.png",dpi=1000) 
np.savetxt("Map/GrandeMap.txt", GrandeMap)
np.savetxt("Map/PetiteMap.txt", PetiteMap)

plt.imshow(PetiteMap)
plt.show()
