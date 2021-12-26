import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from matplotlib.cm import get_cmap

from random import randint 
from random import shuffle


plt.set_cmap('inferno')

####################################################################################################
##########################################Recuperation maps#########################################
#################################################################################################### 


FichierMapDistance = open("Map/mapDistance.txt", "r")
L = FichierMapDistance.readlines()
FichierMapDistance.close()

n = len(L) 
r = 1
R = 3
iDebutSortie = n//2 - n//20
iFinSorite = n//2 + n//20
Sortie = [[i, 0] for i in range(iDebutSortie, iFinSorite + 1)]

MapDistance = np.zeros((n,n))
for i in range(n):
    L[i] = L[i].split()
    for j in range(n):
        if L[i][j] == "inf":
            MapDistance[i, j] = np.inf
        else:
            MapDistance[i, j] = float(L[i][j])
            
FichierMap = open("Map/Map.txt", "r")
L = FichierMap.readlines()
FichierMap.close()

Map = np.zeros((n,n))
for i in range(n):
    L[i] = L[i].split()
    for j in range(n):
        if L[i][j] == "inf":
            Map[i, j] = np.inf
        else:
            Map[i, j] = float(L[i][j])
        
####################################################################################################
####################################################################################################
####################################################################################################   

def affiche(X):
    plt.axis("off")
    plt.imshow(X)
    plt.show()
    
def cond(x):
    for i in x:
        if i<0 or i>=n:
            return False
    if MapDistance[x[0], x[1]] == np.inf:
        return False
    return True

def bougeSansCol(pos, r, R, posOccup):
    x, y = pos[0], pos[1]
    # MovPoss = []                                  #deplacement diag
    # for k in [x + i for i in range(-r, r + 1)]:
    #     for l in [y + i for i in range(-r, r + 1)]:
    #         MovPoss.append([k, l])
    MovPoss = []                                    #deplacement ds un cercle de rayon R
    for k in [x + i for i in range(-R, R + 1)]:
        for l in [y + i for i in range(-R, R + 1)]:
            if (k - x)**2 + (l - y)**2 <= R**2:
                MovPoss.append([k, l])
    #MovPoss = [[x+1, y], [x, y+1], [x-1, y], [x, y-1]] #pas de mov diag
    PosPoss = [i for i in MovPoss if cond(i)]
    r = [np.inf, 0]
    for i in PosPoss:
        if MapDistance[i[0], i[1]] < r[0]:
            r = [MapDistance[i[0], i[1]], i]
    if r[-1] == 0:
        print("Mauvaise CI")
        return -1
    return r[-1]


def bougeAvecCol(pos, r, R, posOccup):
    x, y = pos[0], pos[1]
    colPos = []
    for POS in posOccup:
        D = ((x - POS[0])**2 + (y - POS[1])**2)**0.5
        if D < 2*r + R:
            colPos.append(POS)
        
    # MovPoss = []                                  #deplacement diag
    # for k in [x + i for i in range(-r, r + 1)]:
    #     for l in [y + i for i in range(-r, r + 1)]:
    #         MovPoss.append([k, l])
    MovPoss = []                                    #deplacement ds un cercle de rayon R
    for k in [x + i for i in range(-R, R + 1)]:
        for l in [y + i for i in range(-R, R + 1)]:
            if (k - x)**2 + (l - y)**2 <= R**2:
                MovPoss.append([k, l])
    #MovPoss = [[x+1, y], [x, y+1], [x-1, y], [x, y-1]] #pas de mov diag
    PosPoss = [i for i in MovPoss if cond(i)]
    
    PossPossFinal = []
    
    for i in PosPoss:
        k = True
        for j in colPos:
            D = ((i[0] - j[0])**2 + (i[1] - j[1])**2)**0.5
            if D < 2*r:
                k = False
        if k:
            PossPossFinal.append(i)
        
    r = [MapDistance[pos[0], pos[1]] + 1, pos] #on rajoute 1 car on veut que les gens bouge un max
    for i in PossPossFinal:
        if MapDistance[i[0], i[1]] < r[0]:
            r = [MapDistance[i[0], i[1]], i]
    if r[-1] == 0:
        print("Mauvaise CI")
        return -1
    return r[-1]

def colorRond(pos, col1, col2, r):
    x, y = pos[0], pos[1]
    MovPoss = []
    for k in [x + i for i in range(-r, r + 1)]:
        for l in [y + i for i in range(-r, r + 1)]:
            if (k - x)**2 + (l - y)**2 <= r**2:
                MovPoss.append([k, l])
    PosPoss = [i for i in MovPoss if cond(i)]
    for i in PosPoss:
        Map[i[0], i[1]] = col1
    for i in PosPoss:
        if (i[0] - x)**2 + (i[1] - y)**2 <= (3*r/4)**2:
            Map[i[0], i[1]] = col2
        
def colorCarre(pos, col1, col2, r):
    x, y = pos[0], pos[1]
    MovPoss = []
    for k in [x + i for i in range(-r, r + 1)]:
        for l in [y + i for i in range(-r, r + 1)]:
            MovPoss.append([k, l])
    PosPoss = [i for i in MovPoss if cond(i)]
    shuffle(PosPoss)
    for i in PosPoss:
        Map[i[0], i[1]] = col1
    Map[x, y] = col2
        
def actualise(pos, Lpos, posOccup):
    i = posOccup.index(pos[-1])
    if pos[-1] in Sortie:
        colorRond(pos[-1], 10*n, 10*n, r)
        Lpos.remove(pos)
        posOccup.remove(pos[-1])
    else:
        res = bougeAvecCol(pos[-1], r, R, posOccup)
        Lpos[i].append(res)
        posOccup[i] = res
        colorRond(pos[-1], 2*n, 6*n, r)
    
def init():
    pos = [[randint(1, n-1), randint(1, n-1)]]
    while Map[pos[-1][0], pos[-1][1]] != 10*n:
        pos[-1] = [randint(1, n-1), randint(1, n-1)]
    colorRond(pos[-1], 2*n, 6*n, r)
    return pos

def simu(n):
    Lpos = []
    for i in range(n):
        Lpos.append(init())
    return Lpos

fig = plt.figure()
            
Lpos = simu(1500)
posOccup = [i[-1] for i in Lpos]

ims = []
im = plt.imshow(Map, animated=True)
ims.append([im])

k = 0
while len(Lpos) > 0:
    for i in Lpos:
        colorRond(i[-1], 10*n, 10*n, r)
    for i in Lpos:
        actualise(i, Lpos, posOccup)
    
    im = plt.imshow(Map, animated=True)
    ims.append([im])
    k+=1

ani = animation.ArtistAnimation(fig, ims, interval=1, blit=True)


f = r"Animation/animation.gif" 
writergif = animation.PillowWriter(fps=30) 
ani.save(f, writer=writergif)
plt.show()
