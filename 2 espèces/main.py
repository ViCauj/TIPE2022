import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from random import choices
from random import randint 
from random import shuffle


plt.set_cmap('hot')
#plt.set_cmap('inferno')

#Recup Map
f = open("Map\config.txt") 
L = int(f.readline())
l = int(f.readline())
r = int(f.readline())
pL = int(f.readline())
pl = int(f.readline())
f.close()

PetiteMap = np.loadtxt("Map/PetiteMap.txt")
GrandeMap = np.loadtxt("Map/GrandeMap.txt")
MapPassageFlatten = np.loadtxt("Map/MapPassage.txt")
f = open("Map/MapPassageShape.txt")
MapPassageShape = list(map(int, f.readline().replace("(", "").replace(")", "").split(", ")))
f.close()
MapPassage = np.reshape(MapPassageFlatten, MapPassageShape)

FichierMapDistance = open("Map/mapDistance1.txt", "r")
L = FichierMapDistance.readlines()
FichierMapDistance.close()

MapDistance1 = np.zeros((pL,pl))
for i in range(pL):
    L[i] = L[i].split()
    for j in range(pl):
        if L[i][j] == "inf":
            MapDistance1[i, j] = np.inf
        else:
            MapDistance1[i, j] = float(L[i][j])
            
FichierMapDistance = open("Map/mapDistance2.txt", "r")
L = FichierMapDistance.readlines()
FichierMapDistance.close()

MapDistance2 = np.zeros((pL,pl))
for i in range(pL):
    L[i] = L[i].split()
    for j in range(pl):
        if L[i][j] == "inf":
            MapDistance2[i, j] = np.inf
        else:
            MapDistance2[i, j] = float(L[i][j])
print("fin récup info")

#cste
beta = 10
alpha = 0.7 #proba dispariton dBoson
J0 = 1
Js = 2
Jd = 3
ColPers1 = 0.4
ColPers2 = 0.6
ColPersE = 0.1
ColPersH = 0.9
chgmtEtat = 1
#

# MvtPos = np.ones((3, 3))
# MvtPos[0, 0], MvtPos[0, 2], MvtPos[2, 0], MvtPos[2, 2] = np.inf, np.inf, np.inf, np.inf

MapDyn1 = np.zeros((pL, pl))
MapDyn2 = np.zeros((pL, pl))
LdBosons = []

def S(i, j, k, ETAT):
    MvtPos = np.ones((3, 3))
    if ETAT == "HEUREUX":
        MvtPos[0, 0], MvtPos[0, 2], MvtPos[2, 0], MvtPos[2, 2] = np.inf, np.inf, np.inf, np.inf
    Res = np.zeros((3,3))
    for a in range(-1, 2):
        for b in range(-1, 2):
            if (0 <= a+i < pL) and (0 <= b+j < pl):
                if k == 1:
                    Res[a + 1, b + 1] = MapDistance1[a + i, b + j]*MvtPos[a + 1, b + 1]
                else:
                    Res[a + 1, b + 1] = MapDistance2[a + i, b + j]*MvtPos[a + 1, b + 1]
            else:
                Res[a + 1, b + 1] = np.inf
    return Res

def D(i, j, k, ETAT):
    MvtPos = np.ones((3, 3))
    if ETAT == "HEUREUX":
        MvtPos[0, 0], MvtPos[0, 2], MvtPos[2, 0], MvtPos[2, 2] = np.inf, np.inf, np.inf, np.inf
    Res = np.zeros((3,3))
    for a in range(-1, 2):
        for b in range(-1, 2):
            if (0 <= a+i < pL) and (0 <= b+j < pl):
                if k == 1:
                    Res[a + 1, b + 1] = MapDyn1[a + i, b + j]*(MvtPos[a + 1, b + 1] != np.inf)
                else:
                    Res[a + 1, b + 1] = MapDyn2[a + i, b + j]*(MvtPos[a + 1, b + 1] != np.inf)
            else:
                Res[a + 1, b + 1] = 0
    return Res

def occup(i, j):
    Res = np.zeros((3,3))
    for a in range(-1, 2):
        for b in range(-1, 2):
            if (0 <= a+i < pL) and (0 <= b+j < pl):
                if not (PetiteMap[a + i, b + j] in [0, 1]):
                    Res[a + 1, b + 1] = 1
            else:
                Res[a + 1, b + 1] = 0
    Res[1, 1] = 0
    return Res

def P(pos, posAvant, k, ETAT):
    i, j = pos[0], pos[1]
    d = D(i, j, k, ETAT)
    s = S(i, j, k, ETAT)
    o = occup(i, j)
    p = np.zeros((3,3))
    inertie = np.ones((3,3))
    inertie[posAvant[0]-i+1, posAvant[1]-j+1] = np.exp(-beta*Jd)
    inertie[1 - posAvant[0] + i, 1 - posAvant[1] + j] = np.exp(beta*J0)
    for a in range(3):
        for b in range(3):
            if np.inf in [s[a, b], s[1, 1]]:
                TauS = np.inf
            else:
                TauS = s[a, b] - s[1, 1]
            TauD = d[a, b] - d[1, 1]
            p[a, b] = np.exp(-beta*Js*TauS)*np.exp(beta*Jd*TauD)*inertie[a, b]*(1-o[a, b])
    p = (1/np.sum(p))*p
    
    Res = []
    Weights = []
    for a in range(3):
        for b in range(3):
            Res.append([a, b])
            Weights.append(p[a, b])
    ind = Res.index(choices(Res, Weights)[0])
    return [Res[ind][0]-1+i, Res[ind][1]-1+j], Weights[ind]

def colGrandeMap(i, j, colExt, colInt):
    L = MapPassage[i, j]
    millieuI, millieuJ = (L[0][0] + L[-1][0])/2, (L[0][1] + L[-1][1])/2
    rayon =  13*((((L[0][0] - L[-1][0])**2)/2)**0.5)/20
    for pos in L:
        distCarre = (pos[0] - millieuI)**2 + (pos[1] - millieuJ)**2
        if distCarre <= rayon**2:
            if distCarre <= (3*rayon/4)**2:
                GrandeMap[int(pos[0]), int(pos[1])] = colInt
            else:
                GrandeMap[int(pos[0]), int(pos[1])] = colExt
    
def init(nombrePers):
    if np.count_nonzero(PetiteMap == 1) < nombrePers:
        print("pas assez de place")
        return -1
    res = []
    while len(res) < nombrePers:
        pos = [randint(1, pL-1), randint(1, pl-1)]
        while PetiteMap[pos[0], pos[1]] != 1:
            pos = [randint(1, pL-1), randint(1, pl-1)]
        k = randint(1,2)
        if k == 1:
            PetiteMap[pos[0], pos[1]] = ColPers1
            colGrandeMap(pos[0], pos[1], ColPersH, ColPers1)
            res.append([MapDistance1[pos[0], pos[1]], [pos], ["HEUREUX", 0], k]) 
        else:
            PetiteMap[pos[0], pos[1]] = ColPers2
            colGrandeMap(pos[0], pos[1], ColPersH, ColPers2)
            res.append([MapDistance2[pos[0], pos[1]], [pos], ["HEUREUX", 0], k])
    return res

#Simu
nbrSor = 0

DebutSortie = pL//2 - 1
FinSortie = pL//2 + 1

nombrePers = 1
Lpos = init(150)

fig = plt.figure(dpi = 300)
ax1=fig.add_subplot(3,1,1)
ax2=fig.add_subplot(3,1,2)
ax3=fig.add_subplot(3,1,3)

ax1.set_axis_off()
ax2.set_axis_off()
ax3.set_axis_off()

ims = []

CompteurImmobile = 0
while len(Lpos) > 0 and CompteurImmobile < 3:
    immobile = 0
    LposSuiv = []
    probaPosSuiv = []
    for pos in Lpos:
        if len(pos[1]) == 1:
            posSuiv = P(pos[1][-1], pos[1][-1], pos[-1], pos[-2][0])
        else:
            posSuiv = P(pos[1][-1], pos[1][-2], pos[-1], pos[-2][0])
            
        # etat = []
        # premierPass = True
        # if pos[-1] == 1:
        #     while (posSuiv[0] in LposSuiv) and (posSuiv[0] != pos[1][-1]):
        #         if premierPass:
        #             MapDistance1[pos[1][-1][0], pos[1][-1][1]] += 1
        #             premierPass = False
        #         etat.append([PetiteMap[posSuiv[0][0], posSuiv[0][1]], posSuiv[0]])
        #         PetiteMap[posSuiv[0][0], posSuiv[0][1]] = ColPers1
        #         if len(pos[1]) == 1:
        #             posSuiv = P(pos[1][-1], pos[1][-1], pos[-1], pos[-2][0])
        #         else:
        #             posSuiv = P(pos[1][-1], pos[1][-2], pos[-1], pos[-2][0])
        #     for i in etat:
        #         PetiteMap[i[1][0], i[1][1]] = i[0]
        #     if not premierPass:
        #         MapDistance1[pos[1][-1][0], pos[1][-1][1]] -= 1
        # else:
        #     while (posSuiv[0] in LposSuiv) and (posSuiv[0] != pos[1][-1]):
        #         if premierPass:
        #             MapDistance2[pos[1][-1][0], pos[1][-1][1]] += 1
        #             premierPass = False
        #         etat.append([PetiteMap[posSuiv[0][0], posSuiv[0][1]], posSuiv[0]])
        #         PetiteMap[posSuiv[0][0], posSuiv[0][1]] = ColPers2
        #         if len(pos[1]) == 1:
        #             posSuiv = P(pos[1][-1], pos[1][-1], pos[-1], pos[-2][0])
        #         else:
        #             posSuiv = P(pos[1][-1], pos[1][-2], pos[-1], pos[-2][0])
        #     for i in etat:
        #         PetiteMap[i[1][0], i[1][1]] = i[0]
        #     if not premierPass:
        #         MapDistance2[pos[1][-1][0], pos[1][-1][1]] -= 1
        
        if posSuiv[0] in LposSuiv:
            posSuiv[0][0], posSuiv[0][1] = pos[1][-1][0], pos[1][-1][1]
        LposSuiv.append(posSuiv[0])
        if pos[-2][1] >= chgmtEtat:
            pos[-2][0] = "ENERVE"
        elif pos[-2][1] <= -chgmtEtat:
            pos[-2][0] = "HEUREUX"
        if pos[-2][1] == chgmtEtat + 1:  
            if pos[1][-1] != LposSuiv[-1]:
                pos[-2][1] -= 1
        elif pos[-2][1] == -(chgmtEtat + 1):
            if pos[1][-1] == LposSuiv[-1]:
                pos[-2][1] += 1
        else:
            if pos[1][-1] == LposSuiv[-1]:
                pos[-2][1] += 1
            else:
                pos[-2][1] -= 1
        probaPosSuiv.append(posSuiv[1])
        
        if len(pos[1]) > 1 and pos[1][-1] == pos[1][-2]:
            immobile += 1
            
    if immobile == len(Lpos):
        CompteurImmobile += 1
        print("Bloqué")
    else:
        CompteurImmobile = 0
        
    for pos in Lpos:
        PetiteMap[pos[1][-1][0], pos[1][-1][1]] = 1
        colGrandeMap(pos[1][-1][0], pos[1][-1][1], 1, 1)
        
    BosonsARetirer = []
    for PosdBoson in LdBosons:
        meurt = choices([True,False], [alpha, 1-alpha])[0]

        if meurt:
            if PosdBoson[-1] == 1:
                MapDyn1[PosdBoson[0][0], PosdBoson[0][1]] -= 1
                if MapDyn1[PosdBoson[0][0], PosdBoson[0][1]] == 0:
                    BosonsARetirer.append(PosdBoson)
            else:
                MapDyn2[PosdBoson[0][0], PosdBoson[0][1]] -= 1
                if MapDyn2[PosdBoson[0][0], PosdBoson[0][1]] == 0:
                    BosonsARetirer.append(PosdBoson)
    for PosdBoson in BosonsARetirer:
        LdBosons.remove(PosdBoson)
    for i in range(len(LposSuiv)):
        #actualise MapDyn
        if LposSuiv[i] != Lpos[i][1][-1]: #que si il bouge
            if Lpos[i][-1] == 1:
                if MapDyn1[Lpos[i][1][-1][0], Lpos[i][1][-1][1]] == 0:
                    LdBosons.append([[Lpos[i][1][-1][0], Lpos[i][1][-1][1]], 1])
                    MapDyn1[Lpos[i][1][-1][0], Lpos[i][1][-1][1]] += 1
                else:
                    MapDyn1[Lpos[i][1][-1][0], Lpos[i][1][-1][1]] += 1
            else:
                if MapDyn2[Lpos[i][1][-1][0], Lpos[i][1][-1][1]] == 0:
                    LdBosons.append([[Lpos[i][1][-1][0], Lpos[i][1][-1][1]], 2])
                    MapDyn2[Lpos[i][1][-1][0], Lpos[i][1][-1][1]] += 1
                else:
                    MapDyn2[Lpos[i][1][-1][0], Lpos[i][1][-1][1]] += 1
            
        #Actualise Map
        Lpos[i][1].append(LposSuiv[i])
        if Lpos[i][-1] == 1:
            Lpos[i][0] = MapDistance1[LposSuiv[i][0], LposSuiv[i][1]]
        else :
            Lpos[i][0] = MapDistance2[LposSuiv[i][0], LposSuiv[i][1]]
    for pos in Lpos:
        if pos[0] == 1:
            Lpos.remove(pos)
            nbrSor += 1
            print(len(Lpos))
        else:
            if pos[-1] == 1:
                if pos[-2][0] == "HEUREUX":
                    PetiteMap[pos[1][-1][0], pos[1][-1][1]] = ColPers1
                    colGrandeMap(pos[1][-1][0], pos[1][-1][1], ColPersH, ColPers1)
                else:
                    PetiteMap[pos[1][-1][0], pos[1][-1][1]] = ColPers1
                    colGrandeMap(pos[1][-1][0], pos[1][-1][1], ColPersE, ColPers1)
            else:
                if pos[-2][0] == "HEUREUX":
                    PetiteMap[pos[1][-1][0], pos[1][-1][1]] = ColPers2
                    colGrandeMap(pos[1][-1][0], pos[1][-1][1], ColPersH, ColPers2)
                else:
                    PetiteMap[pos[1][-1][0], pos[1][-1][1]] = ColPers2
                    colGrandeMap(pos[1][-1][0], pos[1][-1][1], ColPersE, ColPers2)
    #Lpos.sort() #pas obligatoire 
    shuffle(Lpos)
    
    # im = plt.imshow(GrandeMap, animated=True)
    # ims.append([im])
    
    im1 = ax1.imshow(GrandeMap, animated = True)
    im2 = ax2.imshow(MapDyn1, animated = True)
    im3 = ax3.imshow(MapDyn2, animated = True)
    ims.append([im1, im2, im3])

print("Fin Simulation")
ani = animation.ArtistAnimation(fig, ims, interval=1)
f = r"animation.gif" 
writergif = animation.PillowWriter(fps=15)
ani.save(f, writer=writergif)
plt.show()
            
        
### TODO :
###     .rajouter bosons dynamique négatif pour les gens bloqué
