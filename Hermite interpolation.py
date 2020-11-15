# %%
import matplotlib.pyplot as plt
import numpy as np
import math
from scipy.misc import derivative



def roznica_dzielona(pseudoargumenty=[],wartosci_dane=[], pochodne=[],poziom=0,wsp=[]):
    """Funkcja zwraca współczynniki wielomianu w bazie Newtona (wartości na przekątnej macierzy)
    Dostaje juz opowiednio zmodyfikowane argumenty, wartosci i liste pochodnych"""

    wartosci=wartosci_dane
    if(poziom==0):
        wsp1=[]
    else:
        wsp1=wsp
    wsp1.append(wartosci[0])
    deriv=pochodne
    wartosci_next=[]
    
    for i in range(0, len(wartosci)-1):
        if pseudoargumenty[i+1+poziom]!=pseudoargumenty[i]:
            a=(wartosci[i+1]-wartosci[i])/(pseudoargumenty[i+1+poziom]-pseudoargumenty[i])
            wartosci_next.append(a)
        else:
            wartosci_next.append(deriv[0])
            deriv.pop(0)
   

    if(len(pseudoargumenty)==len(wsp1)):
        return wsp1
    roznica_dzielona(pseudoargumenty, wartosci_next,deriv, poziom+1, wsp1)
    
    
    return wsp1



def Hermit(wezly):
    """Pierwsza czesc funkcji przygotowuje listy: pseudoargumentow, wartosci i pochodnych.
    Funkcja Hermit dostaje liste list. W kazdej wewnetrznej liscie sa kolejno: argument, wartosc,
    wartosci kolejnych pochodnych...
    """
    arg=[]
    wart=[]
    poch=[]
    for i in range(0, len(wezly)):
        for j in range (0, len(wezly[i])-1):
            arg.append(wezly[i][0])
            wart.append(wezly[i][1])
                 

    maxi=len(wezly[0])
    for i in range(1, len(wezly)):
        if maxi<len(wezly[i]):
            maxi=len(wezly[i])
#n wskauje ktora to pochodna
#h ile razy n ta pochodna ma byc wrzucona
    for n in range(2,maxi):
        for i in range(0, len(wezly)):
            if n<len(wezly[i]):
                for h in range (n, len(wezly[i])):
                    poch.append(wezly[i][n]/math.factorial(n-1))


    return roznica_dzielona(arg, wart, poch)


def argi(wezly):
    """Funkcja z listy argumentow tworzy liste pseudoargumentow"""
    argsy=[]
    for i in range(0,len(wezly)):
        
            for j in range (0,len(wezly[i])-1):
                argsy.append(wezly[i][0])
    return argsy                        


def wartosc_wielomianu(argumenty, wspolczynniki,x):
    """Funkcja zwraca wartosc wielomianiu interpolacyjnego Herimte'a dla zadanego punktu. Argumentami sa 
    wezly i wspolczynniki w bazie Newtona oraz x dla ktorego chcemy policzyc wartosc"""
    wart=wspolczynniki[-1]
    for i in range(len(argumenty)-1,0,-1):
       
        wart=wart*(x-argumenty[i-1])+wspolczynniki[i-1]
        #print(wart)
        
    return wart       

# %%
"""
Przyklad z angielskiej wiki
"""

# %%

wezly=[[-1, 2,-8,56],[0,1,0,0],[1,2,8,56]]

def f(x):
    return ((x**8+1))
wsp=Hermit(wezly)

xvals = np.linspace(-3,3,100)
yvals=[f(x) for x in xvals]

plt.plot(xvals,[wartosc_wielomianu(argi(wezly),wsp,x) for x in xvals])
plt.plot(xvals,yvals,c="green")
plt.scatter([wezly[i][0] for i in range(0, len(wezly))], [wezly[i][1] for i in range(0, len(wezly))], c="red")
plt.ylim(-4,4)
plt.show()

# %%
"""
rys dla podanych wezlow

"""

# %%
wezly=[[0,1], [1,2], [2,5,4]]



xvals = np.linspace(0,5,100)
yvals=[x**2+1 for x in xvals]

plt.plot(xvals, yvals, c="green")
wsp=Hermit(wezly)

plt.plot(xvals,[wartosc_wielomianu(argi(wezly),wsp,x) for x in xvals])

plt.scatter([wezly[i][0] for i in range(0, len(wezly))], [wezly[i][1] for i in range(0, len(wezly))], c="red")
                                                        
plt.show()

# %%
"""
Przyklad
"""

# %%

wezly= [[-1, 0, 0], [0, 1, 1], [1, 0]]

wsp=Hermit(wezly)

xvals = np.linspace(-5,5,100)


plt.plot(xvals,[wartosc_wielomianu(argi(wezly),wsp,x) for x in xvals])
plt.scatter([wezly[i][0] for i in range(0, len(wezly))], [wezly[i][1] for i in range(0, len(wezly))], c="red")
plt.ylim(-4,4)

plt.show()

# %%
"""
Przyklad z liczeniem pochodnej
"""

# %%
argumenty=[-2,-1,0,1,2]
krotnosci=[2,1,2,1,2]

def f(x):
    return ((x**4+x**2-1))
    
wezly=[]
for i in range(0,len(argumenty)):
    wezly.append([])
    wezly[i].append(argumenty[i])
    wezly[i].append(f(argumenty[i]))
                     
    for j in range (1, krotnosci[i]):
        wezly[i].append(derivative(f,argumenty[i],n=j))    


wsp=Hermit(wezly)

xvals = np.linspace(-5,5,100)
yvals=[f(x) for x in xvals]

plt.plot(xvals, yvals, c="green")
plt.plot(xvals,[wartosc_wielomianu(argi(wezly),wsp,x) for x in xvals])
plt.scatter([wezly[i][0] for i in range(0, len(wezly))], [wezly[i][1] for i in range(0, len(wezly))], c="red")
plt.ylim(-1,10)

plt.show()

