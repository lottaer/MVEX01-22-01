import random as rnd
import numpy as np
import matplotlib.pyplot as plt
import time

#22 jan 2022
#Enkel Galton-Watson process
#kandidat1.6.1
#MÅL:Graf över ålder och iterationer

def simGen(hour,p0,p1,p2):
    z = []
    z.append(1)
    #Simulera population över 10 timmar
    for t in range(hour):
        #print(z)
        z.append(0)
        #För varje cell i timme[i]
        #Kolla celldelning
        #här p0 = 0.3, p1 = 0.2, p2 = 0.5
        for j in range(z[t]):
            p = rnd.random()
            if p > p0+p1:
                z[t+1]+=2
            elif p > p0 and p < p0+p1:
                z[t+1]+=1
    return z



def simGenProt(dstart,p,n,k,q,hour):
    d = []
    d.append([dstart, 0])
    for t in range(hour):
        dtid = [] #Hjälpvariabel för varje ny generation
        # cell = object i d[t][i][0]
        for i in range(len(d)):
            if rnd.random() < q: #Kolla om cellen delar samt att den inte har mer än k defekta protein
                #Fördela defekta protein
                d1 = d[i][0]+n
                d2 = 0
                for protein in range(d[i][0]+n):
                    if rnd.random() < p:
                        d1-=1
                        d2+=1
                if d1 <= k:
                    dtid.append([d1,d[i][1]]) # vi vill att denna ska ärva gen från innan
                if d2 <= k:
                    dtid.append([d2,t]) # cell föds vid t

            elif d[i][0]+2 <= k:
                dtid.append([d[i][0]+n,d[i][1]])
        d=dtid

    return d



def calcAVG(totvec,tot, it, hour):
    avgvec = []
    for i in range(len(totvec[0])):
        sum = 0
        for j in range(it):
            sum += totvec[j][i]
        sum = sum/it
        avgvec.append(sum)

    print(avgvec)
    x1 = np.array(avgvec)
    y1= np.linspace(1,hour,num=hour-1)

    #RÄknar ut genomsnitt av antal celler i varje åldersgrupp
    avgquote = []
    for age in avgvec:
        avgquote.append(age/(tot/it))


    # Tar bort de under 1% av åldersfördelningen innan det grafas
    newavg = []
    for i in avgquote:
        if i > 0.01:
          newavg.append(i)
    x2 = np.array(newavg)
    y2= np.linspace(1,len(newavg),num=len(newavg))
    return x2,y2



def main():
    # n = 0 och q = 1.0 ger 1024 celler stämmer alltså!
    iteration = 10000 #Antal simulringar
    totvec = []
    tot = 0
    for it in range(iteration):
        hour = 15 #dstart är vilken typ du start med, hour timmar du kör
                    #(dstart,p,n,k,q,hours):
        d= simGenProt(0,0.1,2,4,0.7,hour)
        agevec = []
        for i in range(hour):
            agevec.append(0)
        for cell in d:

             for i in range(hour):
                 if hour-cell[1] == i:
                     agevec[i]+=1
        tot+=len(d)
        del agevec[0]
        totvec.append(agevec)


    print(tot/iteration)
    #Räknar ut genomsnitt av varje ålder
    x2,y2 = calcAVG(totvec,tot,iteration,hour)





    #Gör grafer
    #fig, axs = plt.subplots(1)
    #axs[0].plot(x1,y1,'ro')
    print(x2)
    print(np.sum(x2))
    fig = plt.figure()
    #ax = fig.add_axes([0,0,1,len(x2)])
    plt.ylabel("Ålder av celler")
    plt.xlabel("Kvot av cell med ålder i mot hela populationen")
    plt.title("Åldersfördelning vid flertyp-GW")

    plt.plot(x2,y2,'bo')


    plt.show()
    #print(d)



main()
