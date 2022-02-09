import random as rnd
#22 jan 2022
#Enkel Galton-Watson process
#kandidat1.5 mål, d-protein egenskaper
#Istället för att bestämma sannolikheter för t.ex död ska det simuleras!
#p och q istället, p odds för att defekt protein överförs till dotter
#q odds för att moder med mindre än k protein förökas.


##############################
#Den här används för flertyp Galton-Watson
#Simgen är functionen för enkel GW
#simGenProt simulerar flertyp
#############################




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



def simGenProt(dstart,p,n,k,q):
    d = []
    d.append([dstart])
    #Simulera population över 10 timmar
    for t in range(10):

        dtid = [] #Hjälpvariabel för varje ny generation
        for cell in d[t]:
            pr = rnd.random()
            if pr < q and cell < k: #Kolla om cellen delar samt att den inte
                                        #har mer än k defekta protein
                #Fördela defekta protein
                d1 = cell+n
                d2 = 0
                for protein in range(cell):
                    pd = rnd.random()
                    if pd < p:
                        d1-=1
                        d2+=1
                dtid.append(d1)
                dtid.append(d2)

            elif cell < k:
                dtid.append(cell+n)
        d.append(dtid)

    return d



def main():
    dstart=0
    p=0.3
    n=2
    k=4
    q=0.1
    hours = 10

    #Här sker själva simueringen
                 #(dstart,p,n,k,q):
    d = simGenProt(0,0.3,2,3,1.0)

    #Analysera d

    #Här räknar vi antal protein just nu bara anpassat för k=3
    #antal defekta protein i sista generationen
    print(len(d[10]))
    prot = [0]*6
    for defe in d[10]:
        if defe == 0:
            prot[0]+=1
        elif defe == 1:
            prot[1]+=1
        elif defe == 2:
            prot[2]+=1
        elif defe == 3:
            prot[3]+=1
        elif defe == 4:
            prot[4]+=1
        else:
            prot[5]+=1

    #Skriver ut hur mnga av varje celltyp och sista generationen
    print(prot)
    print(d[10])


    #print("\n",d)
    #print("\n",d)


main()
