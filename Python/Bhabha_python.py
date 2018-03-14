"""
Created on Sat Feb 17 12:17:42 2018
@autor: Fabio Kopp , fabio.kopp@ufrgs.br, departamento de fisica 
"""

import numpy as np
import scipy.integrate as integrate
import matplotlib.pyplot as plt
from scipy.stats import chi2


data1 = np.loadtxt("experimental-data.txt")
data2=  np.loadtxt("sigtotalexp.dat")
colc0 = data1[:,0] #cos(x)

def dsigma(x,s2):
	alphaqed2=np.power((1.0/137.035999074),2)
	s=np.power(s2,2)
	cf=0.389*np.power(10,6)
	FACTOR1=(np.power((1.0 + x),2) + 4.0)/(np.power((1.0-x),2))
	FACTOR2=(1.0+np.power(x,2))/2.0
	FACTOR3=-(np.power((1.0+x),2))/(1.0-x)	
	DSIGMAR=(alphaqed2/(2.0*(s)))*(FACTOR1+FACTOR2+\
	FACTOR3)*cf
	return DSIGMAR

def dsigm(x,s2):
    res=dsigma(x,s2)*(2.0*np.pi*(np.sqrt( 1.0 - np.power(x,2))))
    return res
def sigtotal(s2):
    res=integrate.quad(dsigm, -0.84,0.84,args=(s2))
    return res[0]

def chisquare(s2,v):
    res=0.0
    s=np.power(s2,2)
    n=len(data1[:,0])
    #Para testar com outros valores de energia de centro de massa s2, modificar
    # as colunas referentes a y e yerror de data1.
    x=data1[:,0]
    y=data1[:,9]/s
    d2y=np.sqrt(np.power(data1[:,10]/s,2))
    yth= dsigma(x,s2)
    res= np.sum(np.power( ((y-yth)/d2y),2))
    ndf=n-v
    resr=res/ndf
    chiprob=(1-chi2.cdf(res,ndf))*100.0
    print("                                       ")
    print("                                       ")
    print("---------Analise estatistica-----------")
    
    print("Chi-quadrado = " , res )
    print("Chi-quadrado reduzido = " , resr)
    print("Numero de graus de liberdade(ndf) = ", ndf)
    print("Numero de variaveis ajustadas = ", v)
    print("A P(chi^2,ndf) =  ", chiprob,  "  %")
    print("---------------------------------------")
    
    
x=[]
y=[]
ecm=float(input("Digite a energia de centro de massa (GeV)   "))
for xcos in np.arange(-0.8,0.8,0.02):
     x.append(xcos)
     y.append(dsigma(xcos,ecm))
     #print (i , dsigma(i,34.8))


x1=[]
y1=[]
for s2 in np.arange(10,50,0.1):
    x1.append(s2)
    y1.append(sigtotal(s2))
    #print(s2, sigtotal(s2))
    


#Avaliando a teoria vs dados experimentais
#         s2: energia de centro de massa,  numero de variaveis ajustadas    
chisquare(34.8,0.0)



#gerar graficos
#grafico 1
plt.figure(figsize=(7.5,6.0))
plt.yscale('log', nonposy='clip')
plt.title(r'Espalhamento Bhabha: secao de choque diferencial',fontsize=15)
plt.xlabel(r'$cos(\theta)$',fontsize=15)
plt.ylabel(r'$\frac{d\sigma }{d\Omega} \ [nb/sterad]$',fontsize=15)
plt.xticks(np.arange(-1.0,1.0,0.1))
#plt.yticks(np.arange(0.001,10.0,0.1))
plt.text(-0.2,0.002, r' Dados obtidos do experimento TASSO ', fontsize=12)
#plt.ticklabel_format(style='sci', axis='x', scilimits=(0,0))
#plt.xlim(10,10E15)
plt.ylim(0.001,10.0,0.05)
plt.plot(x,y,'red',label="Teoria (QED)",linewidth=2.5)
plt.errorbar(colc0,data1[:,3]/np.power(14.0,2),xerr=0,yerr=\
             data1[:,4]/np.power(14.0,2),linestyle='', fmt=".",\
             linewidth=1.0,label=r'$\sqrt{s}=14.0 \ GeV$')
plt.errorbar(colc0,data1[:,6]/np.power(22.0,2),xerr=0,yerr=\
             data1[:,7]/np.power(22.0,2),linestyle='', fmt=".",\
             linewidth=1.0,label=r'$\sqrt{s}=22.0 \ GeV$')
plt.errorbar(colc0,data1[:,9]/np.power(34.8,2),xerr=0,yerr=\
             data1[:,10]/np.power(34.8,2),linestyle='', fmt=".",\
             linewidth=1.0,label=r'$\sqrt{s}=34.8 \ GeV$')
plt.errorbar(colc0,data1[:,12]/np.power(38.3,2),xerr=0,yerr=\
             data1[:,13]/np.power(38.3,2),linestyle='', fmt=".",\
             linewidth=1.0,label=r'$\sqrt{s}=38.3 \ GeV$')

plt.errorbar(colc0,data1[:,15]/np.power(43.6,2),xerr=0,yerr=\
             data1[:,16]/np.power(43.6,2),linestyle='', fmt=".",\
             linewidth=1.0,label=r'$\sqrt{s}=43.6 \ GeV$')
#plt.grid(True)
plt.legend(loc="upper left")
plt.legend()
plt.savefig('plot1.png', dpi = 300)
plt.show()
plt.close()



#grafico 2
plt.figure(figsize=(6,5.5))
plt.title(r'Espalhamento Bhabha: secao de choque total',fontsize=15)
plt.xlabel(r'$\sqrt{s} \ [GeV]$',fontsize=15)
plt.ylabel(r'$\sigma \ [nb]$',fontsize=15)
plt.legend(loc='upper right')
plt.xticks(np.arange(10,55,5.0))
#plt.yticks(np.arange(0.001,10.0,0.1))
#plt.text(0.02,1.3, r'an equation: $E=mc^2$', fontsize=15)
#plt.ticklabel_format(style='sci', axis='x', scilimits=(10,0))
#plt.xlim(10,10E15)
#plt.ylim(0.001,10.0,0.1)
plt.plot(x1,y1,'red',label="Teoria (QED)",linewidth=2.5)
plt.errorbar(data2[:,0],data2[:,1],xerr=0,yerr=data2[:,2],\
	linestyle='', fmt=".",c= 'black',linewidth=1.0,label="Dados: Experimento TASSO")
plt.legend()
plt.savefig('plot2.png', dpi = 300)
plt.show()
plt.close()
