import numpy as np
import numba
from scipy.integrate import ode
import sys
import time
import progressbar
import os
from ggfuncs import *
from pyqtgraph import PlotWidget
import pyqtgraph as pq

from PyQt5.QtWidgets import (QWidget, QMainWindow, QTextEdit, 
    QAction, QFileDialog, QApplication, QPushButton,QCheckBox,
    QProgressBar,QGridLayout,QLabel,QLineEdit,QToolTip,
    QMessageBox,QPlainTextEdit)
import PyQt5.QtWidgets
from PyQt5.QtGui import QIcon,QFont,QPixmap  
from PyQt5.QtCore import Qt, pyqtSignal, QObject,QBasicTimer
from PyQt5 import QtCore


class elements_ggchempy(object):
    def __init__(self):
        self.spec = []
        self.mass = []
        self.units = {'spec':'', 'mass':'amu'}
    def print_attributes(self):
        print("--------------- gas -----------------------------") 
        for attribute, value in self.__dict__.items():
            if attribute!='units': 
                print("%15s %2s %12.5g %15s" %(attribute, '=', value, self.units[attribute]))
                
class species_ggchempy(object):
    def __init__(self):
        self.spec = []
        self.mass = []
        self.abun = []
        self.numb = []
        self.charge = []
        self.freq = []
        self.Ebind =[]
        self.Ediff = []
        self.atoms = []
        self.Natoms = []
        self.idx_CO = -1## -1 mean no index.
        self.idx_H2 = -1
        self.idx = {}
        self.units = {'spec':'', 'mass':'amu', 'abun':'', 'numb':'cm^-3', 'charge':''}
    def print_attributes(self):
        print("--------------- gas -----------------------------") 
        for attribute, value in self.__dict__.items():
            if attribute!='units': 
                print("%15s %2s %12.5g %15s" %(attribute, '=', value, self.units[attribute]))

class reactions_ggchempy(object):
    def __init__(self):
        self.spec = []
        self.itype = []
        self.a = []
        self.b = []
        self.c = []
        self.rc = []
        self.idx = []
        self.units = {'spec':'', 'mass':'amu', 'abun':'', 'numb':'cm^-3'}
    def print_attributes(self):
        print("--------------- gas -----------------------------") 
        for attribute, value in self.__dict__.items():
            if attribute!='units': 
                print("%15s %2s %12.5g %15s" %(attribute, '=', value, self.units[attribute]))
                                
class gas_ggchempy(object):
    def __init__(self):
        self.nH = 2e4
        self.Av = 10.0
        self.Chi= 1.0
        self.Zeta = 1.3e-17
        self.T  = 10.0
        self.nH0=3000.0
        self.Av0=2.0
        self.units = {'nH':'cm^-3', 'Av':'mag', 'Chi':'Chi0', 'Zeta':'s^-1', 'T':'K'}
    def print_attributes(self, tofile):
        self.tofile=tofile
        if self.tofile==None:
            print("--------------- gas -----------------------------") 
            for attribute, value in self.__dict__.items():
                if attribute not in ['units','tofile','nH0','Av0']: 
                    print("%15s %2s %12.5g %15s" %(attribute, '=', value, self.units[attribute]))
        else:
            print("--------------- gas -----------------------------", file=self.tofile) 
            for attribute, value in self.__dict__.items():
                if attribute not in ['units','tofile','nH0','Av0']: 
                    print("%15s %2s %12.5g %15s" %(attribute, '=', value, self.units[attribute]), file=self.tofile)

class surface_ggchempy(object):
    Rdb = 0.5
class mantle_ggchempy(object):
    Rdb = 0.5
class dust_ggchempy(object):
    def __init__(self):
        self.radius = 1.0e-5
        self.rho = 3.0
        self.site_density= 1.5e+15
        self.T  = 10.0
        self.nd = 0.0
        self.mass = 0.0
        self.surface = surface_ggchempy()
        self.mantle = mantle_ggchempy()
        self.nbins = 1
        self.stick0 = 1.0
        self.Tpeak = 70.0
        self.fc = 3.0e-19
        self.site_number=0.0
        self.units = {'radius':'cm', 'rho':'g cm^-3', 'site_density':'sites cm^-2', 'T':'K',\
                      'nd':'cm^-3', 'mass':'g', 'stick0':'', 'Tpeak':"K", 'fc':''}
    def print_attributes(self, tofile=None):
        self.tofile=tofile
        if self.tofile==None:
            print("--------------- dust ----------------------------") 
            for attribute, value in self.__dict__.items():
                if attribute not in ['units','surface','mantle','nbins','site_number','tofile']: 
                    print("%15s %2s %12.5g %15s" %(attribute, '=', value, self.units[attribute]))
        else:
            print("--------------- dust ----------------------------", file=self.tofile) 
            for attribute, value in self.__dict__.items():
                if attribute not in ['units','surface','mantle','nbins','site_number','tofile']: 
                    print("%15s %2s %12.5g %15s" %(attribute, '=', value, self.units[attribute]), file=self.tofile)

class ggpars_ggchempy(object):
    def __init__(self):
        self.ggfiles= ['in/network2.txt','in/ed.txt','in/iabun.txt']
        self.d2gmr = 0.01
        self.ne = 12
        self.ns = 655
        self.nr = 6067
        self.nrnp = 8     # 3 + 5
        self.nstr = 8     # string length
        self.nt = 100
        self.ti = 1.0
        self.tf = 1.0e+8
        self.yr_sec = 365.0*24.0*3600.0
        self.atol = 1.0e-20
        self.rtol = 1.0e-6
        self.pi = 3.1415926
        self.kb =1.38054e-16
        self.mp =1.66054e-24
        self.G  = 6.6720e-8 
        self.aa=0.01  # when iswitch.iNTD==1
        self.EM=100.0 # when iswitch.iNTD==2
        self.outdir='out/'
        self.cputime = 0.0
        self.iprogress=0
    def print_attributes(self, tofile=None):
        self.tofile=tofile
        if self.tofile==None:
            print("--------------- ggpars --------------------------") 
            for attribute, value in self.__dict__.items():
                if attribute=='ggfiles':
                    print('Infiles:', self.ggfiles[0], self.ggfiles[1], self.ggfiles[2])
                if attribute in ['ne','ns','nr','d2gmr','tf', 'atol', 'rtol', 'outdir']:
                    if attribute=='outdir':
                        print("%15s %2s %12s" %(attribute, '=', value))
                    else:
                        print("%15s %2s %12.5g" %(attribute, '=', value))
        else:
            print("--------------- ggpars --------------------------", file=self.tofile) 
            for attribute, value in self.__dict__.items():
                if attribute=='ggfiles':
                    print('Infiles:', self.ggfiles[0], self.ggfiles[1], self.ggfiles[2], file=self.tofile)
                if attribute in ['ne','ns','nr','d2gmr','tf', 'atol', 'rtol', 'outdir']:
                    if attribute=='outdir':
                        print("%15s %2s %12s" %(attribute, '=', value), file=self.tofile)
                    else:
                        print("%15s %2s %12.5g" %(attribute, '=', value), file=self.tofile)

class iswitch_ggchempy(object):
    def __init__(self):
        self.iSS=0
        self.iNTD=0
        self.icollapse=0
    
class ggchempy(object):
    def __init__(self):
        self.species = species_ggchempy()
        self.elements = elements_ggchempy()
        self.reactions = reactions_ggchempy()
        self.gas = gas_ggchempy()
        self.dust = dust_ggchempy()
        self.ggpars = ggpars_ggchempy()
        self.iswitch = iswitch_ggchempy()
        self.tyr = []
        self.ts =  []
        self.modelname = ""
        self.res={}
        
    def init_ggchem(self):
        
        self.dust.mass = (4.0/3.0) * self.ggpars.pi * self.dust.radius**3 * self.dust.rho
        self.dust.nd = self.gas.nH * self.ggpars.mp * self.ggpars.d2gmr / self.dust.mass
        self.dust.site_number = 4.0 * self.ggpars.pi * self.dust.radius**2 * self.dust.site_density
        
        network = open(self.ggpars.ggfiles[0],'r').readlines()
        fEbind = self.ggpars.ggfiles[1]
        fiabun = self.ggpars.ggfiles[2]
    
        ### get number of elements, species and reactions:
        L0 = network[0].split()
        ne, ns, nr = L0[0:3]
        self.ggpars.ne, self.ggpars.ns, self.ggpars.nr = int(ne), int(ns), int(nr)
        
        
        ###self.ggpars.print_attributes()
        self.gas.print_attributes(tofile=None)
        self.dust.print_attributes(tofile=None)
        self.ggpars.print_attributes(tofile=None)
        print('Reading network...')
        
        self.elements.spec = np.zeros( self.ggpars.ne, dtype='U'+str(self.ggpars.nstr))
        self.elements.mass = np.zeros( self.ggpars.ne)
        
        self.species.spec = np.zeros( self.ggpars.ns, dtype='U'+str(self.ggpars.nstr))
        self.species.mass = np.zeros( self.ggpars.ns)
        self.species.abun = np.zeros( self.ggpars.ns)
        self.species.numb = np.zeros( self.ggpars.ns)
        self.species.charge = np.zeros( self.ggpars.ns)
        self.species.freq   = np.zeros(self.ggpars.ns)
        self.species.Ebind  = np.zeros(self.ggpars.ns)
        self.species.Ediff  = np.zeros(self.ggpars.ns)
        self.species.atoms  = np.zeros((self.ggpars.ns, self.ggpars.ne))
        self.species.Natoms = np.zeros(self.ggpars.ns)
        
        self.reactions.spec  = np.zeros( (self.ggpars.nr, self.ggpars.nrnp), dtype='U'+str(self.ggpars.nstr))
        self.reactions.itype = np.zeros(self.ggpars.nr, dtype='int')
        self.reactions.a     = np.zeros(self.ggpars.nr)
        self.reactions.b     = np.zeros(self.ggpars.nr)
        self.reactions.c     = np.zeros(self.ggpars.nr)
        self.reactions.rc    = np.zeros(self.ggpars.nr)
        self.reactions.idx   = np.zeros( (self.ggpars.nr, self.ggpars.nrnp), dtype='int')
        
        ### get elements and mass of elements (amu):
        for i in range(1, self.ggpars.ne+1, 1):
            Lx = network[i].split()
            self.elements.spec[i-1]=Lx[0]
            self.elements.mass[i-1]=float(Lx[1])
            
        ### get species and compute species mass (amu):    
        ii=0
        for i in range(self.ggpars.ne+1+1, self.ggpars.ne+1+1+self.ggpars.ns, 1):
            Lx = network[i].split()
            spx = Lx[0].strip()
            
            spmassx = 0.0; Natomsx=0.0
            for j in range(0, self.ggpars.ne, 1):
                self.species.atoms[ii,j]  = float(Lx[j+2])
                spmassx += self.elements.mass[j]*self.species.atoms[ii][j]
                Natomsx += self.species.atoms[ii][j]
            if spx=='E': spmassx=1.0/1832.0
            if spx[0:1]=='G': spmassx=4.0/3.0*self.ggpars.pi*self.dust.radius**3*self.dust.rho/self.ggpars.mp
            
            self.species.spec[ii]   = spx 
            if spx == 'CO': self.species.idx_CO = ii
            if spx == 'H2': self.species.idx_H2 = ii
            self.species.charge[ii] = float(Lx[1])
            self.species.mass[ii]   = spmassx
            self.species.Natoms[ii] = Natomsx
            self.species.idx[spx]=ii
            #if '+' in species.spec[ii] or '-' in species.spec[ii]:
            #    print(species.spec[ii],species.charge[ii])
            
            ii+=1
            
        ### get reactions and get indices:
        nrnp=self.ggpars.nrnp
        nstr=self.ggpars.nstr
        
        ii=0
        for i in range(self.ggpars.ne+self.ggpars.ns+1+1+1, self.ggpars.ne+1+1+1+self.ggpars.ns+self.ggpars.nr, 1):
            Lx = network[i]
            
            Latom = np.zeros(self.ggpars.ne,dtype='int')
            Ratom = np.zeros(self.ggpars.ne,dtype='int')
            Lcharge, Rcharge = 0, 0
            
            for j in range(0, nrnp, 1):
                reac = Lx[j*nstr : (j+1)*nstr].strip()
                self.reactions.spec[ii,j]=reac.strip()
                for k in range(0, self.ggpars.ns, 1):
                    if reac.strip()==self.species.spec[k].strip():
                        #print(reac.strip(),species.spec[k])
                        self.reactions.idx[ii][j]=k
                        ### check balance:
                        if j<=2:
                            for m in range(0, self.ggpars.ne, 1):
                                Latom[m]= Latom[m]+ self.species.atoms[k][m]
                                Lcharge = Lcharge + self.species.charge[k]
                        else:
                            for m in range(0, self.ggpars.ne, 1):
                                Ratom[m]= Ratom[m]+ self.species.atoms[k][m]
                                Rcharge = Rcharge + self.species.charge[k]
                if reac=='': self.reactions.idx[ii][j]=-1
                
            for m in range(0, self.ggpars.ne, 1):
                if Latom[m]!=Ratom[m]:
                    print(ii,Latom)
                    print(ii,Ratom)
                    print('Un-balance of atoms of reaction',ii)
                    sys.exit()
                if Lcharge!=Rcharge:
                    print('Un-balance of electrons in raction',ii,':',Lcharge,Rcharge)
                    sys.exit()
                            
            self.reactions.a[ii], self.reactions.b[ii], self.reactions.c[ii], self.reactions.itype[ii] =\
                float(Lx[64:73]), float(Lx[73:82]),float(Lx[82:91]), int(Lx[91:93])
            ii+=1
            
        ### get Binding energies:
        #Ebind, Ediff, Freq = np.zeros(ns), np.zeros(ns), np.zeros(ns)
        specx  = np.loadtxt(fEbind, usecols=[0], skiprows=2, dtype='str')
        Ebindx = np.loadtxt(fEbind, usecols=[1], skiprows=2)
        nx = np.size(specx)
        for i in range(0, self.ggpars.ns, 1):
            for j in range(0, nx, 1):
                if self.species.spec[i].strip()=='J'+specx[j].strip():
                    self.species.Ebind[i]= Ebindx[j]
                    self.species.Ediff[i]= Ebindx[j]*self.dust.surface.Rdb
                    self.species.freq[i] = np.sqrt(2.0*self.dust.site_density*self.ggpars.kb*self.species.Ebind[i]/self.ggpars.pi**2/self.species.mass[i]/self.ggpars.mp)
                ## For gas-phase species. 
                ## For use in non-thermal desportion when iswitch.iNTD!=0 (iTYPE=14)
                if self.species.spec[i].strip()==specx[j].strip():
                    self.species.Ebind[i]= Ebindx[j]    
                    
        ### get initial abundance:
        L= open(fiabun,'r').readlines()
        nx = int(L[0].split()[0])
        specx  = np.loadtxt(fiabun, usecols=[0], skiprows=1, dtype='str')
        iabunx = np.loadtxt(fiabun, usecols=[1], skiprows=1)
        specx = specx[0:nx]
        iabunx = iabunx[0:nx]
        
        
        for i in range(0, self.ggpars.ns, 1):
            for j in range(0, nx, 1):
                if self.species.spec[i].strip()==specx[j].strip():   
                    self.species.abun[i] = iabunx[j]
                    self.species.numb[i] = iabunx[j]*self.gas.nH
            if self.species.spec[i].strip()=='G0' or self.species.spec[i].strip()=='GRAIN0': 
                if (self.species.abun[i]==0.0 and self.species.numb[i]==0.0):
                    self.species.abun[i] = self.dust.nd/self.gas.nH
                    self.species.numb[i] = self.dust.nd
        return
    
    def compute_reaction_rate_coefficients(self):
        #print('---compute_reaction_rate_coefficients(): gas.nH',gas.nH)
        for i in range(0, self.ggpars.nr, 1):
            ### determine dust grains:
            if self.dust.nbins<=1:
                rdx = self.dust.radius
                ndx = self.dust.nd
            else:
                idust = 0 ###???
                rdx = rds[idust]
                ndx = nds[idust]
                
            if self.reactions.itype[i]==0:
                ### ion -- dust
                self.reactions.rc[i] = self.reactions.a[i]*(self.gas.T/300.0)**self.reactions.b[i] * self.gas.nH/ndx
                
            elif self.reactions.itype[i]==1:
                ### cosmic ray ionization
                self.reactions.rc[i] = self.reactions.a[i]*self.gas.Zeta
                
            elif 2<=self.reactions.itype[i]<=12:
                ### 2 --- ion-molecule reaction and charge exchange reaction
                ### 3 --- negative ion - neutral species reaction
                ### 4 --- radiative association : ion-neutral
                ### 5 --- electron ejection : negative ion – neutral
                ### 6 --- neutral + neutral --> ion + electron
                ### 7 --- neutral + neutral chemical reaction
                ### 8 --- neutral + neutral radiative association
                ### 9 --- dissociative recombination
                ### 10--- radiative recombination (no dissociation)
                ### 11--- cation anion recombination
                ### 12--- electron attachment
                self.reactions.rc[i] = self.reactions.a[i]*(self.gas.T/300.0)**self.reactions.b[i]*np.exp(-self.reactions.c[i]/self.gas.T)
                
            elif self.reactions.itype[i]==13:
                ### 13--- photo-ionization and photo-dissociation
                if self.iswitch.iSS==1:
                    specx = self.reactions.spec[i,0].strip()
                    if specx=='H2' or specx=='CO':
                        xh2 = self.species.abun[self.species.idx_H2]
                        xco = self.species.abun[self.species.idx_CO]
                        xk  = self.gas.Chi*self.reactions.a[i]*np.exp(-self.reactions.c[i]*self.gas.Av)
                        self.reactions.rc[i] = H2_CO_self_shielding(specx, self.gas.Av, xh2, xco, xk)
                        #print(specx,'self-shielding is considered: rc=',reactions.rc[i],xco,xh2)
                    else:
                        self.reactions.rc[i] = self.gas.Chi*self.reactions.a[i]*np.exp(-self.reactions.c[i]*self.gas.Av)
                else:
                    self.reactions.rc[i] = self.gas.Chi*self.reactions.a[i]*np.exp(-self.reactions.c[i]*self.gas.Av)
                
                
            elif self.reactions.itype[i]==14:
                ### 14--- surface reaction
                k1,k2,k3,k4=self.reactions.idx[i,0], self.reactions.idx[i,1], self.reactions.idx[i,2], self.reactions.idx[i,3]
                if (k1*k2==0): 
                    print('## surface reaction (itype=14): no match!!!',k1,k2)
                    sys.exit()
                pro=0.0
                if self.reactions.c[i]==0.0:
                    if (k1==k2): pro=0.5
                    if (k1!=k2): pro=1.0
                else:
                    pro = np.exp(-self.reactions.c[i]/self.dust.T)
                rdiff1 = self.species.freq[k1]*np.exp(-self.species.Ediff[k1]/self.dust.T)/self.dust.site_number
                rdiff2 = self.species.freq[k2]*np.exp(-self.species.Ediff[k2]/self.dust.T)/self.dust.site_number
                
                self.reactions.rc[i] = self.reactions.a[i]*pro*(rdiff1+rdiff2)/ndx ##dust.nd
                
                ## for network do not consider reaction desorption
                if 'J' not in self.species.spec[k4]: 
                    self.reactions.rc[i]=0.0
                
                ## when consider reactive desorption (b in KJ/mol):
                if self.iswitch.iNTD!=0 and self.reactions.b[i]<0.0:
                    ## JA + JB -> C
                    fx=0.0  ## fraction
                    if 'J' not in self.reactions.spec[i,3] and self.reactions.spec[i,4].strip()=="":
                        if self.iswitch.iNTD==1:
                            ## Garrod et al., 2007
                            
                            fx = reactive_desorption_2007(self.reactions.b[i], \
                                                            self.species.Natoms[k4],\
                                                            self.species.Ebind[k4],
                                                            aa=self.ggpars.aa)
                            #print(fx,i,reactions.spec[i,0:5],reactions.b[i],species.Natoms[k4],species.Ebind[k4])
                        elif self.iswitch.iNTD==2:
                            ## Minissale et al., A&A 585, A24 (2016)
                            enthalpy = np.abs(self.reactions.b[i] * 120.274) ## KJ/mol->K
                            fx = reactive_desorption_2016(enthalpy,\
                                                            self.species.Natoms[k4],\
                                                            self.species.mass[k4],\
                                                            self.species.Ebind[k4],self.ggpars.EM)
                        
                        ## times fraction to the rate coefficient from last reaction (JA + JB -> JC).
                        self.reactions.rc[i]   = fx*self.reactions.rc[i-1]       ## current reaction
                        self.reactions.rc[i-1] = (1.0-fx)*self.reactions.rc[i-1] ## last reaction
                        #print(fx,reactions.rc[i-1],reactions.rc[i])
               
            elif self.reactions.itype[i]==15:
                ### 15--- thermal evaporation
                self.reactions.rc[i] = self.reactions.a[i]*self.species.freq[self.reactions.idx[i,0]]*np.exp(-self.species.Ebind[self.reactions.idx[i,0]]/self.dust.T)
                
            elif self.reactions.itype[i]==16:
                ### 16--- cosmic ray induced general desorption
                self.reactions.rc[i] = self.dust.fc*self.reactions.a[i]*self.species.freq[self.reactions.idx[i,0]]*np.exp(-self.species.Ebind[self.reactions.idx[i,0]]/self.dust.Tpeak)
                
            elif 17<=self.reactions.itype[i]<=18: 
                ### 17--- photodissociation by cosmic rays
                ###      (JAB -> JA + JB corresponds to AB -> A + B)
                ### 18--- photodissociation by cosmic rays
                ###      (JAB -> JA + JB corresponds to AB -> AB+ + E)
                self.reactions.rc[i] = self.reactions.a[i]*self.gas.Zeta
        
            elif 19<=self.reactions.itype[i]<=20:
                ### 19--- photodissociation by background photons
                ###      (JAB -> JA + JB corresponds to AB -> A + B)
                ### 20--- photodissociation by background photons
                ###      (JAB -> JA + JB corresponds to AB -> AB+ + E)
                self.reactions.rc[i] = self.gas.Chi*self.reactions.a[i]*np.exp(-self.reactions.c[i]*self.gas.Av)
                
            elif self.reactions.itype[i]==99:
                ### 99--- adsorption on grain surface
                vth = np.sqrt(8.0*self.ggpars.kb*self.gas.T/self.ggpars.pi/self.species.mass[self.reactions.idx[i,0]]/self.ggpars.mp)
                self.reactions.rc[i] = self.dust.stick0*self.ggpars.pi*rdx**2*vth*ndx
        
        
        
    def run(self, modelname):
        """
        The main code to run
        """
        self.modelname = modelname
        self.init_ggchem()
        self.compute_reaction_rate_coefficients()
        
        ns = self.ggpars.ns
        nr = self.ggpars.nr
        rc = self.reactions.rc
        idx = self.reactions.idx
        
        @numba.jit(nopython=True )   
        def fode(t,y):
            """Function for scipy.integrate.ode:"""
            ydot = np.zeros(ns)
            up, down = np.zeros(ns),np.zeros(ns)
            for i in range(0,nr,1):
                IR1, IR2, IR3, IP1, IP2, IP3, IP4, IP5 = idx[i,:]
                term = 0.0
                if (IR1!=-1 and IR2==-1 and IR3==-1): term = rc[i]*y[IR1]
                if (IR1!=-1 and IR2!=-1 and IR3==-1): term = rc[i]*y[IR1]*y[IR2]
                if (IR1!=-1 and IR2!=-1 and IR3!=-1): term = rc[i]*y[IR1]*y[IR2]*y[IR3]
                ###if(term<0.0): term=0.0
                
                if (IP1!=-1): up[IP1]=up[IP1]+term
                if (IP2!=-1): up[IP2]=up[IP2]+term
                if (IP3!=-1): up[IP3]=up[IP3]+term
                if (IP4!=-1): up[IP4]=up[IP4]+term
                if (IP5!=-1): up[IP5]=up[IP5]+term
                            
                if (IR1!=-1): down[IR1]=down[IR1]+term
                if (IR2!=-1): down[IR2]=down[IR2]+term
                if (IR3!=-1): down[IR3]=down[IR3]+term
            ydot = up-down
            return ydot
        
        @numba.jit(nopython=True)
        def fjac(t,y):
            """Jacobian Function for scipy.integrate.ode:"""
            pd = np.zeros((ns,ns))
            for i in range(0,nr,1):
                IR1, IR2, IR3, IP1, IP2, IP3, IP4, IP5 = idx[i,:]
                term,term1,term2= 0.0,0.0,0.0
                
                if (IR2==-1 and IR3==-1):
                    ### no second reantant:
                    term = rc[i]
                    pd[IR1][IR1] = pd[IR1][IR1] -term
                    pd[IP1][IR1] = pd[IP1][IR1] +term
                    if (IP2!=-1): pd[IP2][IR1] = pd[IP2][IR1] +term
                    if (IP3!=-1): pd[IP3][IR1] = pd[IP3][IR1] +term
                    if (IP4!=-1): pd[IP4][IR1] = pd[IP4][IR1] +term
                    if (IP5!=-1): pd[IP5][IR1] = pd[IP5][IR1] +term
                else:
                    ### with  second reantant:
                    term1 = rc[i]*y[IR1]
                    term2 = rc[i]*y[IR2]
                    pd[IR1][IR1] = pd[IR1][IR1] - term2
                    pd[IR1][IR2] = pd[IR1][IR2] - term1
                    pd[IR2][IR1] = pd[IR2][IR1] - term2
                    pd[IR2][IR2] = pd[IR2][IR2] - term1
                    pd[IP1][IR1] = pd[IP1][IR1] + term2
                    pd[IP1][IR2] = pd[IP1][IR2] + term1
                    if (IP2!=-1): pd[IP2][IR1] = pd[IP2][IR1] + term2
                    if (IP2!=-1): pd[IP2][IR2] = pd[IP2][IR2] + term1
                    if (IP3!=-1): pd[IP3][IR1] = pd[IP3][IR1] + term2
                    if (IP3!=-1): pd[IP3][IR2] = pd[IP3][IR2] + term1
                    if (IP4!=-1): pd[IP4][IR1] = pd[IP4][IR1] + term2
                    if (IP4!=-1): pd[IP4][IR2] = pd[IP4][IR2] + term1
                    if (IP5!=-1): pd[IP5][IR1] = pd[IP5][IR1] + term2
                    if (IP5!=-1): pd[IP5][IR2] = pd[IP5][IR2] + term1
            return pd
            
        #### speed up with numba.jit:
        def wrapper_to_fode(t, y):
            return fode(t, y) 
        def wrapper_to_fjac(t,y):
            return fjac(t,y)
        
        ### open output file:
        if not os.path.isdir(self.ggpars.outdir[0:-1]):
            os.makedirs(self.ggpars.outdir[0:-1])
    
        fileout=open(self.ggpars.outdir+self.modelname+'.dat','w')
        print('Model name:',self.modelname)
        print('----------------- ggchempy 2021 -----------------', file=fileout)
        self.gas.print_attributes(tofile=fileout)
        self.dust.print_attributes(tofile=fileout)
        self.ggpars.print_attributes(tofile=fileout)
        print('--------------- time evolution-------------------', file=fileout)
        
        fmt = '%10s '*(self.ggpars.ns+1)            
        print( fmt%('time',*self.species.spec), file=fileout )
        fmt = '%10.3e '*(self.ggpars.ns+1)
                
        ### time grid
        self.tyr = np.logspace( np.log10(self.ggpars.ti), np.log10(self.ggpars.tf), self.ggpars.nt)
        self.ts  = self.tyr * self.ggpars.yr_sec
        
        ### DO the integration:     
        t0 = 0.0
        y0 = self.species.numb
        r  = ode(wrapper_to_fode, wrapper_to_fjac).set_integrator(\
                 'vode', method='bdf', \
                 rtol=self.ggpars.rtol, atol=self.ggpars.atol, \
                 with_jacobian=True, nsteps=500000)#, ixpr=True)
        r.set_initial_value(y0, t0)
        
        start_time = time.time()
        
        ######## start ODE solver:
        numbt = np.zeros((self.ggpars.nt, self.ggpars.ns))
        
        widgets = [
        '\x1b[33mggchem:\x1b[39m',
        progressbar.Percentage(),
        progressbar.Bar(marker='\x1b[32m#\x1b[39m'),
        ]
        
            
        #loop=tqdm(total=ggpars.nt,position=0,leave=True)
        bar = progressbar.ProgressBar(widgets=widgets, maxval=self.ggpars.nt).start()
        
        it = 0
        while r.successful() and r.t < self.ggpars.tf * self.ggpars.yr_sec:
            self.iprogress = it
            
            ### to update rate coefficients:
            self.compute_reaction_rate_coefficients()
            
            r.integrate(self.ts[it])
            r.y[np.argwhere(r.y<0.0)]=0.0

            numbt[it][:] = r.y
            self.species.numb = r.y
            self.species.abun = r.y/self.gas.nH
            
            nH_save      = np.float(self.gas.nH)
            numb_save    = np.array(r.y)
            
            ########################
            print( fmt%(r.t/self.ggpars.yr_sec, *r.y), file=fileout)
                
            it+=1
            bar.update(it)
            #loop.set_description("ggchem:".format(it))
            #loop.update(1)
            if it==self.ggpars.nt: break
        bar.finish()
        #loop.close()
        
        self.ggpars.cputime = time.time() - start_time
        print("Finished, CPU time = %s seconds ---" % (self.ggpars.cputime) )
        fileout.close()
        
        self.res = {}
        for i in range(0,ns,1):
            self.res[self.species.spec[i]]=numbt[:,i]   ###number density (cm^-3)

        self.res['time']= self.tyr
        self.res['nH']  = self.gas.nH
        self.res['Av']  = self.gas.Av
        return self.res
    
def loadgg(ggfile,NAUTILUS=False):
    """
    *** PURPOSE:
        load existing ggchem model file.
        e.g.: TMC1=loadgg('out/TMC1.dat')
    *** INPUT:
        ggfile --- ggchem model file
        NAUTILUS --- if True, load model file with format of NAUTILUS.
    *** OUTPUT:
        ab --- model results in python dictionary.
    """
    import numpy as np
    
    f=open(ggfile,'r')
    d=f.readlines()
    f.close()
    
    ### count skiprows
    skiprows=0
    for x in d:
        skiprows+=1
        if 'time' in x:
            break
    ### store species
    x=d[skiprows-1].split()
    sp=x[:]
    
    if NAUTILUS==True: 
        skiprows=1
        x=d[skiprows-1].split()
        sp=x[:]
        
    ### load data:
    datax=np.loadtxt(ggfile,unpack=True,skiprows=skiprows,dtype='str')
    ns = np.shape(datax)
    data = np.zeros(ns)
    for i in range(0,ns[0],):
        for j in range(0,ns[1],1):
            if 'e' not in datax[i][j] and 'E' not in datax[i][j]:
                data[i][j] = float(datax[i][j].replace('-','e-'))
            else:
                data[i][j] = float(datax[i][j])
    ab={}
    for i in range(0,np.size(sp),1):
        ab[sp[i]]=data[i][:]
    return ab    

def latex_species(sp):
    """
    *** PURPOSE:
        make a species formula to Latex format:
        e.g. H3O+   =>   r'${\rm H_3O^+}$'
    *** INPUT:
        sp --- species in string. e.g. 'H2O'
    *** OUTPUT:
        xsp --- species in latex format
    """
    xsp=r'${\rm '
    for x in sp:
        if x in '23456789':
            x='_'+x
        elif x in '+-':
            x='^{'+x+'}'
        xsp+=x
    xsp+='}$'
    if '[1_3C]' in xsp:
        xsp=xsp.replace('[1_3C]','^{13}C')
    return xsp


###############################################################

ggchempy = ggchempy()

class GGCHEMGUI(QMainWindow):

    def __init__(self):
        super().__init__()

        self.initUI()


    def initUI(self):   
        
        QToolTip.setFont(QFont('SansSerif', 10))
        self.setToolTip('GGCHEM v1.0') 
        
        ### set common font:
        self.setFont(QFont('SansSerif', 10))
        self.statusBar()
        key_bg_color='#ff6600'

        
        label_bg_color="#ffc266"#"#ffe0cc" #'#e6e6e6'
        
        ### Input files:
        self.functions = QLabel(self)
        self.functions.setGeometry(0, 39, 350, 20)
        self.functions.setAutoFillBackground(True)
        self.functions.setStyleSheet("background-color: "+label_bg_color+";  border: 1px solid "+label_bg_color) 
        self.functions.setText('Functions:')
        self.functions.setToolTip('please select your physical and/or chemical functions.')
        
        
        ### Check boxes to select parameters:
        self.cb_static = QCheckBox('static mode', self)
        #grid.addWidget(self.cb_static)
        #cb_static.move(20, 60)
        self.cb_static.setChecked(False)
        self.cb_static.setGeometry(QtCore.QRect(20, 60, 200, 20))
        self.cb_static.toggle()
        self.cb_static.stateChanged.connect(self.change_to_static_model)
        self.cb_static.setToolTip('A model with fixed physical parameters listed below.')
        
        self.cb_collap = QCheckBox('collapse mode', self)
        #grid.addWidget(self.cb_collap)
        #cb_collap.move(20, 80)
        self.cb_collap.setChecked(True)
        self.cb_collap.setGeometry(QtCore.QRect(20, 80, 200, 20))
        self.cb_collap.toggle()
        self.cb_collap.stateChanged.connect(self.change_to_collap_model)
        self.cb_collap.setToolTip('A model with varying physical parameters of a collapsing core. \
        \nIn this mode, nH should be set to: "initial density(nH0), final density". e.g. "3000 , 1e7" \
        \nAv (extinction) will be used as the base extinction (Av0). e.g. you can set it to 2.0.\
        \n( Av will be updated according to Av = Av0 + (nH/nH0)**(2/3). \
        \nNOTE: It is better to set "Time steps" to be a large number, e.g. >=2000.')
        
        self.cb_ss = QCheckBox('self-sheilding of CO and H2', self)
        #grid.addWidget(self.cb_ss)
        #cb_ss.move(20, 100)
        self.cb_ss.setChecked(False)
        self.cb_ss.setGeometry(QtCore.QRect(20, 100, 250, 20))
        self.cb_ss.toggle()
        self.cb_ss.setToolTip('Ref.: Lee et al., A&A, 311, 690-707')

        
        self.cb_rd_g07 = QCheckBox('reactive desorption of G07', self)
        #grid.addWidget(self.cb_rd_g07)
        #cb_rd_g07.move(20, 120)
        self.cb_rd_g07.setChecked(True)
        self.cb_rd_g07.setGeometry(QtCore.QRect(20, 120, 250, 20))
        self.cb_rd_g07.toggle()
        self.cb_rd_g07.setToolTip('Ref.: Garrod et al., ApJ, 467, 1103 (2017)')
        self.cb_rd_g07.stateChanged.connect(self.change_to_g07)
        
        
        self.cb_rd_m16 = QCheckBox('reactive desorption of M16', self)
        #grid.addWidget(self.cb_rd_m16)
        #cb_rd_m16.move(20, 140)
        self.cb_rd_m16.setChecked(False)
        self.cb_rd_m16.setGeometry(QtCore.QRect(20, 140, 250, 20))
        self.cb_rd_m16.toggle()
        self.cb_rd_m16.setToolTip('Ref.: Minissale et al., A&A 585, A24 (2016)')
        self.cb_rd_m16.stateChanged.connect(self.change_to_m16)
        
        ### Input files:
        self.inputf = QLabel(self)
        self.inputf.setGeometry(0, 180, 350, 20)
        self.inputf.setAutoFillBackground(True)
        self.inputf.setStyleSheet("background-color: "+label_bg_color+";  border: 1px solid "+label_bg_color) 
        self.inputf.setText('Input files:')
        self.inputf.setToolTip('Please type the three necessary input files below.')
        
        xmid = 170
        xwid = 170
        self.inputfnet = QLabel(self)
        self.inputfnet.setGeometry(20, 200, 170, 20)
        self.inputfnet.setText('Reaction network:')
        
        self.finputnet = QLineEdit(self)
        self.finputnet.setGeometry(xmid, 200, xwid, 20)
        self.finputnet.setText('in/network2.txt')
        self.finputnet.setToolTip('file of reaction network in which element and species.  \
                                  \nand reaction lists should be included.')#. For example:\
                                  #\n12 500 5000 #number of elements, species and reactions. #Element list:\
                                  #\nH     1.0\
                                  #\nHE    4.0\
                                  #\n... \
                                  #\n       e  H  He  C  N  O  Si S  Fe Na Mg P  CL # species list:\
                                  #\nH      0  1  0   0  0  0  0  0  0  0  0  0  0\
                                  #\n...\
                                  #\n         #reaction list:\
                                  #\nC     O       CO               alpha beta gamma itype\
                                  #\n...')
        
        #self.finputnet.textChanged[str].connect(self.onChanged)
        
        self.inputfiab = QLabel(self)
        self.inputfiab.setGeometry(20, 220, 170, 20)
        self.inputfiab.setText('Initial abundances:')
        
        self.finputiab = QLineEdit(self)
        self.finputiab.setGeometry(xmid, 220, xwid, 20)
        self.finputiab.setText('in/iabun.txt')
        self.finputiab.setToolTip('file of initial abundances list.')
        #self.finputiab.textChanged[str].connect(self.onChanged)
        
        self.inputfed = QLabel(self)
        self.inputfed.setGeometry(20, 240, 170, 20)
        self.inputfed.setText('Binding energies:')
        
        self.finputed = QLineEdit(self)
        self.finputed.setGeometry(xmid, 240, xwid, 20)
        self.finputed.setText('in/ed.txt')
        self.finputed.setToolTip('file of binding energy list.')
        #self.finputed.textChanged[str].connect(self.onChanged)
        
        
        ### Gas Parameters:
        self.gaspars = QLabel(self)
        self.gaspars.setGeometry(0, 270, 350, 20)
        self.gaspars.setAutoFillBackground(True)
        self.gaspars.setLineWidth(2.0)
        self.gaspars.setText('Gas parameters:')
        self.gaspars.setStyleSheet("background-color: "+label_bg_color+";  border: 1px solid "+label_bg_color) 
        
        xwid_left,xwid = 95,70  ### length for left-QLineEdit
        
        self.gas_nH_label = QLabel(self)
        self.gas_nH_label.setGeometry(20, 290, 100, 20)
        self.gas_nH_label.setText('n<sub>H</sub>=')
        
        self.gas_Av_label = QLabel(self)
        self.gas_Av_label.setGeometry(180, 290, 100, 20)
        self.gas_Av_label.setText('A<sub>V</sub> =')
        
        self.gas_T_label = QLabel(self)
        self.gas_T_label.setGeometry(20, 310, 100, 20)
        self.gas_T_label.setText('T<sub>gas</sub> =')
        
        self.gas_Zeta_label = QLabel(self)
        self.gas_Zeta_label.setGeometry(180, 310, 100, 20)
        self.gas_Zeta_label.setText(u'\u03B6<sub>CR</sub> =')  ##zeta
        
        self.gas_Chi_label = QLabel(self)
        self.gas_Chi_label.setGeometry(20, 330, 100, 20)
        self.gas_Chi_label.setText(u'\u03C7 =')  ##chi
        
        self.gas_nH_value = QLineEdit(self)
        self.gas_nH_value.setGeometry(80, 290, xwid_left, 20)
        self.gas_nH_value.setText('2e+4')
        self.gas_nH_value.setToolTip('gas density. unit: cm^-3\
                                     \n e.g.:\
                                     \n     "1.0e+4" for static mode.\
                                     \n     "3000.0,1e7" for collapse mode.')
    
        self.gas_Av_value = QLineEdit(self)
        self.gas_Av_value.setGeometry(270, 290, xwid, 20)
        self.gas_Av_value.setText('10.0')
        self.gas_Av_value.setToolTip('visual extinction. unit: mag\
                                     \n When "collapse mode" is checked, Av is used as the base extinction.')
        
        self.gas_T_value = QLineEdit(self)
        self.gas_T_value.setGeometry(80, 310, xwid_left, 20)
        self.gas_T_value.setText('10.0')
        self.gas_T_value.setToolTip('gas temperature. unit: K')
    
        self.gas_Zeta_value = QLineEdit(self)
        self.gas_Zeta_value.setGeometry(270, 310, xwid, 20)
        self.gas_Zeta_value.setText('1.3e-17')
        self.gas_Zeta_value.setToolTip('cosmic ionization. unit: s^-1')
        
        self.gas_Chi_value = QLineEdit(self)
        self.gas_Chi_value.setGeometry(80, 330, xwid_left, 20)
        self.gas_Chi_value.setText('1.0')
        self.gas_Chi_value.setToolTip(u'radiation field. unit: \u03C7_0 (Draine 1978)')
        
        self.d2gmr_label = QLabel(self)
        self.d2gmr_label.setGeometry(180, 330, 100, 20)
        self.d2gmr_label.setText(u'd2gmr =')  ##chi
        self.d2gmr_label.setToolTip('dust-to-gas mass ratio')
        self.d2gmr_label.setStyleSheet("background-color: "+label_bg_color+";  border: 1px solid "+label_bg_color) 
        
        self.d2gmr_value = QLineEdit(self)
        self.d2gmr_value.setGeometry(270, 330, xwid, 20)
        self.d2gmr_value.setText('0.01')
        self.d2gmr_value.setToolTip('dust-to-gas mass ratio')
        
        
        ### Dust parameters:
        #_translate = QtCore.QCoreApplication.translate
        #red_color = "<html><head/><body><p><span style=\" color:#cc0000;\">Dust Parameters:</span></p></body></html>"
        self.dustpars = QLabel(self)
        self.dustpars.setGeometry(0, 360, 350, 20)
        self.dustpars.setAutoFillBackground(True)
        self.dustpars.setText('Dust parameters:' )
        self.dustpars.setStyleSheet("background-color: "+label_bg_color+";  border: 1px solid "+label_bg_color) 
        #self.dustpars.setText(_translate(('Dust Parameters:', red_color) )
        
        
        self.dust_T_label = QLabel(self)
        self.dust_T_label.setGeometry(20, 380, 100, 20)
        self.dust_T_label.setText('T<sub>dust</sub>=')
        
        self.dust_R_label = QLabel(self)
        self.dust_R_label.setGeometry(180, 380, 100, 20)
        self.dust_R_label.setText('r<sub>dust</sub>=')
        
        self.dust_rho_label = QLabel(self)
        self.dust_rho_label.setGeometry(20, 400, 100, 20)
        self.dust_rho_label.setText(u'\u03C1<sub>dust</sub> =')
        
        self.dust_Rdb_label = QLabel(self)
        self.dust_Rdb_label.setGeometry(180, 400, 100, 20)
        self.dust_Rdb_label.setText(r'R<sub>db</sub>=')
        
        self.dust_T_value = QLineEdit(self)
        self.dust_T_value.setGeometry(80, 380, xwid_left, 20)
        self.dust_T_value.setText('10.0')
        self.dust_T_value.setToolTip('dust temperature. unit: K')
             
        self.dust_R_value = QLineEdit(self)
        self.dust_R_value.setGeometry(270, 380, xwid, 20)
        self.dust_R_value.setText('1.0e-5')
        self.dust_R_value.setToolTip('dust radius. unit: cm')
        
        self.dust_rho_value = QLineEdit(self)
        self.dust_rho_value.setGeometry(80, 400, xwid_left, 20)
        self.dust_rho_value.setText('3.0')
        self.dust_rho_value.setToolTip('dust mass density. unit: g cm^-3')
             
        self.dust_Rdb_value = QLineEdit(self)
        self.dust_Rdb_value.setGeometry(270, 400, xwid, 20)
        self.dust_Rdb_value.setText('0.77')
        self.dust_Rdb_value.setToolTip('diffusion-to-binding energy ratio.')
        
        ### ODE solver parameters
        self.odepars = QLabel(self)
        self.odepars.setGeometry(0, 430, 350, 20)
        self.odepars.setAutoFillBackground(True)
        self.odepars.setText('ODE solver parameters:' )
        self.odepars.setStyleSheet("background-color: "+label_bg_color+";  border: 1px solid "+label_bg_color) 

        self.ode_atol_label = QLabel(self)
        self.ode_atol_label.setGeometry(20, 450, 100, 20)
        self.ode_atol_label.setText(r'ATOL=')
        
        self.ode_atol_label = QLabel(self)
        self.ode_atol_label.setGeometry(180, 450, 100, 20)
        self.ode_atol_label.setText(r'RTOL=')
        
        self.ode_atol_value = QLineEdit(self)
        self.ode_atol_value.setGeometry(80, 450, xwid_left, 20)
        self.ode_atol_value.setText('1.0e-15')
        self.ode_atol_value.setToolTip('absolute tolerance')
             
        self.ode_rtol_value = QLineEdit(self)
        self.ode_rtol_value.setGeometry(270, 450, xwid, 20)
        self.ode_rtol_value.setText('1.0e-5')
        self.ode_rtol_value.setToolTip('relative tolerance')
        
        
        ### output control
        self.output_control = QLabel(self)
        self.output_control.setGeometry(0, 480, 350, 20)
        self.output_control.setAutoFillBackground(True)
        self.output_control.setText('Output controls:' )
        self.output_control.setStyleSheet("background-color: "+label_bg_color+";  border: 1px solid "+label_bg_color) 
        
        self.output_modelname_label = QLabel(self)
        self.output_modelname_label.setGeometry(20, 500, 100, 20)
        self.output_modelname_label.setText(r'Model name=')
        self.output_modelname_value = QLineEdit(self)
        self.output_modelname_value.setGeometry(120, 500, 220, 20)
        self.output_modelname_value.setText('DC')
        self.output_modelname_value.setToolTip('Simulated model will be saved into "out/Modelname.dat"')
        self.output_modelname_value.textChanged[str].connect(self.ModelnameOnChanged)
        
        self.output_finalage_label = QLabel(self)
        self.output_finalage_label.setGeometry(20, 520, 100, 20)
        self.output_finalage_label.setText(r'Final age=')
        self.output_finalage_value = QLineEdit(self)
        self.output_finalage_value.setGeometry(120, 520, 220, 20)
        self.output_finalage_value.setText('1.0e+7')
        self.output_finalage_value.setToolTip('The final age when the model stops.\
                                            \n ***Will not work when "collapse mode" is checked.')
        
        
        self.output_timestep_label = QLabel(self)
        self.output_timestep_label.setGeometry(20, 540, 100, 20)
        self.output_timestep_label.setText(r'Time steps=')
        self.output_timestep_value = QLineEdit(self)
        self.output_timestep_value.setGeometry(120, 540, 220, 20)
        self.output_timestep_value.setText('100')
        self.output_timestep_value.setToolTip('Total time steps in log-space \
                                              \nat which the simulated results will be saved.\
                                              \nA large number is better for "collapse mode", e.g. >=2000')
        
        ### output information, such as Running, CPUtime. 
        self.output_info = QLabel(self)
        self.output_info.setGeometry(20, 575, 200, 20)
        self.output_info.setAutoFillBackground(True)
        self.output_info.setText('')
        self.output_info.setStyleSheet("background-color: "+key_bg_color+";  border: 1px solid "+key_bg_color) 
        self.output_info.setToolTip('this area is to show the processing information.')

        self.btn_start = QPushButton('run', self)
        #self.btn_start.move(220, 250)
        self.btn_start.setGeometry(QtCore.QRect(220, 575, 100, 20))
        self.btn_start.clicked.connect(self.change_output_info_to_running)
        self.btn_start.clicked.connect(self.runModel)
        self.btn_start.setToolTip('Press it to run and wait for a moment.')
        self.btn_start.setStyleSheet("background-color: cyan;  border: 1px solid "+key_bg_color) 
        
        ### start and progress bar
        self.pbar = QProgressBar(self)
        self.pbar.setGeometry(20, 600, 300, 25)
        self.pbar.setValue(0)
        
        ### Save the key order:
        self.gg_keys = ['icollapse','iSS','iNTD',\
                   'fnet','fiab','fed',\
                   'nH','Av','Tgas','Zeta','Chi',\
                   'd2gmr',\
                   'Tdust','rdust','rhodust','Rdb',\
                   'atol','rtol',\
                   'modelname','finalage','timestep']
        self.n_gg_keys = np.size(self.gg_keys)
        self.gg_dict = {'icollapse':0,'iSS':0,'iNTD':0,\
                   'fnet':'in/network.txt','fiab':'in/iabun.txt','fed':'in/ed.txt',\
                   'nH':2e4,'Av':10.0,'Tgas':10.0,'Zeta':1.3e-17,'Chi':1.0,\
                   'd2gmr':0.01,\
                   'Tdust':10.0,'rdust':1.0e-5,'rhodust':3.0,'Rdb':0.77,\
                   'atol':1.0e-15,'rtol':1.0e-5,\
                   'modelname':'XXX','finalage':1e7,'timestep':100}
                   
        
        self.btn_load_pars = QPushButton('load from file:',self)
        self.btn_load_pars.setGeometry(QtCore.QRect(20, 0, 95, 18))
        self.btn_load_pars.clicked.connect(self.load_pars_from_file)
        self.btn_load_pars.setStyleSheet("background-color: cyan;  border: 1px solid "+key_bg_color) 
        self.btn_load_pars.setToolTip('Load parameters from the following file.')

        self.load_value = QLineEdit(self)
        self.load_value.setGeometry(115, 0, 225, 18)
        self.load_text = 'models/'+ self.output_modelname_value.text().strip() +'.txt'
        self.load_tip = 'A file contains previous parameters:\n'+ self.load_text +\
                        '\n Type your file here and press "load from file".'
        self.load_value.setText( self.load_text )
        self.load_value.setToolTip(self.load_tip)
        
        
        self.btn_save_pars = QPushButton('save into file:',self)
        self.btn_save_pars.setGeometry(QtCore.QRect(20, 20, 95, 18))
        self.btn_save_pars.clicked.connect(self.save_pars_into_file)
        self.btn_save_pars.setStyleSheet("background-color: cyan;  border: 1px solid "+key_bg_color) 
        self.btn_save_pars.setToolTip('Save parameters into the following file.')

        self.save_value = QLineEdit(self)
        self.save_value.setGeometry(115, 20, 225, 18)
        self.save_text = 'models/'+ self.output_modelname_value.text().strip() +'.txt'
        self.save_tip = 'A file to save parameters:\n'+ self.save_text +\
                        '\n Type your file here and press "save into file".'
        self.save_value.setText( self.save_text)
        self.save_value.setToolTip( self.save_tip)
        
        
        # Add PlotWidget control
        #self.plotWidget_ted = PlotWidget(self)
        # Set the size and relative position of the control
        #self.plotWidget_ted.setGeometry(QtCore.QRect(240,80,100,100))
        
        
        window = self.setGeometry(300, 300, 350, 630)
        self.setWindowTitle('GGCHEM v1.0')
        self.setWindowIcon( QIcon('icon/icon.png') )
        #window.setToolTip('A two-phase gas-grain chemical model.')
        #self.setStyleSheet("background-color: gray;  border: 1px solid gray") 
        
        self.show()
    
        
    def collect_info(self):
        if self.cb_static.isChecked(): 
            self.gg_dict['icollapse']=0
        if self.cb_collap.isChecked(): 
            self.gg_dict['icollapse']=1
        if self.cb_ss.isChecked(): 
            self.gg_dict['iSS']=1
        else:
            self.gg_dict['iSS']=0
        if self.cb_rd_g07.isChecked(): 
            self.gg_dict['iNTD']=1
        if self.cb_rd_m16.isChecked(): 
            self.gg_dict['iNTD']=2
        if not self.cb_rd_g07.isChecked() and not self.cb_rd_m16.isChecked():
            self.gg_dict['iNTD']=0
            
        self.gg_dict['fnet']=self.finputnet.text()
        self.gg_dict['fiab']=self.finputiab.text()
        self.gg_dict['fed'] =self.finputed.text()    
        self.gg_dict['nH'] =self.gas_nH_value.text()  
        self.gg_dict['Tgas'] =self.gas_T_value.text() 
        self.gg_dict['Av'] =self.gas_Av_value.text() 
        self.gg_dict['Zeta'] =self.gas_Zeta_value.text() 
        self.gg_dict['Chi'] =self.gas_Chi_value.text() 
        self.gg_dict['d2gmr'] =self.d2gmr_value.text() 
        self.gg_dict['Tdust'] =self.dust_T_value.text() 
        self.gg_dict['rdust'] =self.dust_R_value.text() 
        self.gg_dict['rhodust'] =self.dust_rho_value.text() 
        self.gg_dict['Rdb'] =self.dust_Rdb_value.text() 
        self.gg_dict['atol'] =self.ode_atol_value.text() 
        self.gg_dict['rtol'] =self.ode_rtol_value.text() 
        self.gg_dict['modelname'] =self.output_modelname_value.text() 
        self.gg_dict['finalage'] =self.output_finalage_value.text() 
        self.gg_dict['timestep'] =self.output_timestep_value.text() 
    
    def set_info(self):
        if int(self.gg_dict['icollapse'])==0: 
            self.cb_static.setChecked(True)
            self.cb_collap.setChecked(False)
        elif int(self.gg_dict['icollapse'])==1: 
            self.cb_static.setChecked(False)
            self.cb_collap.setChecked(True)   
             
        if int(self.gg_dict['iSS'])==1: 
            self.cb_ss.setChecked(True)   
        elif int(self.gg_dict['iSS'])==0:
            self.cb_ss.setChecked(False) 
            
        if int(self.gg_dict['iNTD'])==1: 
            self.cb_rd_g07.setChecked(True) 
            self.cb_rd_m16.setChecked(False) 
        if int(self.gg_dict['iNTD'])==2: 
            self.cb_rd_g07.setChecked(False) 
            self.cb_rd_m16.setChecked(True) 
        if int(self.gg_dict['iNTD'])==0: 
            self.cb_rd_g07.setChecked(False) 
            self.cb_rd_m16.setChecked(False) 
            
        self.finputnet.setText( self.gg_dict['fnet']) 
        self.finputiab.setText( self.gg_dict['fiab']) 
        self.finputed.setText( self.gg_dict['fed'] )    
        self.gas_nH_value.setText( self.gg_dict['nH'])
        self.gas_T_value.setText( self.gg_dict['Tgas'])  
        self.gas_Av_value.setText( self.gg_dict['Av']  ) 
        self.gas_Zeta_value.setText( self.gg_dict['Zeta']) 
        self.gas_Chi_value.setText( self.gg_dict['Chi'] ) 
        self.d2gmr_value.setText( self.gg_dict['d2gmr']) 
        self.dust_T_value.setText( self.gg_dict['Tdust'])
        self.dust_R_value.setText( self.gg_dict['rdust'])
        self.dust_rho_value.setText( self.gg_dict['rhodust']) 
        self.dust_Rdb_value.setText( self.gg_dict['Rdb']    ) 
        self.ode_atol_value.setText( self.gg_dict['atol']   ) 
        self.ode_rtol_value.setText( self.gg_dict['rtol']   ) 
        self.output_modelname_value.setText( self.gg_dict['modelname'])
        self.output_finalage_value.setText( self.gg_dict['finalage'] ) 
        self.output_timestep_value.setText( self.gg_dict['timestep'] ) 
        QApplication.processEvents()
        
    def load_pars_from_file(self):
        #fname = QFileDialog.getOpenFileName(self, 'Open file', '') 
        #print(fname)
        #self.load_value.setText(fname)
        try:
            x1,x2=np.loadtxt(self.load_value.text(), usecols=[0,1], unpack=True, dtype='str')
            for ik in range(x1.size):
                if x1[ik].strip() in self.gg_keys:
                    self.gg_dict[ x1[ik].strip() ]=x2[ik]
            #f_in.close()
            self.set_info()
            msg_f_in_ok = QMessageBox()
            msg_f_in_ok.setIcon(QMessageBox.Information)
            msg_f_in_ok.setText("Parameters are updated according to '%s'."%(self.load_value.text()) )
            msg_f_in_ok.exec_()
        except:
            msg_f_in = QMessageBox()
            msg_f_in.setIcon(QMessageBox.Information)
            msg_f_in.setText("'%s' is not found."%(self.load_value.text()) )
            msg_f_in.exec_()
        
                
    def save_pars_into_file(self):
        #self.filename = QFileDialog.getOpenFileName(self)
        f_out = open(self.save_value.text(),'w')
        
        self.collect_info()
        
        for ik in range(self.n_gg_keys):
            print('%s  %s'%( self.gg_keys[ik], self.gg_dict[self.gg_keys[ik]]), file=f_out)
        
        f_out.close()
        
        msg_f_out = QMessageBox()
        msg_f_out.setIcon(QMessageBox.Information)
        msg_f_out.setText("Parameters are saved into '%s'."%(self.save_value.text()) )
        msg_f_out.exec_()
    def update_load_save_content(self):
        self.load_text = 'models/'+ self.output_modelname_value.text().strip() +'.txt'
        self.load_tip = 'A file contains previous parameters:\n'+ self.load_text
        self.load_value.setText( self.load_text )
        self.load_value.setToolTip( self.load_tip )
        self.save_text = 'models/'+ self.output_modelname_value.text().strip() +'.txt'
        self.save_tip = 'A file contains previous parameters:\n'+ self.load_text
        self.save_value.setText( self.load_text )
        self.save_value.setToolTip( self.load_tip )
                    
    def change_to_static_model(self, state):
        if state == Qt.Checked:
            #self.setWindowTitle('Static model is checked.')
            self.cb_collap.setChecked(False)
            self.gas_nH_label.setText('n<sub>H</sub>=')
            self.gas_nH_value.setText('2e+4')
            self.gas_Av_label.setText('A<sub>V</sub>=')
            self.gas_Av_value.setText('10.0')
            self.output_timestep_value.setText('100')
            self.output_finalage_value.setText('1.0e7')
            self.output_modelname_value.setText('DC')
            self.update_load_save_content()

    
    def change_to_collap_model(self, state):
        if state == Qt.Checked:
            #self.setWindowTitle('Collapse model is checked.')
            self.cb_static.setChecked(False)
            self.gas_nH_label.setText('n<sub>Hi</sub>,n<sub>Hf</sub>=')
            self.gas_nH_value.setText('3000,1e+7')
            self.gas_Av_label.setText('A<sub>V0</sub>=')
            self.gas_Av_value.setText('2.0')
            self.output_timestep_value.setText('2000')
            self.output_finalage_value.setText('1.0e+7')
            self.output_modelname_value.setText('DC_collapse')
            self.update_load_save_content()
   
    def ModelnameOnChanged(self):
        self.update_load_save_content()

                    
    def change_to_g07(self, state):
        if state == Qt.Checked:
            #self.setWindowTitle('G07 is checked.')
            self.cb_rd_m16.setChecked(False)
        #else:
        #    self.setWindowTitle(' ')  
    
    def change_to_m16(self, state):
        if state == Qt.Checked:
            #self.setWindowTitle('M16 is checked.')
            self.cb_rd_g07.setChecked(False)
        #else:
        #    self.setWindowTitle('GGCHEM   ') 
    
    """        
    def on_nH_change(self):
        if self.cb_collap.isChecked()==1:
            if ',' in self.gas_nH_value.text():
                pass
            else:
                msg_err = QMessageBox()
                msg_err.setIcon(QMessageBox.Information)
                msg_err.setText("Error: For collapse, a comma is needed: e.g. '3000,1e+7'.")
                msg_err.exec_()
        QApplication.processEvents()    
    """            
                
    def changeTitle(self, state):

        if state == Qt.Checked:
            self.setWindowTitle('Static model is checked.')
        else:
            self.setWindowTitle(' ')        
                      
    #def keyPressEvent(self, e):
    #    """press 'ESC' to exit"""
    #    if e.key() == Qt.Key_Escape:
    #        self.close()
            
    def buttonClicked(self):
        sender = self.sender()
        self.statusBar().showMessage(sender.text() + ' was pressed')
        
    #def mousePressEvent(self, event):
    #    self.c.closeApp.emit()
    
    def change_output_info_to_running(self):
        self.output_info.setText('Running. Please wait.')
        QApplication.processEvents()

    """
    def paintEvent(self,evt):
        painter = QPainter(self)
        painter.setRenderHint(QPainter.Antialiasing, True)
        
        painter.setPen(Qt.back)
        painter.drawRect(self.rect())
        
        realSize = min(self.widt(), self.height()) - 10
    """    

    def runModel(self):
        self.pbar.setValue(0)
        
        ggchempy.iswitch.iSS=0
        ggchempy.iswitch.iNTD=0
        ggchempy.iswitch.icollapse=0
        
        if self.cb_ss.isChecked(): 
            ggchempy.iswitch.iSS=1
        if self.cb_rd_g07.isChecked(): 
            ggchempy.iswitch.iNTD=1
        if self.cb_rd_m16.isChecked(): 
            ggchempy.iswitch.iNTD=2
        if self.cb_static.isChecked(): 
            ggchempy.iswitch.icollapse=0
        if self.cb_collap.isChecked(): 
            ggchempy.iswitch.icollapse=1 
            
            
        input_network = self.finputnet.text()
        input_iabun   = self.finputiab.text()
        input_ed      = self.finputed.text()
        ggchempy.ggpars.ggfiles  = [input_network, input_ed, input_iabun] 
        
        
        ggchempy.gas.T     = float( self.gas_T_value.text() )
        ggchempy.gas.Zeta  = float( self.gas_Zeta_value.text() )
        ggchempy.gas.Chi   = float( self.gas_Chi_value.text() )
        ggchempy.gas.Av    = float( self.gas_Av_value.text() )
        
        if ggchempy.iswitch.icollapse==0:
            ggchempy.gas.nH    = float( self.gas_nH_value.text() )
        elif iswitch.icollapse==1:
            nHxx     = self.gas_nH_value.text().split(',')
            ggchempy.gas.nH0 = float(nHxx[0])
            ggchempy.gas.nH1 = float(nHxx[1])
            ggchempy.gas.Av0 = ggcehmpy.gas.Av
        
        ggchempy.ggpars.d2gmr = float(self.d2gmr_value.text())
        
        ggchempy.dust.T    = float( self.dust_T_value.text() )
        ggchempy.dust.radius = float( self.dust_R_value.text() )
        ggchempy.dust.rho    = float( self.dust_rho_value.text() )
        ggchempy.dust.surface.Rdb  = float( self.dust_Rdb_value.text() )
        
        ggchempy.ggpars.rtol = float(self.ode_rtol_value.text())
        ggchempy.ggpars.atol = float(self.ode_atol_value.text())
        
        modelname = self.output_modelname_value.text()
        ggchempy.ggpars.tf = float(self.output_finalage_value.text())
        ggchempy.ggpars.nt = int(self.output_timestep_value.text())
        
        """
        ##show a summary:
        if iswitch.icollapse==0: 
            coll_msg_nH = "nH = %10.3g"%(gas.nH)
        if iswitch.icollapse==1: 
            coll_msg_nH = 'nH = %10.3g -> %10.3g; Av0=%10.3g'%(gas.nH0,gas.nH1,gas.Av0)
        summary = 'Name=%s \niSS=%d, \niNTD=%d, \nicollapse=%d \n%s'%(modelname,\
                   iswitch.iSS, iswitch.iNTD, iswitch.icollapse, coll_msg_nH)
                   
        msg = QMessageBox()
        msg.setIcon(QMessageBox.Information)
        msg.setText(summary )
        msg.exec_()
        """
        
        self.change_output_info_to_running()
        
        #start_time = time.time()
        #ggchempy.run(modelname)
        #end_time = time.time()
        
        ggchempy.modelname = modelname
        ggchempy.init_ggchem()
        ggchempy.compute_reaction_rate_coefficients()
        
        ns = ggchempy.ggpars.ns
        nr = ggchempy.ggpars.nr
        rc = ggchempy.reactions.rc
        idx = ggchempy.reactions.idx
        
        @numba.jit(nopython=True )   
        def fode(t,y):
            """Function for scipy.integrate.ode:"""
            ydot = np.zeros(ns)
            up, down = np.zeros(ns),np.zeros(ns)
            for i in range(0,nr,1):
                IR1, IR2, IR3, IP1, IP2, IP3, IP4, IP5 = idx[i,:]
                term = 0.0
                if (IR1!=-1 and IR2==-1 and IR3==-1): term = rc[i]*y[IR1]
                if (IR1!=-1 and IR2!=-1 and IR3==-1): term = rc[i]*y[IR1]*y[IR2]
                if (IR1!=-1 and IR2!=-1 and IR3!=-1): term = rc[i]*y[IR1]*y[IR2]*y[IR3]
                ###if(term<0.0): term=0.0
                
                if (IP1!=-1): up[IP1]=up[IP1]+term
                if (IP2!=-1): up[IP2]=up[IP2]+term
                if (IP3!=-1): up[IP3]=up[IP3]+term
                if (IP4!=-1): up[IP4]=up[IP4]+term
                if (IP5!=-1): up[IP5]=up[IP5]+term
                            
                if (IR1!=-1): down[IR1]=down[IR1]+term
                if (IR2!=-1): down[IR2]=down[IR2]+term
                if (IR3!=-1): down[IR3]=down[IR3]+term
            ydot = up-down
            return ydot
        
        @numba.jit(nopython=True)
        def fjac(t,y):
            """Jacobian Function for scipy.integrate.ode:"""
            pd = np.zeros((ns,ns))
            for i in range(0,nr,1):
                IR1, IR2, IR3, IP1, IP2, IP3, IP4, IP5 = idx[i,:]
                term,term1,term2= 0.0,0.0,0.0
                
                if (IR2==-1 and IR3==-1):
                    ### no second reantant:
                    term = rc[i]
                    pd[IR1][IR1] = pd[IR1][IR1] -term
                    pd[IP1][IR1] = pd[IP1][IR1] +term
                    if (IP2!=-1): pd[IP2][IR1] = pd[IP2][IR1] +term
                    if (IP3!=-1): pd[IP3][IR1] = pd[IP3][IR1] +term
                    if (IP4!=-1): pd[IP4][IR1] = pd[IP4][IR1] +term
                    if (IP5!=-1): pd[IP5][IR1] = pd[IP5][IR1] +term
                else:
                    ### with  second reantant:
                    term1 = rc[i]*y[IR1]
                    term2 = rc[i]*y[IR2]
                    pd[IR1][IR1] = pd[IR1][IR1] - term2
                    pd[IR1][IR2] = pd[IR1][IR2] - term1
                    pd[IR2][IR1] = pd[IR2][IR1] - term2
                    pd[IR2][IR2] = pd[IR2][IR2] - term1
                    pd[IP1][IR1] = pd[IP1][IR1] + term2
                    pd[IP1][IR2] = pd[IP1][IR2] + term1
                    if (IP2!=-1): pd[IP2][IR1] = pd[IP2][IR1] + term2
                    if (IP2!=-1): pd[IP2][IR2] = pd[IP2][IR2] + term1
                    if (IP3!=-1): pd[IP3][IR1] = pd[IP3][IR1] + term2
                    if (IP3!=-1): pd[IP3][IR2] = pd[IP3][IR2] + term1
                    if (IP4!=-1): pd[IP4][IR1] = pd[IP4][IR1] + term2
                    if (IP4!=-1): pd[IP4][IR2] = pd[IP4][IR2] + term1
                    if (IP5!=-1): pd[IP5][IR1] = pd[IP5][IR1] + term2
                    if (IP5!=-1): pd[IP5][IR2] = pd[IP5][IR2] + term1
            return pd
            
        #### speed up with numba.jit:
        def wrapper_to_fode(t, y):
            return fode(t, y) 
        def wrapper_to_fjac(t,y):
            return fjac(t,y)
        
        ### open output file:
        if not os.path.isdir(ggchempy.ggpars.outdir[0:-1]):
            os.makedirs(ggchempy.ggpars.outdir[0:-1])
    
        fileout=open(ggchempy.ggpars.outdir+ggchempy.modelname+'.dat','w')
        print('Model name:', ggchempy.modelname)
        print('----------------- ggchempy 2021 -----------------', file=fileout)
        ggchempy.gas.print_attributes(tofile=fileout)
        ggchempy.dust.print_attributes(tofile=fileout)
        ggchempy.ggpars.print_attributes(tofile=fileout)
        print('--------------- time evolution-------------------', file=fileout)
        
        fmt = '%10s '*(ggchempy.ggpars.ns+1)            
        print( fmt%('time',*ggchempy.species.spec), file=fileout )
        fmt = '%10.3e '*(ggchempy.ggpars.ns+1)
                
        ### time grid
        ggchempy.tyr = np.logspace( np.log10(ggchempy.ggpars.ti), np.log10(ggchempy.ggpars.tf), ggchempy.ggpars.nt)
        ggchempy.ts  = ggchempy.tyr * ggchempy.ggpars.yr_sec
        
        ### DO the integration:     
        t0 = 0.0
        y0 = ggchempy.species.numb
        r  = ode(wrapper_to_fode, wrapper_to_fjac).set_integrator(\
                 'vode', method='bdf', \
                 rtol=ggchempy.ggpars.rtol, atol=ggchempy.ggpars.atol, \
                 with_jacobian=True, nsteps=500000)#, ixpr=True)
        r.set_initial_value(y0, t0)
        
        start_time = time.time()
        
        ######## start ODE solver:
        numbt = np.zeros((ggchempy.ggpars.nt, ggchempy.ggpars.ns))
        it = 0
        while r.successful() and r.t < ggchempy.ggpars.tf * ggchempy.ggpars.yr_sec:
            ggchempy.iprogress = it
            
            ### to update rate coefficients:
            ggchempy.compute_reaction_rate_coefficients()
            
            r.integrate(ggchempy.ts[it])
            r.y[np.argwhere(r.y<0.0)]=0.0

            numbt[it][:] = r.y
            ggchempy.species.numb = r.y
            ggchempy.species.abun = r.y/ggchempy.gas.nH
            
            nH_save      = np.float(ggchempy.gas.nH)
            numb_save    = np.array(r.y)
            
            ########################
            print( fmt%(r.t/ggchempy.ggpars.yr_sec, *r.y), file=fileout)
            
            it+=1
            self.pbar.setValue(it)
            
            if it==ggchempy.ggpars.nt: break

        
        ggchempy.ggpars.cputime = time.time() - start_time
        print("Finished, CPU time = %s seconds ---" % (ggchempy.ggpars.cputime) )
        fileout.close()
        
        ggchempy.res = {}
        for i in range(0,ns,1):
            ggchempy.res[ggchempy.species.spec[i]]=numbt[:,i]   ###number density (cm^-3)

        ggchempy.res['time']= ggchempy.tyr
        ggchempy.res['nH']  = ggchempy.gas.nH
        ggchempy.res['Av']  = ggchempy.gas.Av

    
        self.output_info.setText('CPU time =%10.3g (s)'%(ggchempy.ggpars.cputime))
        
        self.setWindowTitle( 'GGCHEMPY version 1.0')
        
        ##self.output_info.setText('Running')
        msg = QMessageBox()
        msg.setIcon(QMessageBox.Information)
        msg.setText("%s model is finished. \nYou can find your model in 'out/%s'.dat. \nPress OK to continue."%(modelname,modelname))
        msg.exec_()
        
        
        ##msg.setInformativeText("This is additional information")
        ##msg.setWindowTitle("MessageBox demo")
        ##msg.setDetailedText("The details are as follows:")
        ##retval = msg.exec_()
        ##msg.show()

def run_GUI():
    app = QApplication(sys.argv)
    ex = GGCHEMGUI()
    sys.exit(app.exec_())        
         
if __name__ == '__main__':

    app = QApplication(sys.argv)
    ex = GGCHEMGUI()
    sys.exit(app.exec_())

