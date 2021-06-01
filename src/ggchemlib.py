"""
This is GGCHEM libiary for gas-grain chemical simulations of interstellar medium.
"""
import numpy as np        
from numba import jit
from scipy.integrate import ode
from scipy.interpolate import interp1d
import sys
import time
import progressbar
from src.ggfuncs import *
import os

try:
    from pyfiglet import Figlet
    ff = Figlet(font='slant')
    print(ff.renderText('ggchem v 1 . 0'))
except:
    pass
    
class gas:
    nH=2.0e+4   # density (cm^-3)
    T=10.0      # temperature (K)
    Av=10.0     # visual extinction (mag)
    Chi=1.0     # radiation field (chi0)
    Zeta=1.0e-17# cosmic ionization (s^-1)
    ### for collapsing core when iswitch.icollapse=1:
    nH0 = 3000.0 # initial density
    nH1 = 1e7    # final density
    Av0 = 2.0    # base visual extinction at initial density
    nHs = []     # array to save density at each step
    Avs = []     # array to save extinction at each step

class dust:
    nbins=0      # nbins>1 means dust size distribution mode. Then sample dust into nbins.
    T=10.0       # temperature (K)
    radius=1.0e-5# dust grain radius (cm)
    rho=3.0      # mass density (g cm^-3)
    Tpeak=70.0   # peak temperature heated by CR particle.
    nd = 0.0     # dust number density. Will be calculated by nH with d2gmr.
    rds = []    # to save radii when nbins>1
    nds = []    # to save nd when nbins>1
    stick0=1.0  # neutral species
    stickn=0.0  # negative species
    stickp=0.0  # positive species
    site_density=1.5e+15
    site_number=0.0
    fc=3.0e-19
    class surface:
        Rdb=0.4
    class mantle:
        Rdb=0.8
        
class ggpars:
    nt=100
    ti=1.0
    tf=1.0e+8
    nstr=8
    nrnp=8
    pi=np.pi
    kb=1.38054e-16
    mp=1.66054e-24
    G = 6.6720e-8  
    B = 0.7  ##<=1.0, braking effects of turbulence and magnetic
    ne=0
    ns=0
    nr=0
    d2gmr  = 0.01
    ggfiles= ['in/network2.txt','in/ed.txt','in/iabun.txt']
    yr_sec = 365.0*24.0*3600.0
    rtol=1.0e-5
    atol=1.0e-15
    
    aa=0.01  # when iswitch.iNTD==1
    EM=100.0 # when iswitch.iNTD==2
    
    cputime=0.0
    
    outdir='out/'

class elements:
    spec=''
    mass=0.0
    
class species:
    spec=''
    mass=0.0   ## mass of species
    freq=0.0   ## frequency (s^-1)
    Ebind=0.0  ## binding energy (K)
    Ediff=0.0  ## Diffusion energy (K)
    charge=0.0 ## charge number 
    atoms=0.0  ## atom number e.g. for H2: 1 0 0 0 0 0 0 0 0 0 0 = H He C ... CL
    Natoms=0.0 ## total number of atom in a species
    abun=0.0   ## abundance
    numb=0.0   ## number density
    idx_CO = -1## -1 mean no index.
    idx_H2 = -1
    idx = {}
class reactions:
    spec = ''
    a=0.0    ## alpha
    b=0.0    ## beta
    c=0.0    ## gamma
    itype=0  ## reaction type
    rc=0.0   ## rate coefficient
    idx=0    ## index of species in reactions

class iswitch:
    iANA=0   ## analysis mode
    iSS =0   ## Self-shielding of CO and H2
    iNTD=0   ## 0==OFF.  *** Non-Thermal Desorption, reactive desorption ***
             ## 1= Garrod et al., 2007;
             ## 2= Minissale et al., A&A 585, A24 (2016).
    icollapse=0
    
def init_ggchem():
    """
    *** PURPOSE:
        (1) initialize reaction network (including elements, species, reactions, 
        binding energies and initial abundances) according to given parameters.
        (2) The balance of each reaction will be checked.
        (3) The index of species will be initialized.
    """
    dust.nd = gas.nH*ggpars.mp*ggpars.d2gmr/(4.0/3.0*ggpars.pi*dust.radius**3*dust.rho)
    dust.site_number = 4.0* ggpars.pi * dust.radius**2 * dust.site_density
    
    #print(modeloutput)
    print('------------------------ggchem 2020-------------------------')
    print('gas.nH=',gas.nH)
    print('gas.T=',gas.T)
    print('gas.Av=',gas.Av)
    print('gas.Chi=',gas.Chi)
    print('gas.Zeta=',gas.Zeta)
    print('dust-to-gas-mass ratio:',ggpars.d2gmr)
    print('dust.nd=',dust.nd)
    print('dust.T=',dust.T)
    print('dust.surface.Rdb=',dust.surface.Rdb)
    print('ti=',ggpars.ti,'---> tf=',ggpars.tf,'(yr), nt=',ggpars.nt)
    
    network = open(ggpars.ggfiles[0],'r').readlines()
    fEbind = ggpars.ggfiles[1]
    fiabun = ggpars.ggfiles[2]
    
    ### get number of elements, species and reactions:
    L0 = network[0].split()
    ne, ns, nr = L0[0:3]
    ggpars.ne, ggpars.ns, ggpars.nr = int(ne), int(ns), int(nr)
    
    ### allocate array:
    elements.spec = np.zeros(ggpars.ne, dtype='U'+str(ggpars.nstr))
    elements.mass = np.zeros(ggpars.ne)
        
    species.spec = np.zeros(ggpars.ns, dtype='U'+str(ggpars.nstr))
    species.mass = np.zeros(ggpars.ns)
    species.freq = np.zeros(ggpars.ns)
    species.Ebind = np.zeros(ggpars.ns)
    species.Ediff = np.zeros(ggpars.ns)
    species.charge= np.zeros(ggpars.ns)
    species.atoms= np.zeros((ggpars.ns,ggpars.ne))
    species.Natoms= np.zeros(ggpars.ns)
    species.abun= np.zeros(ggpars.ns)
    species.numb= np.zeros(ggpars.ns)
    
    reactions.spec = np.zeros((ggpars.nr, ggpars.nrnp),dtype='U'+str(ggpars.nstr))
    reactions.a    = np.zeros(ggpars.nr)
    reactions.b    = np.zeros(ggpars.nr)
    reactions.c    = np.zeros(ggpars.nr)
    reactions.itype= np.zeros(ggpars.nr, dtype='int')
    reactions.rc   = np.zeros(ggpars.nr)
    reactions.idx  = np.zeros((ggpars.nr, ggpars.nrnp), dtype='int')
    
    ### get elements and mass of elements (amu):
    for i in range(1, ggpars.ne+1, 1):
        Lx = network[i].split()
        elements.spec[i-1]=Lx[0]
        elements.mass[i-1]=float(Lx[1])
        
    ### get species and compute species mass (amu):    
    ii=0
    for i in range(ggpars.ne+1+1, ggpars.ne+1+1+ggpars.ns, 1):
        Lx = network[i].split()
        spx = Lx[0].strip()
        
        spmassx = 0.0; Natomsx=0.0
        for j in range(0,ggpars.ne,1):
            species.atoms[ii,j]  = float(Lx[j+2])
            spmassx += elements.mass[j]*species.atoms[ii][j]
            Natomsx += species.atoms[ii][j]
        if spx=='E': spmassx=1.0/1832.0
        if spx[0:1]=='G': spmassx=4.0/3.0*ggpars.pi*dust.radius**3*dust.rho/ggpars.mp
        
        species.spec[ii]   = spx 
        if spx == 'CO': species.idx_CO = ii
        if spx == 'H2': species.idx_H2 = ii
        species.charge[ii] = float(Lx[1])
        species.mass[ii]   = spmassx
        species.Natoms[ii] = Natomsx
        species.idx[spx]=ii
        #if '+' in species.spec[ii] or '-' in species.spec[ii]:
        #    print(species.spec[ii],species.charge[ii])
        
        ii+=1
        
    ### get reactions and get indices:
    nrnp=ggpars.nrnp
    nstr=ggpars.nstr
    
    ii=0
    for i in range(ggpars.ne+ggpars.ns+1+1+1, ggpars.ne+1+1+1+ggpars.ns+ggpars.nr, 1):
        Lx = network[i]
        
        Latom = np.zeros(ggpars.ne,dtype='int')
        Ratom = np.zeros(ggpars.ne,dtype='int')
        Lcharge, Rcharge = 0, 0
        
        for j in range(0, nrnp, 1):
            reac = Lx[j*nstr : (j+1)*nstr].strip()
            reactions.spec[ii,j]=reac.strip()
            for k in range(0, ggpars.ns, 1):
                if reac.strip()==species.spec[k].strip():
                    #print(reac.strip(),species.spec[k])
                    reactions.idx[ii][j]=k
                    ### check balance:
                    if j<=2:
                        for m in range(0,ggpars.ne,1):
                            Latom[m]= Latom[m]+ species.atoms[k][m]
                            Lcharge = Lcharge + species.charge[k]
                    else:
                        for m in range(0,ggpars.ne,1):
                            Ratom[m]= Ratom[m]+ species.atoms[k][m]
                            Rcharge = Rcharge + species.charge[k]
            if reac=='': reactions.idx[ii][j]=-1
            
        for m in range(0,ggpars.ne,1):
            if Latom[m]!=Ratom[m]:
                print(ii,Latom)
                print(ii,Ratom)
                print('Un-balance of atoms of reaction',ii)
                sys.exit()
            if Lcharge!=Rcharge:
                print('Un-balance of electrons in raction',ii,':',Lcharge,Rcharge)
                sys.exit()
                        
        reactions.a[ii],reactions.b[ii],reactions.c[ii],reactions.itype[ii] =\
            float(Lx[64:73]), float(Lx[73:82]),float(Lx[82:91]), int(Lx[91:93])
        ii+=1
        
    ### get Binding energies:
    #Ebind, Ediff, Freq = np.zeros(ns), np.zeros(ns), np.zeros(ns)
    specx  = np.loadtxt(fEbind, usecols=[0], skiprows=2, dtype='str')
    Ebindx = np.loadtxt(fEbind, usecols=[1], skiprows=2)
    nx = np.size(specx)
    for i in range(0, ggpars.ns, 1):
        for j in range(0, nx, 1):
            if species.spec[i].strip()=='J'+specx[j].strip():
                species.Ebind[i]= Ebindx[j]
                species.Ediff[i]= Ebindx[j]*dust.surface.Rdb
                species.freq[i] = np.sqrt(2.0*dust.site_density*ggpars.kb*species.Ebind[i]/ggpars.pi**2/species.mass[i]/ggpars.mp)
            ## For gas-phase species. 
            ## For use in non-thermal desportion when iswitch.iNTD!=0 (iTYPE=14)
            if species.spec[i].strip()==specx[j].strip():
                species.Ebind[i]= Ebindx[j]    
                
    ### get initial abundance:
    L= open(fiabun,'r').readlines()
    nx = int(L[0].split()[0])
    specx  = np.loadtxt(fiabun, usecols=[0], skiprows=1, dtype='str')
    iabunx = np.loadtxt(fiabun, usecols=[1], skiprows=1)
    specx = specx[0:nx]
    iabunx = iabunx[0:nx]
    for i in range(0, ggpars.ns, 1):
        for j in range(0, nx, 1):
            if species.spec[i].strip()==specx[j].strip():   
                species.abun[i] = iabunx[j]
                species.numb[i] = iabunx[j]*gas.nH

        if species.spec[i].strip()=='G0' or species.spec[i].strip()=='GRAIN0': 
            species.abun[i] = dust.nd/gas.nH
            species.numb[i] = dust.nd
        #if species.abun[i]!=0.0: 
        #    print(species.spec[i],species.abun[i])
        

def compute_reaction_rate_coefficients():
    """
    *** PURPOSE:
        compute reaction rate coefficients of all the reactions.
    """
    #print('---compute_reaction_rate_coefficients(): gas.nH',gas.nH)
    for i in range(0, ggpars.nr, 1):
        ### determine dust grains:
        if dust.nbins<=1:
            rdx = dust.radius
            ndx = dust.nd
        else:
            idust = 0 ###???
            rdx = rds[idust]
            ndx = nds[idust]
            
        if reactions.itype[i]==0:
            ### ion -- dust
            reactions.rc[i] = reactions.a[i]*(gas.T/300.0)**reactions.b[i] * gas.nH/dust.nd
            
        elif reactions.itype[i]==1:
            ### cosmic ray ionization
            reactions.rc[i] = reactions.a[i]*gas.Zeta
            
        elif 2<=reactions.itype[i]<=12:
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
            reactions.rc[i] = reactions.a[i]*(gas.T/300.0)**reactions.b[i]*np.exp(-reactions.c[i]/gas.T)
            
        elif reactions.itype[i]==13:
            ### 13--- photo-ionization and photo-dissociation
            if iswitch.iSS==1:
                specx = reactions.spec[i,0].strip()
                if specx=='H2' or specx=='CO':
                    xh2 = species.abun[species.idx_H2]
                    xco = species.abun[species.idx_CO]
                    xk  = gas.Chi*reactions.a[i]*np.exp(-reactions.c[i]*gas.Av)
                    reactions.rc[i] = H2_CO_self_shielding(specx, gas.Av, xh2, xco, xk)
                    #print(specx,'self-shielding is considered: rc=',reactions.rc[i],xco,xh2)
                else:
                    reactions.rc[i] = gas.Chi*reactions.a[i]*np.exp(-reactions.c[i]*gas.Av)
            else:
                reactions.rc[i] = gas.Chi*reactions.a[i]*np.exp(-reactions.c[i]*gas.Av)
            
            
        elif reactions.itype[i]==14:
            ### 14--- surface reaction
            k1,k2,k3,k4=reactions.idx[i,0],reactions.idx[i,1],reactions.idx[i,2],reactions.idx[i,3]
            if (k1*k2==0): 
                print('## surface reaction (itype=14): no match!!!',k1,k2)
                sys.exit()
            pro=0.0
            if reactions.c[i]==0.0:
                if (k1==k2): pro=0.5
                if (k1!=k2): pro=1.0
            else:
                pro = np.exp(-reactions.c[i]/dust.T)
            rdiff1 = species.freq[k1]*np.exp(-species.Ediff[k1]/dust.T)/dust.site_number
            rdiff2 = species.freq[k2]*np.exp(-species.Ediff[k2]/dust.T)/dust.site_number
            
            reactions.rc[i] = reactions.a[i]*pro*(rdiff1+rdiff2)/ndx ##dust.nd
            
            ## for network do not consider reaction desorption
            if 'J' not in species.spec[k4]: 
                reactions.rc[i]=0.0
            
            ## when consider reactive desorption (b in KJ/mol):
            if iswitch!=0 and reactions.b[i]<0.0:
                ## JA + JB -> C
                fx=0.0  ## fraction
                if 'J' not in reactions.spec[i,3] and reactions.spec[i,4].strip()=="":
                    if iswitch.iNTD==1:
                        ## Garrod et al., 2007
                        
                        fx = reactive_desorption_2007(reactions.b[i], \
                                                        species.Natoms[k4],\
                                                        species.Ebind[k4],
                                                        aa=ggpars.aa)
                        #print(fx,i,reactions.spec[i,0:5],reactions.b[i],species.Natoms[k4],species.Ebind[k4])
                    elif iswitch.iNTD==2:
                        ## Minissale et al., A&A 585, A24 (2016)
                        enthalpy = np.abs(reactions.b[i] * 120.274) ## KJ/mol->K
                        fx = reactive_desorption_2016(enthalpy,\
                                                        species.Natoms[k4],\
                                                        species.mass[k4],\
                                                        species.Ebind[k4],ggpars.EM)
                    
                    ## times fraction to the rate coefficient from last reaction (JA + JB -> JC).
                    reactions.rc[i]   = fx*reactions.rc[i-1]       ## current reaction
                    reactions.rc[i-1] = (1.0-fx)*reactions.rc[i-1] ## last reaction
                    #print(fx,reactions.rc[i-1],reactions.rc[i])
            
        elif reactions.itype[i]==15:
            ### 15--- thermal evaporation
            reactions.rc[i] = reactions.a[i]*species.freq[reactions.idx[i,0]]*np.exp(-species.Ebind[reactions.idx[i,0]]/dust.T)
            
        elif reactions.itype[i]==16:
            ### 16--- cosmic ray induced general desorption
            reactions.rc[i] = dust.fc*reactions.a[i]*species.freq[reactions.idx[i,0]]*np.exp(-species.Ebind[reactions.idx[i,0]]/dust.Tpeak)
            
        elif 17<=reactions.itype[i]<=18: 
            ### 17--- photodissociation by cosmic rays
            ###      (JAB -> JA + JB corresponds to AB -> A + B)
            ### 18--- photodissociation by cosmic rays
            ###      (JAB -> JA + JB corresponds to AB -> AB+ + E)
            reactions.rc[i] = reactions.a[i]*gas.Zeta

        elif 19<=reactions.itype[i]<=20:
            ### 19--- photodissociation by background photons
            ###      (JAB -> JA + JB corresponds to AB -> A + B)
            ### 20--- photodissociation by background photons
            ###      (JAB -> JA + JB corresponds to AB -> AB+ + E)
            reactions.rc[i] = gas.Chi*reactions.a[i]*np.exp(-reactions.c[i]*gas.Av)
            
        elif reactions.itype[i]==99:
            ### 99--- adsorption on grain surface
            vth = np.sqrt(8.0*ggpars.kb*gas.T/ggpars.pi/species.mass[reactions.idx[i,0]]/ggpars.mp)
            reactions.rc[i] = dust.stick0*ggpars.pi*rdx**2*vth*ndx
        
        
class ggchem:
    """
    GGCHEM model for gas-grain chemical simulations of interstellar medium.
    Author: Jixing Ge
    Email: gejixing666@gmail.com
    
    Usage:
        
        useGUI=False ## if True, PyQt5 is needed.
        
        if useGUI:
            ### Use the GUI to set model:
            from src.ggchemGUI import *
            app = QApplication(sys.argv)
            ex = GGCHEMGUI()
            sys.exit(app.exec_())
        else:
            from ggchemlib import *
            gas.nH =2e4
            gas.T =10.0
            gas.Zeta = 1.3e-17
            gas.Chi = 1.0
            gas.Av = 10.0
            dust.T=10.0
            dust.radius = 1.0e-5
            ggpars.d2gmr = 0.01
            ggpars.ggfiles = ['network2.txt', 'ed.txt', 'iabun.txt']
            ggchem.run('DarkCloud')
    """
    
    version='GGCHEM V1.0 (2020)'
        
    def run(modeloutput):
        """
        *** PURPOSE:
            run GGCEHM by calling functions as the following order:
            init_ggchem()
            compute_reaction_rate_coefficients()
            doing the integration of ODEs.
        """

        
        ### If dust.nbins>1, extend network
        if dust.nbins>1:
            print('### The dust-size-distribution function is not finished yet.')
            sys.exit()
        """
            fextended = ggpars.ggfiles[0][0:-4]+'_'+str(dust.nbins)+'bins'+ggpars.ggfiles[0][-4:]
            extend_to_dust_bins(ggpars.ggfiles[0], fextended, \
                                nbins=dust.nbins, nstr_new=12)
            ggpars.ggfiles[0] = fextended
            ggpars.nstr = 12
        """
        
        ### initialize network and parameters:
        init_ggchem()
        
        ### 
        compute_reaction_rate_coefficients()
        
        ns = ggpars.ns
        nr = ggpars.nr
        rc = reactions.rc
        idx = reactions.idx
        
        
        #print(rc,ns)
        
        @jit(nopython=True )   
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
        
        @jit(nopython=True)
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
        
        if not os.path.isdir(ggpars.outdir):
            os.makedirs(ggpars.outdir)
    
        fileout=open(ggpars.outdir + modeloutput + '.dat','w')
        print('Model name:',modeloutput)
        print('-------------ggchem 2020----------------', file=fileout)
        print('gas.nH=',gas.nH, file=fileout)
        print('gas.T=',gas.T, file=fileout)
        print('gas.Av=',gas.Av, file=fileout)
        print('gas.Chi=',gas.Chi, file=fileout)
        print('gas.Zeta=',gas.Zeta, file=fileout)
        print('dust-to-gas-mass ratio:',ggpars.d2gmr, file=fileout)
        print('dust.nd=',dust.nd, file=fileout)
        print('dust.T=',dust.T, file=fileout)
        #print('ti=',ggpars.ti,'---> tf=',ggpars.tf,'(yr), nt=',ggpars.nt)
            
        if iswitch.icollapse==1:
            fmt = '%10s '*(ggpars.ns+1+2)            
            print( fmt%('time',*species.spec,'nH','Av'), file=fileout )
            fmt = '%10.3e '*(ggpars.ns+1+2)
        else:
            fmt = '%10s '*(ggpars.ns+1)            
            print( fmt%('time',*species.spec), file=fileout )
            fmt = '%10.3e '*(ggpars.ns+1)
        
        
        ######## initialize ODE solver:
        tyr  = np.logspace(np.log10(ggpars.ti), np.log10(ggpars.tf), ggpars.nt) 
        t  = tyr * ggpars.yr_sec
        
        ###########################################
        ###  then generate varying density with time
        def fden(t,y):
            pi= ggpars.pi
            mH= ggpars.mp*1.4
            G = ggpars.G   
            B = ggpars.B
            term1 = (y[0]**4/gas.nH0)**(1.0/3.0)
            term2 = 24.0*pi*mH*G*gas.nH0
            term3 = (y[0]/gas.nH0)**(1.0/3.0)-1.0
            return B*term1*(term2*term3)**(1.0/2.0)
        
        if iswitch.icollapse==1:  
            ### estimate free-fall collapse time scale and reset time grid:
            #tff  = np.sqrt(3.0*ggpars.pi/ggpars.G/(ggpars.mp*gas.nH1) ) / ggpars.yr_sec
            #tyr  = np.logspace(np.log10(ggpars.ti), np.log10(tff), ggpars.nt) 
            #t  = tyr * ggpars.yr_sec
              
            gas.nHs = np.zeros(ggpars.nt)
            gas.Avs = np.zeros(ggpars.nt)
            
            y0=[gas.nH0*1.3]
            t0=0.0
            r=ode(fden).set_integrator('vode', method='bdf', with_jacobian=False)
            r.set_initial_value(y0, t0)
            for it in range(0, ggpars.nt, 1):
                r.integrate(t[it])
                gas.nHs[it] = r.y
                gas.Avs[it] = gas.Av0 + (gas.nHs[it]/gas.nH0)**(2.0/3.0)  ###2018ApJ...854...13P
                if r.y>=gas.nH1: break
        
            told = t
            ### reset time grid according to the final density gas.nH1
            tff  = interp1d(gas.nHs,told/ggpars.yr_sec)(gas.nH1)
            tyr  = np.logspace(np.log10(ggpars.ti), np.log10(tff), ggpars.nt) 
            t    = tyr * ggpars.yr_sec
            
            ### reset gas.nHs and gas.Avs
            gas.nHs = interp1d(told, gas.nHs)(t)
            gas.Avs = interp1d(told, gas.Avs)(t)
            print('### Collapse mode: density from %10.3g to %10.3g around %10.3g yr'%(gas.nH0,gas.nH1,tff))
         
            
        ################################################################
        ### DO the integration:     
        t0 = 0.0
        y0 = species.numb
        r  = ode(wrapper_to_fode, wrapper_to_fjac).set_integrator(\
                 'vode', method='bdf', \
                 rtol=ggpars.rtol, atol=ggpars.atol, \
                 with_jacobian=True, nsteps=500000)#, ixpr=True)
        r.set_initial_value(y0, t0)
        
        start_time = time.time()
        
        ######## start ODE solver:
        numbt = np.zeros((ggpars.nt,ggpars.ns))
        
        widgets = [
        '\x1b[33mggchem:\x1b[39m',
        progressbar.Percentage(),
        progressbar.Bar(marker='\x1b[32m#\x1b[39m'),
        ]
        
        #loop=tqdm(total=ggpars.nt,position=0,leave=True)
        bar = progressbar.ProgressBar(widgets=widgets, maxval=ggpars.nt).start()
        
        it = 0
        while r.successful() and r.t < ggpars.tf*ggpars.yr_sec:
            ### to update rate coefficients:
            compute_reaction_rate_coefficients()
            ########################
            if iswitch.icollapse==1:
                gas.nH = gas.nHs[it]
                gas.Av = gas.Avs[it] 
                
                ### to ensure that abundances from last step as 
                ### the input abundance for current step:
                if it>0:
                    for i in range(0,ggpars.ns,1):
                        r.y[i] = numb_save[i]/nH_save * gas.nH
            ########################
            
            r.integrate(t[it])
            r.y[np.argwhere(r.y<0.0)]=0.0

            numbt[it][:] = r.y
            species.numb = r.y
            species.abun = r.y/gas.nH
            
            nH_save      = np.float(gas.nH)
            numb_save    = np.array(r.y)
            
            ########################
            if iswitch.icollapse==1:
                print( fmt%(r.t/ggpars.yr_sec, *r.y, gas.nHs[it],gas.Avs[it]), file=fileout)
            else:
                print( fmt%(r.t/ggpars.yr_sec, *r.y), file=fileout)
                
            it+=1
            bar.update(it)
            #loop.set_description("ggchem:".format(it))
            #loop.update(1)
            if it==ggpars.nt: break
        bar.finish()
        #loop.close()
        
        ggpars.cputime = time.time() - start_time
        print("Finished, CPU time = %s seconds ---" % (ggpars.cputime) )
        fileout.close()
        
        res = {}
        for i in range(0,ns,1):
            res[species.spec[i]]=numbt[:,i]   ###number density (cm^-3)

        res['time']= tyr
        if iswitch.icollapse==1:
            res['nH']  = gas.nHs
            res['Av']  = gas.Avs
        else:
            res['nH']  = gas.nH
            res['Av']  = gas.Av
        #self.results = res
        return res
    
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
    
def ggreport(age, model, idx_sp, spec, reac, idx, itype, rc, rct, \
           nrc, ireac_con, nrp, ireac_pro, \
           nr, nrnp, nt, nfirst, fout):
    """
    *** PURPOSE:
        Report the most important production and consumption reactions of a given species.
    *** PARAMETERS:
        age    --- the age at which the important reactions of species will be found.
        model  --- the GGCHEM model file loaded by ggchemlib.loadgg().
        idx_sp --- the index of the species will be analyzed.
        spec   --- the list of all species ( spec=species.spec ).
        reac   --- the list of all reactions ( reac = reactions.spec ).
        rc     --- the array of reaction rate ceofficients (rc = reactions.rc).
        rct    --- the arrary of time-dependent reaction rate coefficient (rct(nt,nr)).
        nrc       --- number of consumption reaction of the species with idx_sp.
        ireac_con --- index of consumption reaction of the species with idx_sp.
        nrp       --- number of production reaction of the species with idx_sp.
        ireac_pro --- index of production reaction of the species with idx_sp.
        nr     --- number of total reactions in the network (nr=ggpars.nr).
        nrnp   --- number of species in a reaction (nrnp=ggpars.nrnp).
        nt     --- number of time steps used in the model.
        nfirst --- number of the most important reactions to be reported.
                   if nfirst>nrc or nfirst>nrp, nrc/nrp will be used as nfrist.
        fout   --- the file name to store the output reactions.
    """
    if iswitch.icollapse==1: print('--- Collapse mode ---')
    if iswitch.icollapse==0: print('--- Static mode ---')
    
    ### rates for all reactions at all time step (nt):
    rate_con = np.zeros((nrc,nt))
    rate_pro = np.zeros((nrp,nt))
    
    ### total rates:
    total_rate_con = np.zeros(nt)
    total_rate_pro = np.zeros(nt)
    
    ### interpolate rate at given age:  
    rate_con_age = np.zeros(nrc)
    rate_pro_age = np.zeros(nrp)
    
    ### Loop all times:
    for it in range(nt):
        ### Calculate rates for all consumption reactions:
        icx=0
        for ic in ireac_con:
            idx1  = -1;  idx2 = -1
            idx1  = idx[ic,0]
            idx2  = idx[ic,1]
            numb1 = 0.0;  numb2 = 0.0
            if idx1!=-1:
                sp1   = spec[idx1]
                numb1 = model[sp1][it]
                
            if idx2!=-1:
                sp2   = spec[idx2]
                numb2 = model[sp2][it]
                
            if iswitch.icollapse==1:
                if it>0:
                    numb1 = numb1/model['nH'][it-1] * model['nH'][it]
                    if idx2!=-1:
                        numb2 = numb2/model['nH'][it-1] * model['nH'][it]
                if idx1!=-1 and idx2==-1: rate_con[icx,it]= rct[it,ic] * numb1
                if idx1!=-1 and idx2!=-1: rate_con[icx,it]= rct[it,ic] * numb1 * numb2
            
            if iswitch.icollapse==0 :   
                if idx1!=-1 and idx2==-1: rate_con[icx,it]= rc[ic] * numb1
                if idx1!=-1 and idx2!=-1: rate_con[icx,it]= rc[ic] * numb1 * numb2
            
            total_rate_con[it] += rate_con[icx,it]
            icx += 1
        
           
        ### Calculate all rates for production reactions:    
        ipx=0
        for ip in ireac_pro:
            idx1  = -1; idx2 = -1
            idx1  = idx[ip,0]
            idx2  = idx[ip,1]
            
            numb1 = 0.0;  numb2 = 0.0
            if idx1!=-1:
                sp1   = spec[idx1]
                numb1 = model[sp1][it]
                
            if idx2!=-1:
                sp2   = spec[idx2]
                numb2 = model[sp2][it]
                
            if iswitch.icollapse==1:
                if it>0:
                    numb1 = numb1/model['nH'][it-1] * model['nH'][it]
                    if idx2!=-1:
                        numb2 = numb2/model['nH'][it-1] * model['nH'][it]
                if idx1!=-1 and idx2==-1: rate_pro[ipx,it]= rct[it,ip] * numb1
                if idx1!=-1 and idx2!=-1: rate_pro[ipx,it]= rct[it,ip] * numb1 * numb2
            
            if iswitch.icollapse==0:        
                if idx1!=-1 and idx2==-1: rate_pro[ipx,it]= rc[ip] * numb1
                if idx1!=-1 and idx2!=-1: rate_pro[ipx,it]= rc[ip] * numb1 * numb2
                
            total_rate_pro[it] += rate_pro[ipx,it]    
            ipx += 1
    
    ### interpolate at given age:            
    icx=0
    for ic in ireac_con:
        rate_con_age[icx] = interp1d(model['time'], rate_con[icx,:])(age)            
        icx+=1
    ipx=0
    for ip in ireac_pro:
        rate_pro_age[ipx] = interp1d(model['time'], rate_pro[ipx,:])(age)
        ipx+=1
                   
    ### report the most important reactions (numbe = nfirst):
    fmt = '%5d    %3d ' + '%10.3g  %10.3g%% ' + 3*' %8s' + '->' + 5*' %8s'
    sum_con_age = np.sum(rate_con_age)
    sum_pro_age = np.sum(rate_pro_age)
    print('Total production rate  = %10.3e at age of %10.3g yr. '% (sum_pro_age,age))
    print('Total consumption rate = %10.3e at age of %10.3g yr. '% (sum_con_age,age))
    
    idx_pro=[]; per_pro=[]
    iorder_pro = np.argsort(rate_pro_age)
    nfirst_pro = nfirst
    if nfirst>nrp: nfirst_pro=nrp
    print('nfirst:',nfirst_pro, file=fout)
    print('Total production rate  = %10.3e at age of %10.3g yr. '% (sum_pro_age,age), file=fout)
    print("  idx    itype   rate         percent ======= reactants ========->========== prodcuts ======== ",file=fout)
    for i in range(nrp-1, nrp-nfirst_pro-1, -1):
        idx_i  = ireac_pro[iorder_pro][i]
        rate_i = rate_pro_age[iorder_pro][i]
        reac_i = reac[ireac_pro[iorder_pro][i],:]
        type_i = itype[ireac_pro[iorder_pro][i]]
        percent = rate_i/sum_pro_age*100.0
        print( fmt % (idx_i, type_i, rate_i, percent, *reac_i[0:3], *reac_i[3:nrnp] ) , file=fout)
        idx_pro.append(idx_i)
        per_pro.append(percent)
    print('----------------------------------------------------------------------------------------------', file=fout)
    
    idx_con=[]; per_con=[]
    iorder_con = np.argsort(rate_con_age)
    nfirst_con = nfirst
    if nfirst>nrc: nfirst_con=nrc
    print('nfirst:',nfirst_con, file=fout)
    print('Total consumption rate = %10.3e at age of %10.3g yr. '% (sum_con_age,age), file=fout)
    print("  idx    itype   rate         percent ======= reactants ========->========== prodcuts ======== ",file=fout)
    for i in range(nrc-1, nrc-nfirst_con-1, -1):
        idx_i  = ireac_con[iorder_con][i]
        rate_i = rate_con_age[iorder_con][i]
        reac_i = reac[ireac_con[iorder_con][i],:]
        type_i = itype[ireac_con[iorder_con][i]]
        percent = rate_i/sum_con_age*100.0
        print( fmt % (idx_i, type_i, rate_i, percent, *reac_i[0:3], *reac_i[3:nrnp] ) , file=fout)
        idx_con.append(idx_i)
        per_con.append(percent)
    print('----------------------------------------------------------------------------------------------', file=fout)
    
    return total_rate_pro, total_rate_con, sum_pro_age, sum_con_age, idx_pro, idx_con, per_pro, per_con

"""    
def wrapper_to_report(age, model, idx_sp, spec, reac, idx, itype, rc, \
                      nrc, reac_con, nrp, reac_pro, \
                      nr, nrnp, nt, nfirst, fout):
    return report(age, model, idx_sp, spec, reac, idx, itype, rc, nrc, \
                  reac_con, nrp, reac_pro, \
                  nr, nrnp, nt, nfirst, fout)
"""

def ggtransform(reacx,percent):
    """
    *** PURPOSE:
        transfer reaction to latex format with appended contribution percent.
    *** INPUT:
        reacx --- reaction e.g. ['A','B','','C','D','','','']
        percent --- percent to the total rate
    *** OUTPUT:
        reacy --- reaction in latex format e.g. r"(percent)$A$ + $B$ -> $C$ + $D$"
    """
    ### replace E with e-:
    for i in range(0,ggpars.nrnp,1):
        if reacx[i]=='E':
            reacx[i]='e-'
            
    reacy = r"(%4.1f%%) " % percent
    reacy = reacy + latex_species(reacx[0])
    for i in range(1,3,1):
        if reacx[i]!="":
            reacy = reacy + ' + ' + latex_species(reacx[i])
    reacy = reacy + ' -> '
    reacy = reacy + latex_species(reacx[3])
    for i in range(4,ggpars.nrnp,1):
        if reacx[i]!="":
            reacy = reacy + ' + ' + latex_species(reacx[i])
    return reacy 
                    
def analyze(model, spec='CO', age=1e5, nfirst=10, fout=None, plot_nfirst=0,fn='x.pdf'):
    """
    *** PURPOSE: 
        Analyze the reaction network for a given species to find the most important 
        reactions at a given age.
        
    *** PARAMETERS:
        model --- ggchem model loaded by ggchemlib.loadgg()
        spec  --- species [str]
        age   --- the age (yr) to do the analysis [float]
        nfirst --- The nfirst most important reactions will be reported. [int]
        plot_nfirst --- if not plot_nfrist!=0, then plot a rate figure 
                        to show the most important reactions (N=plot_nfirst) [int]
        fout  --- output file. e.g. fout=open('analysis.out','w') in your main code.
    
    *** EXAMPLE:
        from src.ggchemlib import *
        init_ggchem()
        compute_reaction_rate_coefficients()
        dc = loadgg('out/test.dat')
        fout = open('analysis.out','w')
        age=5e5
        analyze(dc, spec='CO', age=age, fout=fout)
        analyze(dc, spec='N2H+', age=age, fout=fout)
        fout.close()
    
    """
    
    specx  = species.spec
    reac   = reactions.spec
    rc     = reactions.rc
    
            
    idx    = reactions.idx
    itype  = reactions.itype
    
    nr    = ggpars.nr
    nrnp  = ggpars.nrnp
    nt    = ggpars.nt
    idx_sp = species.idx[spec]
    
    ### rate coefficients for the collapse mode:
    rct = np.zeros( (nt,nr) )  ### time-dependent rc
    if iswitch.icollapse==1:
        for it in range(nt):
            gas.nH = model['nH'][it]
            gas.Av = model['Av'][it]
            compute_reaction_rate_coefficients()
            rct[it,0:nr]=reactions.rc
    
    ### first, find and save reactions with index:
    ireac_con = []#np.array(1, dtype='int')
    ireac_pro = []#np.array(1, dtype='int')
    for i in range(0,nr,1):
        if idx_sp in idx[i,0:3]:
            ireac_con.append(i)
        elif idx_sp in idx[i,3:nrnp]:
            ireac_pro.append(i)
    ireac_con = np.array(ireac_con)
    ireac_pro = np.array(ireac_pro)
    nrc = ireac_con.size
    nrp = ireac_pro.size
    
    print(iswitch.icollapse,ggpars.nt)
    print('')
    print('Analyzing',spec,'...')
    print(spec, file=fout)
    print('Production reactions:',nrp, file=fout)
    print('Consumption reactions:',nrc, file=fout)
    #print('Time steps:',nt, file=fout)
    print('Age:',age, file=fout)
    
    total_rate_pro, total_rate_con, sum_pro, sum_con, idx_pro, idx_con, per_pro, per_con = \
                      ggreport(age, model, idx_sp, specx,\
                      reac, idx, itype, rc, rct, \
                      nrc, ireac_con, nrp, ireac_pro, \
                      nr, nrnp, nt, nfirst, fout)
                      
    ### plot rates:  
    max_rate = np.max([sum_pro, sum_con])                
    if plot_nfirst!=0:
        import pylab as plt
        
        plt.figure( figsize=(10,5) )
        plt.figtext(0.51, 0.95, latex_species(spec), fontsize=14)
        plt.subplot(1,2,1)
        plt.title('Production')
        plt.plot(model['time'], total_rate_pro,'-',lw=4,color='#A9A9A9',label='total')
        for i in range(plot_nfirst):
            idxi = idx_pro[i]
            idx1 = idx[idxi,0]
            idx2 = idx[idxi,1]
            sp1 = species.spec[idx1]
            sp2 = species.spec[idx2]
            if idx1!=-1 and idx2==-1: ratex = rc[idxi]  * model[sp1] 
            if idx1!=-1 and idx2!=-1: ratex = rc[idxi]  * model[sp1] * model[sp2]
            reacy = ggtransform(reactions.spec[idxi,0:nrnp], per_pro[i] )
            plt.plot(model['time'], ratex, '-', label=reacy)
        plt.plot(age, sum_pro, 'ks')
        plt.xscale('log')
        plt.yscale('log')
        plt.axvline(age, color='#A9A9A9')
        plt.axhline(sum_pro, color='#A9A9A9')
        plt.text(1.2e3,sum_pro*1.2,'total rate at %8.1e yr'%age)
        plt.xlim(1e3,1e7)
        plt.ylim(max_rate/1e7, max_rate*10)
        plt.legend(fontsize=9,loc='lower right')
        plt.ylabel(r'Reaction rate (cm$^{-3}$ s$^{-1}$)')
        plt.xlabel('Time (yr)')
        
        plt.subplot(1,2,2)
        plt.title('Consumption')
        plt.plot(model['time'], total_rate_con,'-',lw=4,color='#A9A9A9',label='total')
        for i in range(plot_nfirst):
            idxi = idx_con[i]
            idx1 = idx[idxi,0]
            idx2 = idx[idxi,1]
            sp1 = species.spec[idx1]
            sp2 = species.spec[idx2]
            if idx1!=-1 and idx2==-1: ratex = rc[idxi]  * model[sp1] 
            if idx1!=-1 and idx2!=-1: ratex = rc[idxi]  * model[sp1] * model[sp2]
            reacy = ggtransform(reactions.spec[idxi,0:nrnp], per_con[i] )
            plt.plot(model['time'], ratex, '-', label=reacy)
        plt.plot(age, sum_con, 'ks')
        plt.xscale('log')
        plt.yscale('log')
        plt.axvline(age, color='#A9A9A9')
        plt.axhline(sum_con, color='#A9A9A9')
        plt.text(1.2e3,sum_con*1.2,'total rate at %8.1e yr'%age)
        plt.xlim(1e3,1e7)
        plt.ylim(max_rate/1e7, max_rate*10)
        plt.legend(fontsize=9,loc='lower right')
        plt.ylabel(r'Reaction rate (cm$^{-3}$ s$^{-1}$)')
        plt.xlabel('Time (yr)')
        
        plt.subplots_adjust(wspace=0.3,left=0.1,bottom=0.15,right=0.97,top=0.9)
        plt.savefig(fn)
        plt.close('all')
        print('Figure of "%s" is saved.'%fn)
    
    print(spec,'is finished.')
    
def extend_to_dust_bins(fnetwork,fout,nbins=9,nstr_new=12):
    """
    *** PURPOSE:
        This function is used to extend the reaction network to 
        include dust size distributions.
        For example, to extend 
        JH + JH -> JH2 
        to (n+1) bins dust grains. It becomes
        JH_0 + JH_0 -> JH2_0
        ...
        JH_n + JH_n -> JH2_n
    
    *** INPUT:
        fnetwork --- network file name
        fout ------- output file name
        nbins ------ bin number
        nstr_new --- new string length in the extended reaction netwrok.
                     It is important for a large species. 
                     For exalpme, JCH3OCH3 will be extended to JCH3OCH3_0, ..., 
                     JCH3OCH3_n. Thus, nstr=8 is not enough.
    *** OUTPUT: 
        A new file with name of fout in folder 'in/'    
    """
    nrnp=8  
    nstr=8
    print('ggchem: ggpars.nstr is set to %d.' %(nstr_new))
    
    fx = open(fout,'w')
    network = open(fnetwork,'r').readlines()
    ### get number of elements, species and reactions:
    L0 = network[0].split()
    ne, ns, nr = L0[0:3]
    ne, ns, nr = int(ne), int(ns), int(nr)
    print(ne,ns,nr,file=fx)
    
    ### get elements and mass of elements (amu):
    elements,emass=[],[]
    for i in range(1, ne+1, 1):
        Lx = network[i].split()
        elements.append(Lx[0])
        emass.append(float(Lx[1]))
        fmt='%'+str(nstr_new)+'s  %s'
        print(fmt%( Lx[0].ljust(nstr_new), Lx[1]),file=fx)
    
    print('',file=fx)   ### print a blank line
         
    ### get species and compute species mass (amu):    
    species=[]
    charge=np.zeros(ns,dtype='int')
    spatom=np.zeros((ns,ne),dtype='int')
    ii=0
    ns_new=0
    for i in range(ne+1+1, ne+1+1+ns, 1):
        Lx = network[i].split()
        spx = Lx[0].strip()
        species.append(spx) 
        charge[ii]   = int(Lx[1])
        for j in range(0,ne,1):
            spatom[ii][j]= int(Lx[j+2])
        
        ### copy species to dust grain bins:
        fmt='%'+str(nstr_new)+'s '+(ne+1)*'%3d'
        if 'J' in spx or 'GRAIN' in spx:
            for ib in range(0,nbins,1):
                spy = spx+'_'+str(ib)
                print(fmt%( spy.ljust(nstr_new),charge[ii], *spatom[ii,:]),file=fx)
                ns_new += 1
        else:
            spy = spx
            print(fmt%( spy.ljust(nstr_new),charge[ii], *spatom[ii,:]),file=fx)
            ns_new+=1
        ii+=1
        
    ### get reactions and get indices:
    nr_new=0
    print('',file=fx)
    a,b,c,itype = np.zeros(nr),np.zeros(nr),np.zeros(nr),np.zeros(nr,dtype='int')
    ii=0
    for i in range(ne+ns+1+1+1, ne+1+1+1+ns+nr, 1):
        reaction=''
        Lx = network[i]
        ### reaction rate parameters: alpha, beta, gamma                
        a[ii],b[ii],c[ii],itype[ii] = \
            float(Lx[nstr*nrnp   : nstr*nrnp+9]),\
            float(Lx[nstr*nrnp+9 : nstr*nrnp+18]),\
            float(Lx[nstr*nrnp+18: nstr*nrnp+27]),\
            int(Lx[nstr*nrnp+27  : nstr*nrnp+29])
        
        fmtx='%'+str(nstr_new*nrnp)+'s%9.2e%9.2e%9.2e%2d'  
        fmt='%'+str(nstr_new)+'s'   
        
        if itype[ii] in [0,14,15,16,17,18,19,20,99]:
            ### for dust-related reactions:
            for ib in range(0,nbins,1):
                nr_new+=1
                reaction=''
                for j in range(0, nrnp, 1):
                    reac = Lx[j*nstr : (j+1)*nstr].strip()
                    if 'J' in reac or 'GRAIN' in reac:
                        reacx=reac.strip()+'_'+str(ib)
                        reaction += fmt%(reacx.ljust(nstr_new))
                    else:
                        reaction += fmt%(reac.ljust(nstr_new))
                print(fmtx%(reaction.ljust(nstr_new*nrnp),a[ii],b[ii],c[ii],\
                            itype[ii]), file=fx)
        else:
            nr_new+=1
            reaction=''
            for j in range(0, nrnp, 1):
                reac = Lx[j*nstr : (j+1)*nstr].strip()
                reaction += fmt%(reac.ljust(nstr_new))
            print(fmtx%(reaction.ljust(nstr_new*nrnp),a[ii],b[ii],c[ii],\
                        itype[ii]), file=fx)
            
        ii+=1
    print('ggchem: "%s" is extended to contains %d species and %d reactions.' %\
         ( fnetwork, ns_new, nr_new))
    print('ggchem: New network is saved into %s with %d lines' %(fout,\
           1+ne+1+ns_new+1+nr_new) )
    fx.close()
    
    ### update the frist line with new ns,nr:
    network=open(fout,'r').readlines()
    network[0]= '%d  %d  %d  %d'%(ne, ns_new, nr_new, nbins)
    fx=open(fout,'w')
    for i in range(0,1+ne+1+ns_new+1+nr_new,1):
        print(network[i].replace('\n',''),file=fx)
    fx.close()
    
    return 
