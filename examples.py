from src.ggchemlib import *

### The following three switches (logical) are only used for this script:
useGUI=False       ## If True, PyQt5 is needed.
runBenchmark=False  ## If True, NAUTILUS result will be downloaded from the bechmarking webpage.
                   ## Benchmark files will be saved into "benchmark/"
runAnalysis=False  ## If True, analysis mode is activated. See the end of this script.

### Main code:
if useGUI:
    ### Use the GUI to set model:
    from src.ggchemGUI import *
    app = QApplication(sys.argv)
    ex = GGCHEMGUI()
    sys.exit(app.exec_())
else:
    
    #################################### the static models:
    GGCHEM={}
    modelnames = ['TMC1','HOTCORE','DISK1','DISK2','DISK3']

    ### define common parameters:
    ggpars.d2gmr = 0.01
    ggpars.ti = 1.0
    ggpars.tf = 1.0e+9
    ggpars.ggfiles= ['in/network2.txt','in/ed.txt','in/iabun.txt']
    gas.Zeta = 1.3e-17
    dust.surface.Rdb=0.77
    
    ### the five static models:
    pars={}        # nH,     Tgas, Av,   chi,  Tdust
    pars['TMC1']   =[2.00e4, 10.0, 10.0, 1.0,  10.0]
    pars['HOTCORE']=[2.00e7, 100.0, 10.0, 1.0, 100.0]
    pars['DISK1']  =[5.41e8, 11.4, 37.1, 428.3, 11.4]
    pars['DISK2']  =[2.59e7, 45.9, 1.94, 393.2, 45.9]
    pars['DISK3']  =[3.67e6, 55.2, 0.22, 353.5, 55.2]
    
    
    for model in modelnames:
        gas.nH, gas.T, gas.Av, gas.Chi, dust.T = pars[model]
        GGCHEM[model] = ggchem.run(model)  ## diectory "out" will be created to save model. e.g. "out/TMC1.dat"

    if runBenchmark:
        import shutil
        import glob
        import pylab as pl
        ################################### For benchmark with Semenov's models:
        dir_bench = 'benchmark'
        if not os.path.isdir(dir_bench):
            os.makedirs(dir_bench)
        
        ### Download file:
        res_NAU_zip = "Results_NAUTILUS.zip"
        zipfile = dir_bench + "/"+ res_NAU_zip 
        if not os.path.isfile(dir_bench + "/" + res_NAU_zip):   
            from urllib.request import urlretrieve
            URL_NAUTILUS = "https://www2.mpia-hd.mpg.de/homes/semenov/Chemistry_benchmark/"+res_NAU_zip
            print('Downloading', res_NAU_zip, "into", zipfile)
            urlretrieve (URL_NAUTILUS, zipfile)
        
        ### extract files:
        
        shutil.unpack_archive(zipfile, dir_bench+"/")
        
        allzip = os.listdir(dir_bench+"/")
        for xzip in allzip:
            if '.zip' in xzip:
                shutil.unpack_archive(dir_bench+"/"+xzip, dir_bench+"/")
        
        allzip = os.listdir(dir_bench+"/")
        for xzip in allzip:
            if res_NAU_zip not in xzip and '.zip' in xzip:
                os.remove(dir_bench+"/"+xzip)
        
        ### load NAUTILUS results:
        
        alldat = glob.glob(dir_bench+'/*.dat')
        NAUTILUS = {}
        for datfile in alldat:
           if "TMC1" in datfile: modelname='TMC1'
           if "HOT_CORE" in datfile: modelname='HOTCORE'
           if "DISK" in datfile:
               if '_1.dat' in datfile: modelname='DISK1'
               if '_2.dat' in datfile: modelname='DISK2'
               if '_3.dat' in datfile: modelname='DISK3'    
           NAUTILUS[modelname] = loadgg(datfile,NAUTILUS=True)
        
        ### plot:   
        
        spec = ['CO','JCO','H2O','JH2O','HCO+','C2S','HNCO','JHNCO']
        cols = ['k', 'r', 'g',     'b',  'c',    'y', 'm', '#828282']
        ip=0
        pl.rc('font', size=9)
        pl.figure(figsize=(5,4))
        pl.subplots_adjust(top=0.95,hspace=0.5,right=0.98)
        for md in modelnames:
            ip+=1
            ax=pl.subplot(2,3,ip)
            ax.set_title(md)
            j=0
            leg=[]
            for sp in spec:
                p1,=ax.plot(np.log10(NAUTILUS[md]['#TIME(yr)']), np.log10(NAUTILUS[md][sp.replace('J','g')]) ,'o',ms=5,mfc='None',color=cols[j])
                p2,=ax.plot(np.log10(GGCHEM[md]['time']), np.log10(GGCHEM[md][sp]/GGCHEM[md]['nH']) ,'-',color=cols[j],label=latex_species(sp))
                j+=1
                leg.append(p2)
            ax.set_xlim(2,7)
            ax.set_ylim(-22,-3)
            if ip in [2,3,5]: pl.yticks([-20,-15,-10,-5],[])
            if ip in [1,4]: pl.ylabel('log(x)')
            pl.xlabel('log(t)')
            pl.xticks([2,3,4,5,6,7])
            if ip==5:
                pl.legend(loc=(1.2,-0.1))  
        pl.savefig(dir_bench+'/benchmark.pdf')
        print("See ",dir_bench+'/benchmark.pdf')
   

        if runAnalysis:
            ### For analysis of a species at a given age in the TMC1 model:
            init_ggchem()
            compute_reaction_rate_coefficients()
            dc = loadgg('out/TMC1.dat')
            fout = open('out/TMC1_analysis.out','w')
            age=5e5
            analyze(dc, spec='CO', age=age, fout=fout)
            analyze(dc, spec='N2H+', age=age, fout=fout)
            fout.close()


