
# Introduction:

GGCHEMPY: Gas-Grain CHEMical code for interstellar medium in Python3.

Author: Jixing Ge

E-mail: gejixing666@gmail.com


# Required Python Verion: 
    > Python 3.0

# Required packages:
    numba
    numpy
    scipy
    matplotlib
    progressbar
    pyfiglet (optional)   -> to create ASCII art text of "GGCHEM"


# Usage:
    <1> install GGCHEM:
    python setup.py build
    python setup.py install
    <2> Prepare your input files. See (2).
    <3> Run models. 
    
    Example:
    ```python
    ### Use the GUI to set model:
    from src.ggchemGUI import *
    app = QApplication(sys.argv)
    ex = GGCHEMGUI()
    sys.exit(app.exec_())
    
    ### Python script
    from src.ggchemlib import *
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
        GGCHEM[model] = ggchem.run(model) 
    
    ```
    See more in examples.py. 
# Benchmark with the five models of Semenov et al., (2010):
<img src="benchmark.png" alt="benchmark" style="width:120px;"/>

***Solid line***: GGCHEMPY

***Points***: model of Semenov et al., (2010).

# References:
    <1> for basic rate equation method:
    Hasegawa T. I., Herbst E., Leung C. M., 1992, ApJS, 82, 167
    Semenov D., et al., 2010, A&A, 522, A42
    
    <2> for reactive desorption:
    Garrod R. T., Wakelam V., Herbst E., 2007, A&A, 467, 1103
    Minissale M., Dulieu F., Cazaux S., Hocuk S., 2016, A&A, 585, A24
    
    <3> for GGCHEM in Fortran:
    Ge J. X., He J. H., Yan H. R., 2016, MNRAS, 455, 3570
    Ge J. X., He J. H., Li A., 2016, MNRAS, 460, L50
    Ge J., Mardones D., Inostroza N., Peng Y., 2020, MNRAS, 497, 3306
    Ge J. X., et al., 2020, ApJ, 891, 36

    <4> for reaction network:
    Semenov D., et al., 2010, A&A, 522, A42
    [KIDA database:](http://kida.astrophy.u-bordeaux.fr/)
    [UDFA database:](http://udfa.ajmarkwick.net/)
