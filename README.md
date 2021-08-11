
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
    pyqt5                 -> for using GUI
    pyfiglet (optional)   -> to create ASCII art text of "GGCHEM"


# Usage:
    <1> install GGCHEM:
    python setup.py build
    python setup.py install
    <2> Prepare your models.
    <3> Run models. 
    
For example, please type the following command on your terminal:
```
Python benchmark.py
```
<img src="figures/code-and-gui.png" alt="code ang GUI" style="width:120px;"/>
  
# Benchmark with the five models of Semenov et al., (2010):
<img src="figures/benchmark.png" alt="benchmark" style="width:120px;"/>

***Solid line***: GGCHEMPY

***Point***: model of Semenov et al., (2010).

# CPU times (s):
<img src="figures/cputime.png" alt="cpu time" style="width:120px;"/>

# Note:
There is a file "in/network2.txt" with updated reactions of HNCO.
Reactive desorption is only available for this updated reaction network.

# For citing GGCHEMPY:
A paper has been submitted. A link is coming soon.

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
    
# Acknowledgements

This work is accomplished with the support from the Chinese Academy of Sciences
(CAS) through a Postdoctoral Fellowship administered by the CAS South America Center for Astronomy
(CASSACA) in Santiago, Chile. 

Jixing  thanks Dr. Jinhua He and Dr. Tie Liu for their constructive suggestions
on improving the code and the paper

# Links:
[KIDA database](http://kida.astrophy.u-bordeaux.fr/)
  
[UDFA database](http://udfa.ajmarkwick.net/)

[GitFront page](https://gitfront.io/r/Jixing-Ge/7b3c064da2aa89e7c1d34633a426340610338338/ggchempy/)
