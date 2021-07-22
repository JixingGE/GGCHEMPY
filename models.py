from ggchempylib import ggchempy, run_GUI, latex_species, loadgg

"""
gg = ggchempy

gg.gas.nH=2e4
gg.gas.Av=10.0
gg.gas.T=10.0
gg.gas.Zeta=1.3e-17
gg.gas.Chi = 1.0
gg.dust.T=10.0
gg.dust.surface.Rdb=0.77
gg.ggpars.tf=1e9
gg.ggpars.nt=100
TMC1 = gg.run('xxx')
"""
"""
import matplotlib.pyplot as plt
for sp in ['C+','C','CO','C2H2N']:
    plt.loglog(TMC1['time'], TMC1[sp]/2e4, label=sp)
plt.legend()
plt.show()
"""


run_GUI()
print(ggchempy.res['time'])




