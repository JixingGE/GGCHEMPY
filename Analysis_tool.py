# -*- coding: utf-8 -*-
"""
This script provides a GUI using pyqt5 to analyze a given species at a age.
A subwindow will show a matplotlib plot of the species as follows:

------------------------------|
                              |
Abundance as function of time |
                              |
------------------------------|
             |                |
Formation    |  Destruction   |
rates        |  rates         |
             |                |
------------------------------|

Sorry for the delayed reply.
I do not frequenly check the jxg@uchile.cl
This is my gmail.
I would like to contribute a talk in July or Agust.
However, we are encouraged to public places 
due to the current rules from my current affilation (CASSACA).
Is it possible to give the talk online if I cann't get there at that time?

Best regards,
jixing
"""
from PyQt5 import QtCore, QtGui, QtWidgets
from PyQt5.Qt import Qt
from pandas import read_table
import matplotlib
matplotlib.use('Qt5Agg')
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg
from matplotlib.backends.backend_qt5agg import NavigationToolbar2QT as NavigationToolbar
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.figure import Figure
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
from ggchempylib import ggchempy
import os
#from pyvalem.formula import Formula
import networkx as nkx
import platform
import numba

default_Rdb = ggchempy.dust.surface.Rdb

OS = platform.system()
font_size=9
if OS=='Linux': 
    font_size=9
    plt.rc('font', family='Times New Roman', size=font_size)
if OS=='Darwin': 
    font_size=4
    plt.rc('font', size=font_size)

nclick=0
#ggchempy=ggchempy()
elements={'HE':'He','FE':'Fe','MG':'Mg','NA':'Na','SI':'Si','CL':'Cl'}
def fill_click(msg):
    return f'<span><a style="color: blue" href="{msg}">{msg}</a></span>'
def format_reaction(r8, outtype='text'):
    """
    Format a reaction (array) to a string with clicked species.
    """
    for i in range(0, r8.size, 1):
        for el in elements:
            if el in r8[i]:
                r8[i] = r8[i].replace(el, elements[el])
        if r8[i]=='E': r8[i]='e-'
    reac=""
    reac_click=""
    if outtype=='latex':
        ireac=0
        if r8[1]!="":
            reac = r'$'+Formula(r8[0]).latex+'$' + ' + ' + r'$'+Formula(r8[1]).latex + '$ -> '
        else:
            reac = r'$'+Formula(r8[0]).latex + '$ -> ' 
        for i in range(3,8,1):
            if r8[i]!="":
                reac+= r'$'+Formula(r8[i]).latex + '$ + '
        reac = reac[0:-3]
    
    if outtype=='text':
        ireac=0
        if r8[1]!="":
            reac = r8[0] + ' + ' + r8[1] + ' -> '
        else:
            reac = r8[0] + ' -> ' 
        
        if r8[1]!="":
            reac_click = fill_click(r8[0]) + ' + ' + fill_click(r8[1]) + ' -> '
        else:
            reac_click = fill_click(r8[0]) + ' -> ' 
                
        for i in range(3,8,1):
            if r8[i]!="":
                reac+= r8[i] + ' + '
                
        for i in range(3,8,1):
            if r8[i]!="":
                reac_click+= fill_click(r8[i]) + ' + '
                
        reac = reac[0:-3]
        reac_click = reac_click[0:-3]
        return reac, reac_click

@numba.jit(nopython=True)
def analyze_numba(age, ns, nr, nt, timeyr, isp, ireac_map, rc_t, sp_numb):
    """
    Called by analyze() by passing numpy.array().
    ns,nr,nt --- number of species, reactions, time steps
    isp      --- index of species to be analyze [int]
    isp_map  --- index map of all species       [array(dtype='int')]
    ireac_map--- index map of all reactions     [array(dtype='int')]
    rc_t     --- rc as function of t (time)     [array(dtyle='float')]
                 rc --- reaction rate coefficient
    sp_numb  --- number density of all species  [array(dtype='float')]
    """
    ratet = np.zeros(nt)
    idx_for = [] #np.empty(0, dtype='int')
    rates_for = np.zeros((nr, nt))
    rate_at_age_for = [] #np.array([])
    idx_des = [] #np.empty(0, dtype='int')
    rates_des = np.zeros((nr, nt))
    rate_at_age_des = [] #np.array([])  
    
    ### compute all rates at each time step  
    nfor, ndes=0,0      
    for ir in np.arange(nr):
        isp1 = ireac_map[ir,0]
        isp2 = ireac_map[ir,1]
        for it in np.arange(nt):
            if isp2!=-1:
                ratet[it] = rc_t[it,ir]* sp_numb[isp1,it] * sp_numb[isp2,it]
            else:
                ratet[it] = rc_t[it,ir]* sp_numb[isp1,it]
    
    #for ir in np.arange(0,nr,1):    
        ifound=0
        for irp in range(8):
            if irp<=2 and isp==ireac_map[ir,irp]: 
                ifound=1
            elif irp>2 and isp==ireac_map[ir,irp]: 
                ifound=2
            
        if ifound==1:
            ### des 
            rates_des[ndes:] = ratet
            idx_des.append(ir)
            rate_at_age_des.append(np.interp(age, timeyr, ratet))
            ndes+=1
        elif ifound==2:
            ### for
            rates_for[nfor,:] = ratet
            idx_for.append(ir)
            rate_at_age_for.append(np.interp(age, timeyr, ratet))
            nfor+=1
            
    rate_at_age_des = np.array(rate_at_age_des)
    rate_at_age_for = np.array(rate_at_age_for)    
    #print(rates_for, rates_des)
    return rate_at_age_for, rate_at_age_des, rates_for, rates_des, idx_for, idx_des
           
# To handle model file from ggchempy and to do analysis of a given species at a age.
class ggmodel(object):
    def __init__(self):
        self.version=0.1
        self.reactions=[]
        self.species=[]
        self.model={}
        self.sp='CO'
        self.age=1e5
        self.nH=1e4
        self.default_Rdb=0.5
    def load_a_model(self, mf, fnetwork, fed, Rdb):
        self.model = read_table(mf, skiprows=28, sep='\s+')
        ### Get parameters from the model file and setup ggchempy.ggpars, ggchempy.gas and ggchempy.dust
        df = open(mf, 'r')
        lines = df.readlines()
        df.close()
        
        for i in range(0,8,1):
            L= lines[i]
            if ' nH ' in L: ggchempy.gas.nH = float(L.split()[2]); self.nH= ggchempy.gas.nH 
            if ' Av ' in L: ggchempy.gas.Av = float(L.split()[2])
            if ' Chi ' in L: ggchempy.gas.Chi = float(L.split()[2])
            if ' Zeta ' in L: ggchempy.gas.Zeta = float(L.split()[2])
            if ' T ' in L: ggchempy.gas.T = float(L.split()[2])
        for i in range(8,28,1): 
            L=lines[i]   
            if ' radius ' in L: ggchempy.dust.radius = float(L.split()[2]); self.nH= ggchempy.gas.nH 
            if ' site_density ' in L: ggchempy.dust.site_density = float(L.split()[2])
            if ' rho ' in L: ggchempy.dust.rho = float(L.split()[2])
            if ' nd ' in L: ggchempy.dust.nd = float(L.split()[2])
            if ' T ' in L: ggchempy.dust.T = float(L.split()[2])
            if ' mass ' in L: ggchempy.dust.mass = float(L.split()[2])
            if ' stick0 ' in L: ggchempy.dust.stick0 = float(L.split()[2])
            if ' fc ' in L: ggchempy.dust.fc = float(L.split()[2])
            if ' Tpeak ' in L: ggchempy.dust.Tpeak = float(L.split()[2])
            
            if ' d2gmr ' in L: ggchempy.ggpars.d2gmr = float(L.split()[2])
            if ' ne ' in L: ggchempy.ggpars.ns = int(L.split()[2])
            if ' ns ' in L: ggchempy.ggpars.ne = int(L.split()[2])
            if ' nr ' in L: ggchempy.ggpars.nr = int(L.split()[2])
            
            if 'Infiles:' in L:
                ggchempy.ggpars.ggfiles= L.split()[1:]
                #print(L.split()[1:])
            
        #ggchempy.ggpars.ggfiles=[fnetwork, fed, 'in/iabun.txt']
        if Rdb!=0.0: 
            ggchempy.dust.surface.Rdb=Rdb
        else:
            ggchempy.dust.surface.Rdb=self.default_Rdb
        ggchempy.init_ggchem()
        ggchempy.compute_reaction_rate_coefficients()
            
    def analyze(self, sp, age, ntop=5, mode='numba'):
        self.sp = sp
        self.age= age
        print('Analyzing %s at age of %10.3e'%(self.sp, self.age) )
        
        ns = ggchempy.ggpars.ns
        nr = ggchempy.ggpars.nr
        nt = ggchempy.ggpars.nt
        sp_idx = ggchempy.species.idx[sp]
        
        reactions = ggchempy.reactions
        
        if mode=='static':
            #reac_idx= np.array([])
            idx_for = np.array([])
            rates_for = np.zeros((ggchempy.ggpars.nr, ggchempy.ggpars.nt))
            rate_at_age_for = np.array([])
            idx_des = np.array([])
            rates_des = np.zeros((ggchempy.ggpars.nr, ggchempy.ggpars.nt))
            rate_at_age_des = np.array([])
            
            nfor=0
            ndes=0
            for ir in range(ggchempy.ggpars.nr):
                ### formation:
                if sp_idx in reactions.idx[ir,3:]:
                    #reac_idx = np.append(reac_idx, ir)
                    rate=0.0
                    if reactions.idx[ir,1]==-1:
                        numb1 = self.model[reactions.spec[ir,0]]
                        rate = ggchempy.reactions.rc[ir] * numb1
                    else:
                        numb1 = self.model[reactions.spec[ir,0]]
                        numb2 = self.model[reactions.spec[ir,0]]
                        rate = ggchempy.reactions.rc[ir] * numb1 * numb2
                    idx_for = np.append(idx_for,ir)
                    rates_for[nfor,:]=rate
                    rate_at_age_for = np.append(rate_at_age_for, np.interp(age, self.model['time'], rate) )
                    nfor+=1
                ### destruction:
                if sp_idx in reactions.idx[ir,0:3]:
                    #print(reactions.spec[ir,0:3])
                    rate=0.0
                    if reactions.idx[ir,1]==-1:
                        numb1 = self.model[reactions.spec[ir,0]]
                        rate = ggchempy.reactions.rc[ir] * numb1
                    else:
                        numb1 = self.model[reactions.spec[ir,0]]
                        numb2 = self.model[reactions.spec[ir,0]]
                        rate = ggchempy.reactions.rc[ir] * numb1 * numb2
                    idx_des = np.append(idx_des, ir)
                    rates_des[ndes,:]=rate
                    rate_at_age_des = np.append(rate_at_age_des, np.interp(age, self.model['time'], rate) )
                    ndes+=1
            
        if mode=='numba':
            ireac_map = ggchempy.reactions.idx
            sp_numb = np.zeros((ns,nt)) 
            for i in range(0, ns, 1):
                sp_numb[i,:] = self.model[ggchempy.species.spec[i]]
            
            rc_t = np.zeros((nt,nr))    
            for it in range(nt):
                ## update parameters 
                """
                ggchempy.gas.nH = self.model['nH'][it]
                gcghempy.gas.Av = self.model['Av'][it]
                ggchempy.gas.T  = self.model['Tgas'][it]
                ggchempy.dust.T = self.model['Tdust'][it]
                ggchempy.dust.nd= self.model['GRAIN0'][it]
                """
                ggchempy.compute_reaction_rate_coefficients()
                rc_t[it,:] = ggchempy.reactions.rc
            
            timeyr = self.model['time'].to_numpy()
                    
            rate_at_age_for, rate_at_age_des, rates_for, rates_des, idx_for, idx_des, \
             = analyze_numba(age, ns, nr, nt, timeyr, sp_idx, ireac_map, rc_t, sp_numb)
             
            
        print( 'found ', np.size(rate_at_age_for), ' formation reactions.')
        print( 'found ', np.size(rate_at_age_des), ' destruction reactions.')
            
        percent_for_at_age = 100.0* rate_at_age_for/np.sum(rate_at_age_for)
        percent_des_at_age = 100.0* rate_at_age_des/np.sum(rate_at_age_des)
            
        rates_for = np.array(rates_for)
        rates_des = np.array(rates_des)
        ifor = np.argsort(rate_at_age_for)
        ides = np.argsort(rate_at_age_des)
        
        return ifor, ides, rates_for, rates_des, reactions, idx_for, idx_des, \
               percent_for_at_age,  percent_des_at_age
        
        
    
"""
class MplCanvas(FigureCanvasQTAgg):
    def __init__(self, parent=None, width=5, height=4, dpi=100):
        fig = Figure(figsize=(width, height), dpi=dpi)
        self.axes = fig.add_subplot(111)
        super(MplCanvas, self).__init__(fig)
"""

# subwindow to show a plot from matplotlib        
class subwindow(QtWidgets.QWidget):
    def createWindow(self,WindowWidth,WindowHeight):
        parent=None
        super(subwindow,self).__init__(parent)
        l = QtWidgets.QGridLayout(self)
        self.setWindowFlags(QtCore.Qt.WindowStaysOnTopHint)
        self.resize(WindowWidth,WindowHeight)
        self.figure = plt.figure(figsize=(14,4))    
        self.canvas = FigureCanvas(self.figure) 
        l.addWidget(self.canvas, 0, 0, 9, (100-4))
        toolbar = NavigationToolbar(self.canvas, self)
        l.addWidget(toolbar)

# transfer QtextEdit to be clickable and add "clicked"        
class ClickableTextEdit(QtWidgets.QTextEdit):
    clicked= QtCore.pyqtSignal()
    
    def mousePressEvent(self, e):
        #self.viewport().setCursor(Qt.PointingHandCursor)
        self.link = self.anchorAt(e.pos())
    
    def mouseReleaseEvent(self, e):
        #self.viewport().setCursor(Qt.IBeamCursor)
        if self.link:
            print(f"Clicked on {self.link}")
            self.clicked.emit()
        
# Main window of Analysis_tool (AT):            
class Ui_MainWindow(object):
    def setupUi(self, MainWindow):
        MainWindow.setObjectName("MainWindow")
        MainWindow.resize(600, 816)
            
        self.centralwidget = QtWidgets.QWidget(MainWindow)
        self.centralwidget.setObjectName("centralwidget")
        self.gridLayout = QtWidgets.QGridLayout(self.centralwidget)
        self.gridLayout.setObjectName("gridLayout")
        
        self.file_model = QtWidgets.QLineEdit(self.centralwidget)
        self.file_model.setObjectName("file_model")
        self.gridLayout.addWidget(self.file_model, 0, 0, 1, 2)
        
        self.btn_load_model = QtWidgets.QPushButton(self.centralwidget)
        self.btn_load_model.setObjectName("btn_load_model")
        #self.btn_load_model.clicked.connect(self.open_model)
        self.gridLayout.addWidget(self.btn_load_model, 0, 2, 1, 2)
        
        self.file_network = QtWidgets.QLineEdit(self.centralwidget)
        self.file_network.setObjectName("file_network")
        self.gridLayout.addWidget(self.file_network, 1, 0, 1, 2)
        
        self.btn_load_network = QtWidgets.QPushButton(self.centralwidget)
        self.btn_load_network.setObjectName("btn_load_network")
        self.gridLayout.addWidget(self.btn_load_network, 1, 2, 1, 2)
        
        
        self.file_binding = QtWidgets.QLineEdit(self.centralwidget)
        self.file_binding.setObjectName("file_binding")
        self.gridLayout.addWidget(self.file_binding, 2, 0, 1, 2)
        
        self.btn_load_binding = QtWidgets.QPushButton(self.centralwidget)
        self.btn_load_binding.setObjectName("btn_load_binding")
        self.gridLayout.addWidget(self.btn_load_binding, 2, 2, 1, 2)
        
        
        self.label_par = QtWidgets.QLabel(self.centralwidget)
        self.label_par.setObjectName("label_par")
        self.gridLayout.addWidget(self.label_par, 3, 0, 1, 1)
        
        self.label_mpl = QtWidgets.QLabel(self.centralwidget)
        self.label_mpl.setObjectName("label_par")
        self.gridLayout.addWidget(self.label_mpl, 3, 2, 1, 3)
        
        
        self.parameters = QtWidgets.QTextEdit(self.centralwidget)
        self.parameters.setPlaceholderText("")
        self.parameters.setObjectName("parameters")
        self.gridLayout.addWidget(self.parameters, 4, 0, 1, 2)
        
        self.mpl_xlim = QtWidgets.QLineEdit(self.centralwidget)
        self.mpl_xlim.setObjectName("xlim")
        self.mpl_xlim.setToolTip("xlim of plot.")
        self.mpl_xlim.setText('xlim=(1,1e7)')
        self.gridLayout.addWidget(self.mpl_xlim, 4, 2, 1, 2)
        
        self.mpl_xobs = QtWidgets.QLineEdit(self.centralwidget)
        self.mpl_xobs.setObjectName("xobs")
        self.mpl_xobs.setToolTip("Observed abundance. '0' means no observation.")
        self.mpl_xobs.setText('xobs=0')
        self.gridLayout.addWidget(self.mpl_xobs, 4, 2, 3, 1)
        
        self.mpl_Rdb = QtWidgets.QLineEdit(self.centralwidget)
        self.mpl_Rdb.setObjectName("Rdb")
        self.mpl_Rdb.setToolTip("""
        Ratio between diffusion and binding energy.
        '0' means using value from model file.
        Note: for model file from old version, the 'Rdb' is missed.
        So set here.
        """)
        self.mpl_Rdb.setText('Rdb=0')
        self.gridLayout.addWidget(self.mpl_Rdb, 4, 3, 3, 1)
        
        self.label_species = QtWidgets.QLabel(self.centralwidget)
        self.label_species.setObjectName("label_species")
        self.gridLayout.addWidget(self.label_species, 5, 0, 1, 1)
        
        self.species = QtWidgets.QLineEdit(self.centralwidget)
        self.species.setObjectName("species")
        self.gridLayout.addWidget(self.species, 5, 1, 1, 1)
        
        self.btn_analyze = QtWidgets.QPushButton(self.centralwidget)
        self.btn_analyze.setAutoFillBackground(False)
        self.btn_analyze.setObjectName("btn_analyze")
        self.btn_analyze.setStyleSheet("background-color : orange")
        self.btn_analyze.clicked.connect(self.do_analysis)
        self.gridLayout.addWidget(self.btn_analyze, 5, 2, 1, 2)
        
        self.ntop = QtWidgets.QLineEdit(self.centralwidget)
        self.ntop.setObjectName("ntop")
        self.ntop.setToolTip("The number of reactions with top rates.")
        self.gridLayout.addWidget(self.ntop, 5, 2, 3, 2)
        
        self.age = QtWidgets.QLineEdit(self.centralwidget)
        self.age.setObjectName("age")
        self.gridLayout.addWidget(self.age, 6, 1, 1, 1)
        
        
        self.label_age = QtWidgets.QLabel(self.centralwidget)
        self.label_age.setObjectName("label_age")
        self.gridLayout.addWidget(self.label_age, 6, 0, 1, 1)
    
                    
        self.label_formation = QtWidgets.QLabel(self.centralwidget)
        self.label_formation.setObjectName("label_formation")
        self.gridLayout.addWidget(self.label_formation, 7, 0, 1, 2)
        
        self.formation_textEdit = ClickableTextEdit(self.centralwidget) #QtWidgets.QTextEdit(self.centralwidget)
        self.formation_textEdit.setPlaceholderText("")
        self.formation_textEdit.setObjectName("formation")
        self.formation_textEdit.setToolTip("Click blue species to analyze it.")
        self.formation_textEdit.clicked.connect(self.update_species_formation)
        self.formation_textEdit.clicked.connect(self.do_analysis)
        self.gridLayout.addWidget(self.formation_textEdit, 8, 0, 1, 4)
        
        self.label_destruction = QtWidgets.QLabel(self.centralwidget)
        self.label_destruction.setObjectName("label_destruction")
        self.gridLayout.addWidget(self.label_destruction, 9, 0, 1, 1)
        
        self.destruction_textEdit = ClickableTextEdit(self.centralwidget) #QtWidgets.QTextEdit(self.centralwidget)
        self.destruction_textEdit.setObjectName("destruction_textEdit")
        self.destruction_textEdit.setToolTip("Click blue species to analyze it.")
        self.destruction_textEdit.clicked.connect(self.update_species_destruction)
        self.destruction_textEdit.clicked.connect(self.do_analysis)
        self.gridLayout.addWidget(self.destruction_textEdit, 10, 0, 1, 4)
        
        MainWindow.setCentralWidget(self.centralwidget)
        self.menubar = QtWidgets.QMenuBar(MainWindow)
        self.menubar.setGeometry(QtCore.QRect(0, 0, 800, 22))
        self.menubar.setObjectName("menubar")
        MainWindow.setMenuBar(self.menubar)
        self.statusbar = QtWidgets.QStatusBar(MainWindow)
        self.statusbar.setObjectName("statusbar")
        MainWindow.setStatusBar(self.statusbar)

        self.retranslateUi(MainWindow)
        QtCore.QMetaObject.connectSlotsByName(MainWindow)

    def retranslateUi(self, MainWindow):
        _translate = QtCore.QCoreApplication.translate
        MainWindow.setWindowTitle(_translate("MainWindow", "Analysis tool for GGCHEMPY"))
        MainWindow.setWindowIcon(QtGui.QIcon('icon/icon_AT.png'))
        self.file_binding.setText(_translate("MainWindow", "Will loaded from model file. Press 'Analyze'"))
        self.file_network.setText(_translate("MainWindow", "Will loaded from model file. Press 'Analyze'"))
        self.file_model.setText(_translate("MainWindow", "out/TMC1.dat"))
        self.parameters.setText("Will be loaded from model file.")
    
        self.label_age.setText(_translate("MainWindow", "Age:"))
        self.label_par.setText(_translate("MainWindow", "Model parameters:"))
        self.label_mpl.setText(_translate("MainWindow", "Other parameters:"))
        self.label_destruction.setText(_translate("MainWindow", "Destruction (clickable):"))
        
        self.btn_load_model.setText(_translate("MainWindow", "Load model file"))
        self.btn_load_network.setText(_translate("MainWindow", "Load network file"))
        
        self.label_formation.setText(_translate("MainWindow", "Formation (clickable):"))
        self.btn_load_binding.setText(_translate("MainWindow", "Load binding energy file"))
        self.label_species.setText(_translate("MainWindow", "Species:"))
        self.species.setText(_translate("MainWindow", "CO"))
        self.age.setText(_translate("MainWindow", "1e5"))
        self.btn_analyze.setToolTip(_translate("MainWindow", "Click to analyze"))
        self.btn_analyze.setText(_translate("MainWindow", "Analyze"))
        
        self.formation_textEdit.setText('Please press "Analyze"')
        self.destruction_textEdit.setText('Please press "Analyze"')
        self.ntop.setText('ntop=5')
        
    def AT_message(self, message):
        """
        Analysis tool (AT) message when some errors occur
        """
        #QtWidgets.QMessageBox.about(self, "title", message)
        msg = QtWidgets.QMessageBox()
        msg.setIcon(QtWidgets.QMessageBox.Critical)
        msg.setText(message)
        #msg.setInformativeText(e)
        msg.setWindowTitle("Error")
        msg.exec_()
        
    def do_analysis(self):
        """"
        Collect all information and do the analysis
        """
        Rdb= float(self.mpl_Rdb.text().split('=')[1])
        gg = ggmodel()
        
        gg.load_a_model( self.file_model.text(), self.file_network.text(), self.file_binding.text(), Rdb )
        
        pars = 'Current parameters from "%s" are listed below:\n'% self.file_model.text()
        pars+= ' ggchempy.gas.nH = %10.3g \n' % ggchempy.gas.nH
        pars+= ' ggchempy.gas.Av = %10.3g \n '% ggchempy.gas.Av
        pars+= ' ggchempy.gas.Chi= %10.3g \n '% ggchempy.gas.Chi
        pars+= ' ggchempy.gas.T = %10.3g \n '% ggchempy.gas.Zeta
        pars+= ' ggchempy.gas.T = %10.3g \n '% ggchempy.gas.T   
        pars+= ' ggchempy.dust.surface.Rdb=%10.3f\n'%ggchempy.dust.surface.Rdb
        pars+= ' ggchempy.dust.radius = %10.3g \n '%ggchempy.dust.radius
        pars+= ' ggchempy.dust.dite_density = %10.3g \n '%ggchempy.dust.site_density
        pars+= ' ggchempy.dust.rho = %10.3g \n '%ggchempy.dust.rho
        pars+= ' ggchempy.dust.nd = %10.3g \n '%ggchempy.dust.nd
        pars+= ' ggchempy.dust.T = %10.3g \n '%ggchempy.dust.T
        pars+= ' ggchempy.dust.mass = %10.3g \n '%ggchempy.dust.mass 
        pars+= ' ggchempy.dust.stick0 = %10.3g \n '%ggchempy.dust.stick0 
        pars+= ' ggchempy.dust.fc = %10.3g \n '%ggchempy.dust.fc 
        pars+= ' ggchempy.dust.Tpeak  = %10.3g \n '%ggchempy.dust.Tpeak 
        pars+= ' ggchempy.ggpars.d2gmr = %10.3g \n '%ggchempy.ggpars.d2gmr 
        pars+= ' ggchempy.ggpars.ns = %d \n '%ggchempy.ggpars.ns 
        pars+= ' ggchempy.ggpars.ne = %d \n '%ggchempy.ggpars.ne 
        pars+= ' ggchempy.ggpars.nr  = %d '%ggchempy.ggpars.nr 
        
        self.parameters.setText(pars)
        self.file_network.setText(ggchempy.ggpars.ggfiles[0])
        self.file_binding.setText(ggchempy.ggpars.ggfiles[1])
        

        sp = self.species.text().upper()
        if sp=='E-': sp='E'
   
        if sp in gg.model.keys():
            ### define the number of reactions with top rates
            ntop_init= int(self.ntop.text().split('=')[1]) ### inital number of top reactions to show
            ntop = ntop_init  ## show on the plot
            ntop_to_text=10   ## show on the QTextEdit
            
            ### get rates of a given species (sp) at a age
            ifor, ides, rates_for, rates_des, reactions, idx_for, idx_des, percent_for_at_age,  percent_des_at_age = \
            gg.analyze( sp, float(self.age.text()), ntop=ntop)
            """
            Notes for results:
            ifor     ---- resorted index of formation reactions
            ides     ---- resorted index of destruction reactions
            rates_for --- rates for formation
            rates_des --- rates for destruction
            reactions --- ggchempy.reactions (all)
            idx_for   --- index of formation reactions
            idx_des   --- index of destruction reactions
            """
            
            nsfor = np.shape(rates_for)
            
            ### open a child window to show plot:
            self.mySubwindow=subwindow()
            self.mySubwindow.createWindow(1000, 600)
            
            if OS=='Windows':
                desktop=QtWidgets.QApplication.desktop()
                self.mySubwindow.move( int(0.01*desktop.width()), int(0.1*desktop.height()))
            
            xlim_temp = self.mpl_xlim.text().split('=')[1]
            xlim_temp = (xlim_temp.replace('(','').replace(')','')).split(',')
            xlim = [float(xlim_temp[0]), float(xlim_temp[1])]
            
            xobs = float( self.mpl_xobs.text().split('=')[1] )

            self.formation_textEdit.clear()
            self.destruction_textEdit.clear()
            
            spec = sp
            age  = float(self.age.text())
            if sp in gg.model.keys():
                #### Abundance as function of time
                axes=self.mySubwindow.figure.add_subplot(221)
                axes.plot(gg.model['time'], gg.model[spec]/gg.nH, label=self.species.text())
                axes.legend(loc='upper left')
                axes.set_xlim(xlim[0],xlim[1])
                axes.set_yscale('log')
                axes.set_xscale('log')
                axes.axvline(age, color='r', linestyle=':')
                if xobs!=0.0: axes.axhline(xobs, color='r', linestyle=':')
                axes.set_xlabel('Time (yr)')
                axes.set_ylabel('Abundance')
            
                def plot_rates(axes, isorted, rates_x, idx_reac, ntop=5, xtype='for'):
                    note=''
                    if ntop>isorted.size: 
                        ntop=isorted.size
                        note='Found reactions with number < ntop='+str(ntop_init)+'. Adjusted to '+str(ntop)+'.'
                    percent=[]
                    reac_list=''
                    reac_list_click=''
                    ### loop to plot n-top reactions:
                    for i in range(ntop):
                        ii=isorted[- (i+1) ]  ### incresing order. so using -index
                        reac, reac_click = format_reaction(reactions.spec[int(idx_reac[ii]),:])
                        if xtype=='for':
                            reac_list += str(i+1)+': ('+ "%8.3f"%(percent_for_at_age[ii]) +'% ) ' + reac + '\n'
                            reac_list_click += str(i+1)+': ('+ "%8.3f"%(percent_for_at_age[ii]) +'% ) ' + reac_click + '\n'
                            self.formation_textEdit.append(str(i+1)+': ('+ "%8.3f"%(percent_for_at_age[ii]) +'% ) ' + reac_click + '\n')
                            axes.plot(gg.model['time'], rates_x[ii,:], label= '('+ "%6.3f"%(percent_for_at_age[ii]) +'% )' +reac)
                            percent.append(percent_for_at_age[ii])
                        if xtype=='des':
                            reac_list += str(i+1)+': ('+ "%8.3f"%(percent_des_at_age[ii]) +'% ) ' + reac + '\n'
                            reac_list_click += str(i+1)+': ('+ "%8.3f"%(percent_des_at_age[ii]) +'% ) ' + reac_click + '\n'
                            self.destruction_textEdit.append(str(i+1)+': ('+ "%8.3f"%(percent_des_at_age[ii]) +'% ) ' + reac_click + '\n')
                            axes.plot(gg.model['time'], rates_x[ii,:], label= '('+ "%6.3f"%(percent_des_at_age[ii]) +'% )' +reac)
                            percent.append(percent_des_at_age[ii])
                        if i==0:
                            ytop = np.interp(age, gg.model['time'], rates_x[ii,:])
                            axes.set_ylim(ytop*10/1e6, ytop*10)
                    
                    ### continue to loop to add more reactions (<=ntop_to_text) on the QTextEdit        
                    if ntop_to_text>ntop:
                        ntop_loop = ntop_to_text
                        if ntop_to_text>isorted.size: ntop_loop=isorted.size
                        for i in range(ntop, ntop_loop, 1):
                            ii=isorted[- (i+1) ]  ### incresing order. so using -index
                            reac, reac_click = format_reaction(reactions.spec[int(idx_reac[ii]),:])
                            if xtype=='for':
                                reac_list += str(i+1)+': ('+ "%8.3f"%(percent_for_at_age[ii]) +'% ) ' + reac + '\n'
                                reac_list_click += str(i+1)+': ('+ "%8.3f"%(percent_for_at_age[ii]) +'% ) ' + reac_click + '\n'
                                self.formation_textEdit.append(str(i+1)+': ('+ "%8.3f"%(percent_for_at_age[ii]) +'% ) ' + reac_click + '\n')
                                percent.append(percent_for_at_age[ii])
                            if xtype=='des':
                                reac_list += str(i+1)+': ('+ "%8.3f"%(percent_des_at_age[ii]) +'% ) ' + reac + '\n'
                                reac_list_click += f'<a '+str(i+1)+': ('+ "%8.3f"%(percent_des_at_age[ii]) +'% ) ' + reac_click + '/a>'
                                self.destruction_textEdit.append(str(i+1)+': ('+ "%8.3f"%(percent_des_at_age[ii]) +'% ) ' + reac_click + '\n')
                                percent.append(percent_des_at_age[ii])
                            
                    if note!='':
                        reac_list += note
                        
                    if xtype=='for':
                        #self.formation_textEdit.setText(reac_list)    
                        axes.set_ylabel('Formation rate')
                    if xtype=='des':
                        #self.destruction_textEdit.setText(reac_list)  
                        axes.set_ylabel('Destruction rate')
                    axes.legend(loc='lower left')
                    axes.set_yscale('log')
                    axes.set_xscale('log')
                    axes.set_xlim(xlim[0],xlim[1])
                    axes.set_xlabel('Time (yr)')
                    
                    axes.axvline(age, color='r', linestyle=':')
                    
                    return reac_list, percent
                
                ### formation    
                axes=self.mySubwindow.figure.add_subplot(223)
                reac_list_for, percent_for = plot_rates(axes, ifor, rates_for, idx_for, ntop=ntop, xtype='for')
                
                ### destruction
                axes=self.mySubwindow.figure.add_subplot(224)
                reac_list_des, percent_des = plot_rates(axes, ides, rates_des, idx_des, ntop=ntop, xtype='des')
                
                ### simple network using networkx
                axes= self.mySubwindow.figure.add_subplot(222)
                G = nkx.DiGraph()
                
                edge_labels={}
                def add_egde_to_G(G, reac_list_xxx, percent_xxx):
                    ireac=0
                    irx=0
                    for reacx in reac_list_xxx.split('\n'):
                        try:
                            if reacx.strip()!="":
                                reac = reacx.split(')')
                                reacy = reac[1].split('->')
                                for r in reacy[0].split(' + '):
                                    for p in reacy[1].split(' + '):
                                        r=r.replace('e-','E').upper(); p=p.replace('e-','E').upper()
                                        if G.has_edge(r.strip(),p.strip())==False: 
                                            if p.strip()==sp:
                                                G.add_edge(r.strip(),p.strip(), color='#0597ED', weight=3.0)
                                                edge_labels[(r.strip(),p.strip())]="%6.3f"%(percent_xxx[irx])
                                            if r.strip()==sp:
                                                G.add_edge(r.strip(),p.strip(), color='#F4C2C2')
                                                edge_labels[(r.strip(),p.strip())]= "%6.3f"%(percent_xxx[irx])
                            ireac+=1
                        except:
                            pass
                        if ireac>2:
                            break
                        irx+=1
                
                add_egde_to_G(G, reac_list_for, percent_for)
                add_egde_to_G(G, reac_list_des, percent_des)
                
                ### transfer upper case to normal case e.g. SI -> Si
                spx = sp
                for el in elements:
                    if el in spx:
                        spx = spx.replace(el, elements[el])
                
                
                re_labels={}
                for node in G:
                    for el in elements:
                        if el in node:
                            re_labels[node] = node.replace(el, elements[el])
                nkx.relabel_nodes(G, re_labels, copy=False)
                if 'E' in G:
                    nkx.relabel_nodes(G, {'E':'e-'}, copy=False)     
                
                if spx=='E': spx='e-'
                color_map=[]
                for node in G:
                    if node==spx:
                        color_map.append('#0597ED')
                    else: 
                        color_map.append('#F4C2C2')
                                       
                pos = nkx.circular_layout(G) #nkx.spring_layout(G)#
                if OS=='Linux':
                    node_size = 500
                else:
                    node_size = 200
                
                edge_colors = nkx.get_edge_attributes(G,'color').values()
                #weights = nkx.get_edge_attributes(G, "weight").values()
                #edges = G.edges()
                #weights = [G[u][v][weight'] for u,v in edges]
                
                nkx.draw(G, pos, node_color=color_map, with_labels = True, node_shape="o", \
                         node_size=node_size, font_size=font_size, edge_color=edge_colors)#, width=weights )#'#0597ED', font_size=8, font_color="k")
                         
                edge_labels_2={}
                for lab in edge_labels:
                    A,B=lab
                    for el in elements:
                        if el in A:
                            A = A.replace(el, elements[el])
                        if el in B:
                            B = B.replace(el, elements[el])
                    if A=='E': A='e-'
                    if B=='E': B='e-'
                    edge_labels_2[(A,B)]=edge_labels[lab]
                    
                nkx.draw_networkx_edge_labels(G, pos, edge_labels=edge_labels_2, font_color='k', font_size=6)
                """
                edge_labels=={('A', 'B'): 'AB', \
                              ('B', 'C'): 'BC', \
                              ('B', 'D'): 'BD'},\
                """

                self.mySubwindow.canvas.draw()
            else:
                axes=self.mySubwindow.figure.add_subplot(111)
                axes.set_title(spec+' does not exist in current model.', color='red')

            #plt.tight_layout(self.mySubwindow.axes)
            plt.subplots_adjust(left=0.1, top=0.98,right=0.98, bottom=0.12, hspace=0.25, wspace=0.25)
            self.mySubwindow.show()
        else:
            self.AT_message(sp+" does not exist in current model. Try another one.")
        
    def update_species_formation(self):
        """
        To update the species when a clickable text is pressed.
        """
        self.species.setText(self.formation_textEdit.link)
    def update_species_destruction(self):
        """
        To update the species when a clickable text is pressed.
        """
        self.species.setText(self.destruction_textEdit.link)
            
### ui to window          
class Window(Ui_MainWindow, QtWidgets.QMainWindow):
    def __init__(self):
        super(Window, self).__init__()
        self.setupUi(self)
        self.btn_load_model.clicked.connect(self.open_model)
    def open_model(self):
        filepath, filetype = self.model_file_Name = QtWidgets.QFileDialog.getOpenFileName(self, 'Open file', 'All Files (*.*)')
        self.file_model.setText(filepath.replace(os. getcwd()+'/',''))
    
if __name__=='__main__':
    import sys
    app=QtWidgets.QApplication(sys.argv)
    if OS=='Darwin': 
        path = os.path.join(os.path.dirname(sys.modules[__name__].__file__), 'icon/icon_AT.png')
        app.setWindowIcon(QtGui.QIcon(path))
    else:
        app.setWindowIcon(QtGui.QIcon('icon/icon_AT.png'))
    mywindow = Window()
    MainWindow = QtWidgets.QMainWindow()
    if OS in ['Linux', 'Darwin']:
        desktop=QtWidgets.QApplication.desktop()
        mywindow.move(desktop.width(), desktop.height()) ### position on the screen.
    else:
        desktop=QtWidgets.QApplication.desktop()
        mywindow.move( int(0.6*desktop.width()), int(0.1*desktop.height())) ### position on the screen.
    mywindow.show()
    sys.exit(app.exec_())    

