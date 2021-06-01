
# -*- coding: utf-8 -*-

from PyQt5.QtWidgets import (QWidget, QMainWindow, QTextEdit, 
    QAction, QFileDialog, QApplication, QPushButton,QCheckBox,
    QProgressBar,QGridLayout,QLabel,QLineEdit,QToolTip,
    QMessageBox,QPlainTextEdit)
import PyQt5.QtWidgets
from PyQt5.QtGui import QIcon,QFont,QPixmap  
from PyQt5.QtCore import Qt, pyqtSignal, QObject,QBasicTimer
from PyQt5 import QtCore
import sys
import numpy as np
from src.ggchemlib import gas,dust,ggchem,ggpars,iswitch


#class Communicate(QObject):
#    closeApp = pyqtSignal() 

class GGCHEMGUI(QMainWindow):

    def __init__(self):
        super().__init__()

        self.initUI()


    def initUI(self):   
        
        QToolTip.setFont(QFont('SansSerif', 10))
        self.setToolTip('GGCHEM version 1.0') 
        
        ### set common font:
        self.setFont(QFont('SansSerif', 10))
        
        """
        grid = QGridLayout()
        self.setLayout(grid)  

        self.c = Communicate()
        self.c.closeApp.connect(self.close)
        
        self.textEdit = QTextEdit()
        self.setCentralWidget(self.textEdit)
        """
        
        self.statusBar()
        key_bg_color='#ff6600'

        """
        openFile = QAction(QIcon('open.png'), 'Open', self)
        openFile.setShortcut('Ctrl+O')
        openFile.setStatusTip('Open new File')
        openFile.triggered.connect(self.showDialog)

        menubar = self.menuBar()
        fileMenu = menubar.addMenu('&File')
        fileMenu.addAction(openFile)    
        
        toolMenu = menubar.addMenu('&Tools')
        helpMenu = menubar.addMenu('&Help')
        """
        
        
        
        
        """
        #logo_pix = QPixmap('logo.png')
        self.gglogo = QLabel(self)
        self.gglogo.setGeometry(0, 0, 350, 30)
        self.gglogo.setScaledContents(True)
        self.gglogo.setStyleSheet("background-color: "+key_bg_color+"; border: 1px solid red")
        #border: 1px solid "+label_bg_color) 
        #self.gglogo.setPixmap(logo_pix)
        #self.gglogo.setAutoFillBackground(True)
        self.gglogo.setText('Gas-Grain CHEMistry model for interstellar medium.')
        self.gglogo.setToolTip('This is a two-phase (gas + dust grain surface) chemical model.')
        """
        
        """
        ### choose mode at top:
        self.btn_smode = QPushButton('single model', self)  
        #grid.addWidget(self.btn_smode) 
        #btn_smode.move(30,30)
        self.btn_smode.setGeometry(QtCore.QRect(20, 30, 150, 20))
        self.btn_smode.clicked.connect(self.buttonClicked) 
        self.btn_smode.setToolTip('Please set parameters below.')
        
        
        
        self.btn_mmode= QPushButton('multiple models', self) 
        #grid.addWidget(self.btn_mmode) 
        #btn_mmode.move(150,30)
        self.btn_mmode.setGeometry(QtCore.QRect(180, 30, 150, 20))
        self.btn_mmode.clicked.connect(self.buttonClicked)
        self.btn_mmode.setToolTip('Please use a Python script.')
        """
        
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
        self.finputnet.setText('/in/network2.txt')
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
        self.finputiab.setText('/in/iabun.txt')
        self.finputiab.setToolTip('file of initial abundances list.')
        #self.finputiab.textChanged[str].connect(self.onChanged)
        
        self.inputfed = QLabel(self)
        self.inputfed.setGeometry(20, 240, 170, 20)
        self.inputfed.setText('Binding energies:')
        
        self.finputed = QLineEdit(self)
        self.finputed.setGeometry(xmid, 240, xwid, 20)
        self.finputed.setText('/in/ed.txt')
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
        
        
        ### start and progress bar
        #self.pbar = QProgressBar(self)
        #self.pbar.setGeometry(30, 550, 200, 25)

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
        
        
        window = self.setGeometry(300, 300, 350, 600)
        self.setWindowTitle('GGCHEM version 1.0')
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
        #    self.setWindowTitle('GGCHEM versio ') 
    
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
        
        iswitch.iSS=0
        iswitch.iNTD=0
        iswitch.icollapse=0
        
        if self.cb_ss.isChecked(): 
            iswitch.iSS=1
        if self.cb_rd_g07.isChecked(): 
            iswitch.iNTD=1
        if self.cb_rd_m16.isChecked(): 
            iswitch.iNTD=2
        if self.cb_static.isChecked(): 
            iswitch.icollapse=0
        if self.cb_collap.isChecked(): 
            iswitch.icollapse=1 
            
            
        input_network = self.finputnet.text()
        input_iabun   = self.finputiab.text()
        input_ed      = self.finputed.text()
        ggpars.files  = [input_network, input_iabun, input_ed] 
        
        
        gas.T     = float( self.gas_T_value.text() )
        gas.Zeta  = float( self.gas_Zeta_value.text() )
        gas.Chi   = float( self.gas_Chi_value.text() )
        gas.Av    = float( self.gas_Av_value.text() )
        
        if iswitch.icollapse==0:
            gas.nH    = float( self.gas_nH_value.text() )
        elif iswitch.icollapse==1:
            nHxx     = self.gas_nH_value.text().split(',')
            gas.nH0 = float(nHxx[0])
            gas.nH1 = float(nHxx[1])
            gas.Av0 = gas.Av
        
        ggpars.d2gmr = float(self.d2gmr_value.text())
        
        dust.T    = float( self.dust_T_value.text() )
        dust.radius = float( self.dust_R_value.text() )
        dust.rho    = float( self.dust_rho_value.text() )
        dust.surface.Rdb  = float( self.dust_Rdb_value.text() )
        
        ggpars.rtol = (self.ode_rtol_value.text())
        ggpars.atol = (self.ode_atol_value.text())
        
        modelname = self.output_modelname_value.text()
        ggpars.tf = float(self.output_finalage_value.text())
        ggpars.nt = int(self.output_timestep_value.text())
        
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
        ggchem.run(modelname)
        #end_time = time.time()
        
        self.output_info.setText('CPU time =%10.3g (s)'%(ggpars.cputime))
        
        self.setWindowTitle( 'GGCHEM version 1.0')
        
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

        
"""            
if __name__ == '__main__':

    app = QApplication(sys.argv)
    ex = GGCHEMGUI()
    sys.exit(app.exec_())
"""
