import sys
import csv
from PyQt5.QtCore import *
from PyQt5.QtGui import *
from PyQt5.QtWidgets import *

from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
#import matplotlib.pyplot as plt
#import numpy as np
import math

from navierstokes import *

# **** Implement Algorithms ****
class FF (FigureCanvas):
    def __init__(self,parent=None):
        self.fig=plt.figure(figsize=(4,3),dpi=80,facecolor='white')
         
        FigureCanvas.__init__(self,self.fig)
        self.setParent(parent)
        
        self.us=np.zeros((ny,nx))
        self.vs=np.zeros((ny,nx))
        self.ps=np.zeros((ny,nx))
        self.X,self.Y=np.meshgrid(x,y)
        # fit plot size to figure frame
        self.compute()
        self.plotfig()
        
    # compute the results of u,v,p
    def compute(self): 
        self.us,self.vs,self.ps=navierstokes(u,v,p)
    
    # prepare fig plot to show
    def plotfig(self):
        
        plt.contour(self.X,self.Y,self.us,cmap=cm.coolwarm,vmin=0,vmax=3.5)
        plt.contourf(self.X,self.Y,self.us,cmap=cm.coolwarm,vmin=0,vmax=3.5)
        plt.colorbar()
        plt.clim(0,3.5)
        plt.quiver(self.X[::3,::3],self.Y[::3,::3],self.us[::3,::3],self.vs[::3,::3])
        plt.tight_layout()

    # responding to the click of button,accept new parameters, update plots
    def updatefig(self):
        self.compute()
        plt.contour(self.X,self.Y,self.us,cmap=cm.coolwarm,vmin=0,vmax=3.5)
        plt.contourf(self.X,self.Y,self.us,cmap=cm.coolwarm,vmin=0,vmax=3.5)       
        plt.quiver(self.X[::3,::3],self.Y[::3,::3],self.us[::3,::3],self.vs[::3,::3])        
        self.draw() 

class App(QMainWindow):
    def __init__(self):
        QMainWindow.__init__(self)
        self.setWindowTitle('Channel Flow - Navier-Stokes')
        self.setGeometry(50,50,800,600)

        # **** Create Menu ****
        self.menu=self.menuBar()
        self.file=self.menu.addMenu('File')
        self.edit=self.menu.addMenu('Edit')
        self.help=self.menu.addMenu('Help')

        self.exit=QAction(QIcon(''),'Exit',self)
        self.exit.setShortcut('Ctrl+Q')
        self.exit.triggered.connect(self.close)
        self.file.addAction(self.exit)

        self.about=QAction(QIcon(''),'About',self)
        self.about.setShortcut('Ctrl+T')
        self.about.triggered.connect(self.aboutact)
        self.help.addAction(self.about)

        # **** Create 4 Tabs to show Equations, Codes and GUI ****
        self.mainW1=QWidget(self)
        self.layout1=QVBoxLayout(self.mainW1)
        
        # initiate tabs
        self.tabs=QTabWidget()
        self.tab1=QWidget()
        self.tab2=QWidget()
        self.tab3=QWidget()
        self.tab4=QWidget()
        # add tab1 and tab2 to tabs
        self.tabs.addTab(self.tab1,'Introduction')
        self.tabs.addTab(self.tab2,'Code 1: Navier-Stokes')
        self.tabs.addTab(self.tab3,'Code 2: GUI')
        self.tabs.addTab(self.tab4,'GUI')
	
        # create tab1 for intro with pdf file scroll bar      
        
        self.tab1_layout=QVBoxLayout(self)
        self.scrollArea=QScrollArea()
        self.pix1=QPixmap('eq.jpg')
        self.pix1r=self.pix1.scaledToWidth(700)            # rescale pixmap
        self.label1=QLabel(self)
        self.label1.setAlignment(Qt.AlignCenter)            # Center Alignment of Label
        self.label1.setPixmap(self.pix1r)
        self.scrollArea.setWidget(self.label1)
        self.tab1_layout.addWidget(self.scrollArea)
        self.tab1.setLayout(self.tab1_layout)
        
        # create tab2 for Code1 Scrollbar
        self.tab2_layout=QVBoxLayout(self)
        self.textedit1=QPlainTextEdit()
        self.text1=open('navierstokes.py').read()
        self.textedit1.setPlainText(self.text1)
        self.tab2_layout.addWidget(self.textedit1)
        self.tab2.setLayout(self.tab2_layout)

        # create tab3 for Code2 Scrollbar
        self.tab3_layout=QVBoxLayout(self)
        self.textedit2=QPlainTextEdit()
        self.text2=open('gui.py').read()
        self.textedit2.setPlainText(self.text2)
        self.tab3_layout.addWidget(self.textedit2)
        self.tab3.setLayout(self.tab3_layout)

        # create tab4 for GUI  (input textboxs and button and results)
        self.tab4_layout=QVBoxLayout(self)
        self.fig=FF()
        self.tab4_layout.addWidget(self.fig)
        self.label3=QLabel('Enter Number of Elements',self)
        self.tab4_layout.addWidget(self.label3)
        self.entry1=QLineEdit('11',self)
        self.entry1.resize(50,20)
        self.entry1.move(200,10)
        self.tab4_layout.addWidget(self.entry1) 
        self.but1=QPushButton('Plot',self)
        self.but1.clicked.connect(self.plotting)
        self.but1.resize(40,20)
        self.but1.move(270,300)
        self.tab4_layout.addWidget(self.but1) 
        self.tab4.setLayout(self.tab4_layout)

        # add tabs to layout1 and add layout1 to mainW1
        self.layout1.addWidget(self.tabs)
        self.layout1.setAlignment(Qt.AlignCenter)
        self.mainW1.setLayout(self.layout1)
        # adjust size and position of mainW1
        self.mainW1.resize(760,560)
        self.mainW1.move(10,30)


    @pyqtSlot()
    def aboutact(self):
        resp=QMessageBox.question(self,'About','This GUI is designed to demonstrate Navier-Stokes-based Channel Flow Analysis',QMessageBox.Ok)
     
    @pyqtSlot()
    def plotting(self):
        self.fig.updatefig()   # update plots
                                
if __name__=="__main__":
	qApp = QApplication(sys.argv)
	aw = App()
	aw.show()
	sys.exit(qApp.exec_())


