import sys
import csv
from PyQt5.QtWidgets import *
from PyQt5.QtCore import *
from PyQt5.QtGui import *

from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
import matplotlib.pyplot as plt
import numpy as np
import math
from scipy.linalg import eigh

# **** Implement Algorithms ****
class FF (FigureCanvas):
    def __init__(self,parent=None,num_elems=10):
        self.fig=plt.figure(figsize=(8,3),dpi=80,facecolor='white')
         
        FigureCanvas.__init__(self,self.fig)
        self.setParent(parent)

        self.nelem=num_elems
        # store M, K, freq, evecs
        self.result=[]
        
        for self.counter in range(1,self.nelem+1):
             M, K, freq, evecs=self.compute()
             exact_freq=math.pi/2
             error=(freq[0]-exact_freq)/exact_freq *100
             self.result.append((self.counter,error,freq[0]))

        self.elements = np.array([x[0] for x in self.result])
        self.err = np.array([x[1] for x in self.result])
        self.freq = np.array([x[2] for x in self.result])
        
        # convert result(element no, error, frequency) to an array
        self.storearr=np.array(self.result)       

	
        # plot elements vs error
        self.ax1=self.fig.add_subplot(121)
        self.ax1.plot(self.elements,self.err,linestyle='--',color='blue')
        plt.xlabel('number of elements')
        plt.ylabel('error (%)')
        # fit plot size to figure frame
        plt.tight_layout()

        self.ax2=self.fig.add_subplot(122)
        self.ax2.plot(self.elements,self.freq,color='red')
        plt.xlabel('number elements')
        plt.ylabel('frequency (Hz)')
        plt.tight_layout()
        

    # compute M, K, freq, evecs
    def compute(self):

        # element mass and stiffness matrices for a bar
        m = np.array([[2,1],[1,2]])/(6. * self.counter)        # define mass matrix
        k = np.array([[1,-1],[-1,1]])*float(self.counter)      # define stiffness matrix

        # construct global mass and stiffness matrices
        M = np.zeros((self.counter+1,self.counter+1))
        K = np.zeros((self.counter+1,self.counter+1))

        # assembly of elements
        for i in range(self.counter):
            M_temp = np.zeros((self.counter+1,self.counter+1))
            K_temp = np.zeros((self.counter+1,self.counter+1))
            M_temp [i:i+2,i:i+2]=m
            K_temp [i:i+2,i:i+2]=k
            M += M_temp
            K += K_temp

        # reomove the fixed degree of freedom
        for i in range(2):
            M = np.delete(M, 0, axis=i)  #axis=0 row, 1 column
            K = np.delete(K, 0, axis=i)

        # eigenvalue problem
        evals, evecs = eigh(K,M)
        freq = np.sqrt(evals)
        return M, K, freq, evecs
    
    # generate plots
    def drawing(self):
        self.result=[]
        for self.counter in range(1,self.nelem+1):
            M, K, freq, evecs=self.compute()
            exact_freq=math.pi/2
            error=(freq[0]-exact_freq)/exact_freq *100
            self.result.append((self.counter,error,freq[0]))

        self.elements = np.array([x[0] for x in self.result])
        self.err = np.array([x[1] for x in self.result])
        self.freq = np.array([x[2] for x in self.result])
        # convert result(element no, error, frequency) to an array
        self.storearr=np.array(self.result)       

    
    # responding to the click of button,accept new number of elements, update plots
    def updatefig(self):
        self.drawing()
        self.ax1.plot(self.elements,self.err,linestyle='--',color='blue')
        self.ax2.plot(self.elements,self.freq, color='red')
        self.draw() 

class App(QMainWindow):
    def __init__(self):
        QMainWindow.__init__(self)
        self.setWindowTitle('Bar Element')
        self.setGeometry(50,50,500,500)

        # **** Create 2 Tabs to show Schematic and Equations ****
        self.mainW1=QWidget(self)
        self.layout1=QVBoxLayout(self.mainW1)
        # initiate tabs
        self.tabs=QTabWidget()
        self.tab1=QWidget()
        self.tab2=QWidget()
        # add tab1 and tab2 to tabs
        self.tabs.addTab(self.tab1,'Schematic')
        self.tabs.addTab(self.tab2,'Equations')
        # create tab1
        self.tab1_layout=QVBoxLayout(self)
        self.pix1=QPixmap('BarSchematic.jpg')
        self.pix1r=self.pix1.scaledToHeight(180)            # rescale pixmap
        self.label1=QLabel(self)
        self.label1.setAlignment(Qt.AlignCenter)            # Center Alignment of Label
        self.label1.setPixmap(self.pix1r)
        self.tab1_layout.addWidget(self.label1)
        self.tab1.setLayout(self.tab1_layout)
        # create tab2
        self.tab2_layout=QVBoxLayout(self)
        self.pix2=QPixmap('Eq.jpg')
        self.pix2r=self.pix2.scaledToHeight(180)
        self.label2=QLabel(self)
        self.label2.setAlignment(Qt.AlignCenter)
        self.label2.setPixmap(self.pix2r)
        self.tab2_layout.addWidget(self.label2)
        self.tab2.setLayout(self.tab2_layout)
        # add tabs to layout1 and add layout1 to mainW1
        self.layout1.addWidget(self.tabs)
        self.mainW1.setLayout(self.layout1)
        # adjust size and position of mainW1
        self.mainW1.resize(260,240)
        self.mainW1.move(10,30)

        # **** Create Menu ****
        self.menu=self.menuBar()
        self.file=self.menu.addMenu('File')
        self.edit=self.menu.addMenu('Edit')
        self.help=self.menu.addMenu('Help')

        self.exit=QAction(QIcon(''),'Exit',self)
        self.exit.setShortcut('Ctrl+Q')
        self.exit.triggered.connect(self.close)
        self.file.addAction(self.exit)

        self.exit=QAction(QIcon(''),'Save Table',self)
        self.exit.setShortcut('Ctrl+S')
        self.exit.triggered.connect(self.saveTable)
        self.file.addAction(self.exit)

        self.about=QAction(QIcon(''),'About',self)
        self.about.setShortcut('Ctrl+T')
        self.about.triggered.connect(self.aboutact)
        self.help.addAction(self.about)

        # **** Create Plot *****
        self.mainW2=QWidget(self)
        self.layout2=QVBoxLayout(self)
        self.fig=FF()
        self.layout2.addWidget(self.fig)
        self.mainW2.setLayout(self.layout2)
        self.mainW2.move(5,260)
        self.mainW2.resize(495,240)

        # **** Create Table to show element no./ frequency / error *****
        #   self.mainW3=QWidget(self)
        self.mainW3=QWidget(self)
        self.createTable(self.fig.nelem,2)
        self.layout3=QVBoxLayout(self)
        self.layout3.addWidget(self.table)
        self.mainW3.setLayout(self.layout3)
        self.mainW3.move(270,40)
        self.mainW3.resize(220,230)

        # **** Create QLabel to show instruction ****
        self.mainW4=QWidget(self)
        self.layout4=QVBoxLayout(self)
        label3=QLabel('Enter Number of Elements',self)
        self.layout4.addWidget(label3)
        self.mainW4.setLayout(self.layout4)
        self.mainW4.resize(200,40)
        self.mainW4.move(20,1)
        
        # **** Create TextInput to input number of elements ****
        self.entry1=QLineEdit('11',self)
        self.entry1.resize(50,20)
        self.entry1.move(200,10) 
        
        # **** Create Plot Button to update charts
        self.but1=QPushButton('Plot',self)
        self.but1.clicked.connect(self.plotting)
        self.but1.resize(40,20)
        self.but1.move(270,10)


    @pyqtSlot()
    def aboutact(self):
        resp=QMessageBox.question(self,'About','This APP is designed to set different number of finite elements to compute the vibration frequency of a bar element',QMessageBox.Ok)

    @pyqtSlot()
    def saveTable(self):
        path,_=QFileDialog.getSaveFileName(self,'Save File','untitled.csv','CSV(*.csv)')
        if path:
             with open (path,'w') as stream:
                  writer=csv.writer(stream)     # set up stream to be passed to csv object writer
                  rowdata=[]
                  for row in range(self.table.rowCount()):
                      for column in range(self.table.columnCount()):
                          itm=self.table.item(row,column)              # retrieve item
                          rowdata.append(str(itm.text()))
                      writer.writerow(rowdata)
                      rowdata=[]

    def createTable(self,r,c):
        self.table=QTableWidget()
        self.table.setRowCount(r)
        self.table.setColumnCount(c)
        temparr=np.delete(self.fig.storearr,0,1)
        self.table.setHorizontalHeaderLabels(('Frequency (Hz)','Error (%)'))
        temparr[:,0],temparr[:,1]=temparr[:,1],temparr[:,0].copy()
        for i in range(r):
            for j in range(c):
                self.table.setItem(i, j, QTableWidgetItem(str(round(temparr[i,j],2))))
         
    @pyqtSlot()
    def plotting(self):
        self.fig.nelem=int(self.entry1.text()) 
        self.fig.updatefig()   # update plots
        #update table
        temparr=np.delete(self.fig.storearr,0,1)
        self.table.setRowCount(self.fig.nelem)       
        for i in range(self.fig.nelem):
            for j in range(2):
               self.table.setItem(i,j,QTableWidgetItem(str(round(temparr[i,j],2))))
                        
qApp = QApplication(sys.argv)
aw = App()
aw.show()
sys.exit(qApp.exec_())


