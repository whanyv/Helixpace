from tkinter import SEL_FIRST
from PySide6.QtWidgets import *
from PySide6.QtCore import *
from PySide6.QtGui import *
from primer import *
import sys
from encode import *
import os
import pydicom
import numpy as np
from tqdm import tqdm
import sys
import glob
import matplotlib.pyplot as plt
#from all_function import *


class entrainWindow(QMainWindow):
    def __init__(self): 
        super(entrainWindow, self).__init__()

        palette = QPalette()
        palette.setColor(QPalette.ColorRole.Button, QColor(100, 53, 53))
        #palette.setC
        # color(QPalette.ColorRole.WindowText, Qt.GlobalColor.white)
        self.setStyleSheet("background-image: url(/Users/pro2.0/Desktop/ec552_project/HELIXPACE/blue_bg.png)")
        # self.setAutoFillBackground(True)
        self.setWindowIcon(QIcon('icon.png'))
        self.setWindowTitle("DNA Storage")
        self.setFixedSize(QSize(1000, 500))

        label = QLabel(self)
        icon = QPixmap('icon.png')
        label.setPixmap(icon)
        label.move(425,100)
        label.resize(icon.width(),icon.height())

        encoding = QPushButton(self)
        encoding.setText("Start to encode your data!")
        encoding.setStyleSheet("background-color : grey")
        encoding.move(440,350)
        encoding.adjustSize()
        encoding.clicked.connect(self.window2)

        decoding = QPushButton(self)
        decoding.setText("Start to decode your data!")
        decoding.setStyleSheet("background-color : grey")
        decoding.move(440,400)
        decoding.adjustSize()
        decoding.clicked.connect(self.decode_Window)

        exit = QPushButton(self)
        exit.setText("Exit")
        exit.setStyleSheet("background-color : grey")
        exit.move(50,400)
        exit.adjustSize()
        exit.clicked.connect(self.closeIt)

        self.createMenuBar()

    
    def window2(self):                                             # <===
        self.w = MainWindow()
        self.w.show()
        self.hide()

    def decode_Window(self):                                             # <===
        self.w = decodeWindow()
        self.w.show()
        self.hide()

    def createMenuBar(self):
        menuBar = self.menuBar()
        fileMenu = QMenu("&File", self)
        menuBar.addMenu(fileMenu)
        editMenu = menuBar.addMenu("&Edit")
        helpMenu = menuBar.addMenu("&Help")
    
    def closeIt(self): 
        self.close()


class MainWindow(QMainWindow):
    def __init__(self):

        super(MainWindow, self).__init__()


        self.setWindowTitle("Encoding...")
        self.setFixedSize(QSize(1000, 500))
        self.setStyleSheet("background-image: url(/Users/pro2.0/Desktop/ec552_project/HELIXPACE/blue_bg.png)")
        # self.setAutoFillBackground(True)
        #self.setWindowIcon(QIcon('icon.png'))


        usernameBox = QLabel(self)
        usernameBox.setText("Username")
        usernameBox.move(165,70)

        passwordBox = QLabel(self)
        passwordBox.setText("Password")
        passwordBox.move(165,120)

        moreBox = QLabel(self)
        moreBox.setText("More personal info")
        moreBox.setFixedSize(200,40)
        moreBox.move(140,170)

        hospital_info = QLabel(self)
        hospital_info.setText("Enter your public primer section")
        hospital_info.setFixedSize(200,40)
        hospital_info.move(100,330)

        self.IDinput = QLineEdit(self)
        self.IDinput.move(300,70)
        self.IDinput.setMaxLength(12)
        self.IDinput.setFixedWidth(500)
        self.IDinput.setPlaceholderText("Choose a username that is 8-12 characters long and contains letters or numbers.")
        self.IDinput.setClearButtonEnabled(True)

        self.pwinput = QLineEdit(self)
        self.pwinput.move(300,120)
        self.pwinput.setMaxLength(12)
        self.pwinput.setFixedWidth(500)
        self.pwinput.setPlaceholderText("Choose a password that is 8-12 characters long and contains letters and numbers.")
        self.pwinput.setClearButtonEnabled(True)

        self.moreinfo = QTextEdit(self)
        self.moreinfo.move(300,170)
        self.moreinfo.setFixedWidth(500)
        self.moreinfo.setFixedHeight(150)
        self.moreinfo.setStyleSheet("QLineEdit"
                                "{"
                                "background : grey;"
                                "}")
        self.moreinfo.textChanged.connect(self.written_info)

        self.pub_primer = QTextEdit(self)
        self.pub_primer.move(300,340)
        self.pub_primer.setFixedSize(350,40)
        self.pub_primer.textChanged.connect(self.written_primer)

        #self.moreinfo.textBackgroundColor("background : lightblue")
        # self.moreinfo.setClearButtonEnabled(True)
        # self.moreinfo.setDragEnabled(True)

        enter = QPushButton(self)
        enter.setText("Enter!")
        button2 = QPushButton(self)
        button2.setText("Return")
        button3 = QPushButton(self)
        button3.setText("Upload DCM image!")
        button4 = QPushButton(self)
        button4.setText("Upload from a local file!")

        enter.move(800,400)
        enter.adjustSize()
        enter.clicked.connect(self.enter_clicked)
        button2.move(50,400)
        button2.adjustSize()
        button2.clicked.connect(self.homepage)
        button3.move(500,400)
        button3.adjustSize()
        button3.clicked.connect(self.Upload_DCM)
        button4.move(660,350)
        button4.clicked.connect(self.upload_primer)
        button4.adjustSize()


    def homepage(self):                                             # <===
        self.w = entrainWindow()
        self.w.show()
        self.hide()

    def upload_primer(self):
        path = QFileDialog.getOpenFileName(self, 'Open a file', '','All Files (*.*)')
        if path != ('', ''):
            print(path[0]) 
        self.hos_primer = os.path.basename(path[0])

    def written_primer(self):
        f = open('public_primer.txt', 'w')
        f.write(str(self.pub_primer.toPlainText))
        f.close()
        self.hos_primer = 'public_primer.txt'
    
    def written_info(self):
        self.personal_info = self.moreinfo.toPlainText()

    def enter_clicked(self):
        username = self.IDinput.text()
        password = self.pwinput.text()
        #personalinfo = self.moreinfo.text()
        seq = primerEncoder.primerencode(username, password)
        print(f"USERNAME:{self.IDinput.text()} \nPASSWORD:{self.pwinput.text()}")
        print(f"seq:{seq}")
        id_and_pw_check = primerEncoder.primercheck(seq)
        if id_and_pw_check == False:
            QMessageBox.information(self,'Warning!','Your username and password combination can not encode to effecient primer. Please try another combination.', QMessageBox.Yes)
        elif id_and_pw_check == True:
            if self.in_file.endswith('.dcm'):
                outfile_name, ok = QInputDialog.getText(self, 'encoding...','Name your DNA sequence file: \n(ex. if you want to name your DNA sequence file \'DNA\', enter : ./DNA.txt)')
                encoder = DNAEncoder()
                encoder.load_dcm(self.in_dir, self.in_file)
                encoder.load_hospital_primer(self.hos_primer)
                encoder.encode(outfile_name,username,password,self.personal_info)
            else:
                 QMessageBox.information(self,'Warning!','Please select .DCM files.', QMessageBox.Yes)

        


        #print(f"USERNAME:{self.IDinput.text()} \nPASSWORD:{self.pwinput.text()}")
        #ADD FUNCTIONS THAT TAKE USER ID AND PASSWORD

    # def closeIt(self): 
    #     self.close()

    def Upload_DCM(self):
        path = QFileDialog.getOpenFileName(self, 'Open a file', '','All Files (*.*)')
        if path != ('', ''):
            print(path[0]) 
        self.in_dir = os.path.dirname(path[0])
        self.in_file = os.path.basename(path[0])


class decodeWindow(QMainWindow):
    def __init__(self): 
        super(decodeWindow, self).__init__()

        self.setWindowTitle("Decoding...")
        self.setFixedSize(QSize(1000, 500))
        self.setStyleSheet("background-image: url(/Users/pro2.0/Desktop/ec552_project/HELIXPACE/blue_bg.png)")
        # self.setAutoFillBackground(True)
        #self.setWindowIcon(QIcon('icon.png'))


        usernameBox = QLabel(self)
        usernameBox.setText("Username")
        usernameBox.move(165,70)

        passwordBox = QLabel(self)
        passwordBox.setText("Password")
        passwordBox.move(165,120)

        # moreBox = QLabel(self)
        # moreBox.setText("Primer seqeunce")
        # moreBox.setFixedSize(200,40)
        # moreBox.move(140,170)

        self.IDinput = QLineEdit(self)
        self.IDinput.move(300,70)
        self.IDinput.setMaxLength(12)
        self.IDinput.setFixedWidth(500)
        #self.IDinput.setPlaceholderText("Choose a username that is 8-12 characters long and contains letters or numbers.")
        self.IDinput.setClearButtonEnabled(True)

        self.pwinput = QLineEdit(self)
        self.pwinput.move(300,120)
        self.pwinput.setMaxLength(12)
        self.pwinput.setFixedWidth(500)
        #self.pwinput.setPlaceholderText("Choose a password that is 8-12 characters long and contains letters and numbers.")
        self.pwinput.setClearButtonEnabled(True)

        # self.moreinfo = QTextEdit(self)
        # self.moreinfo.move(300,170)
        # self.moreinfo.setFixedWidth(500)
        # self.moreinfo.setFixedHeight(200)
        # self.moreinfo.setStyleSheet("QLineEdit"
        #                         "{"
        #                         "background : grey;"
        #                         "}")
        #self.moreinfo.textBackgroundColor("background : lightblue")
        # self.moreinfo.setClearButtonEnabled(True)
        # self.moreinfo.setDragEnabled(True)

        enter = QPushButton(self)
        enter.setText("Enter!")
        button2 = QPushButton(self)
        button2.setText("Return")
        button3 = QPushButton(self)
        button3.setText("Upload DNA sequence!")

        enter.move(800,400)
        enter.adjustSize()
        enter.clicked.connect(self.enter_clicked)
        button2.move(50,400)
        button2.adjustSize()
        button2.clicked.connect(self.homepage)
        button3.move(500,400)
        button3.adjustSize()
        button3.clicked.connect(self.Upload_DNA)

    def homepage(self):                                             # <===
        self.w = entrainWindow()
        self.w.show()
        self.hide()
        

    def Upload_DNA(self):
        path = QFileDialog.getOpenFileName(self, 'Open a file', '','All Files (*.*)')
        # if path != ('', ''):
        #     print(path[0])
        self.in_dir = os.path.dirname(path[0])
        self.in_file = os.path.basename(path[0])
        #return in_dir,in_file
    
    def enter_clicked(self):
        username = self.IDinput.text()
        password = self.pwinput.text()
        #IDPW = ID + password
        decoder = DNADecoder()
        check = decoder.decode_check(username, password,self.in_file)
        if check == True:
            #outfile_name= QInputDialog.getText(self, 'encoding...','out directory')
            # prefix= QInputDialog.getText(self, 'encoding...','prefix')
            decoder.decode(self.in_file,'./decode_output','result')
        elif check == False:
            QMessageBox.information(self,'Warning!','Your username and password do not match! Plesde try again!', QMessageBox.Yes)


        #print(f"USERNAME:{self.IDinput.text()} \nPASSWORD:{self.pwinput.text()}")

    
    def closeIt(self): 
        self.close()

    def start_decode(self):
        encoder = DNADecoder()
        encoder.decode(self.in_file,'./decode_output','result')




app = QApplication(sys.argv)

window = entrainWindow()
window.show()  
app.exec()