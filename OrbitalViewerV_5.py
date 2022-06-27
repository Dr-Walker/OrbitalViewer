## I built this program to help students connect the visual 3d hydrogenic orbitals with the equations which govern their
## size and shape. The program needs to have mayavi, PyQt5, matplotlib, scipy, and numpy installed in order to run.
## If you have questions please don't hesitate to contact me at swalker96@bcit.ca
# I am by no means a "coder" so please forgive any formatting and convention violations.


import matplotlib.pyplot
import sys
from mayavi import mlab
import matplotlib.pyplot as plt
from matplotlib import style
import numpy as np
from scipy.special import lpmv as lpmv
from math import factorial as factorial
from PyQt5 import QtCore, QtGui, QtWidgets
from PyQt5.QtSvg import QSvgWidget
from io import BytesIO
from PyQt5.QtWidgets import QApplication



class Ui_MainWindow(object):
    def setupUi(self, MainWindow):
        MainWindow.setObjectName("MainWindow")
        MainWindow.resize(1080, 1080)
        self.centralwidget = QtWidgets.QWidget(MainWindow)
        self.centralwidget.setObjectName("centralwidget")
        self.titlelabel = QtWidgets.QLabel(self.centralwidget)
        self.titlelabel.setGeometry(QtCore.QRect(20, 10, 671, 41))
        font = QtGui.QFont()
        font.setPointSize(18)
        self.titlelabel.setFont(font)
        self.titlelabel.setObjectName("titlelabel")
        self.nlabel = QtWidgets.QLabel(self.centralwidget)
        self.nlabel.setGeometry(QtCore.QRect(6, 81, 91, 29))
        font = QtGui.QFont()
        font.setPointSize(18)
        self.nlabel.setFont(font)
        self.nlabel.setAlignment(QtCore.Qt.AlignCenter)
        self.nlabel.setObjectName("nlabel")
        self.llabel = QtWidgets.QLabel(self.centralwidget)
        self.llabel.setGeometry(QtCore.QRect(8, 116, 81, 29))
        font = QtGui.QFont()
        font.setPointSize(18)
        self.llabel.setFont(font)
        self.llabel.setAlignment(QtCore.Qt.AlignCenter)
        self.llabel.setObjectName("llabel")
        self.mlabel = QtWidgets.QLabel(self.centralwidget)
        self.mlabel.setGeometry(QtCore.QRect(13, 151, 91, 29))
        font = QtGui.QFont()
        font.setPointSize(18)
        self.mlabel.setFont(font)
        self.mlabel.setAlignment(QtCore.Qt.AlignCenter)
        self.mlabel.setObjectName("mlabel")
        self.nBox = QtWidgets.QComboBox(self.centralwidget)
        self.nBox.setGeometry(QtCore.QRect(110, 80, 60, 22))
        self.nBox.setObjectName("nBox")
        self.nBox.addItem("")
        self.nBox.addItem("")
        self.nBox.addItem("")
        self.nBox.addItem("")
        self.nBox.addItem("")
        self.nBox.addItem("")
        self.lBox = QtWidgets.QComboBox(self.centralwidget)
        self.lBox.setGeometry(QtCore.QRect(110, 120, 60, 22))
        self.lBox.setObjectName("lBox")
        self.lBox.addItem("")
        self.lBox.addItem("")
        self.lBox.addItem("")
        self.lBox.addItem("")
        self.lBox.addItem("")
        self.lBox.addItem("")
        self.mBox = QtWidgets.QComboBox(self.centralwidget)
        self.mBox.setGeometry(QtCore.QRect(110, 160, 60, 22))
        self.mBox.setObjectName("mBox")
        self.mBox.addItem("")
        self.mBox.addItem("")
        self.mBox.addItem("")
        self.mBox.addItem("")
        self.mBox.addItem("")
        self.mBox.addItem("")
        self.mBox.addItem("")
        self.mBox.addItem("")
        self.mBox.addItem("")
        self.mBox.addItem("")
        self.mBox.addItem("")
        self.RadialBut = QtWidgets.QPushButton(self.centralwidget)
        self.RadialBut.setGeometry(QtCore.QRect(45, 260, 181, 141))
        self.RadialBut.setObjectName("RadialBut")
        self.ContourBut = QtWidgets.QPushButton(self.centralwidget)
        self.ContourBut.setGeometry(QtCore.QRect(545, 260, 181, 141))
        self.ContourBut.setObjectName("ContourBut")
        self.Orbbutton = QtWidgets.QPushButton(self.centralwidget)
        self.Orbbutton.setGeometry(QtCore.QRect(295, 260, 181, 141))
        self.Orbbutton.setObjectName("Orbbutton")
        self.Clearbutton = QtWidgets.QPushButton(self.centralwidget)
        self.Clearbutton.setGeometry(QtCore.QRect(795, 525, 181, 141))
        self.Clearbutton.setObjectName("Clearbutton")
        self.frame = QtWidgets.QFrame(self.centralwidget)
        self.frame.setGeometry(QtCore.QRect(460, 80, 261, 141))
        font = QtGui.QFont()
        font.setPointSize(10)
        self.frame.setFont(font)
        self.frame.setStyleSheet("background-color: rgb(255, 61, 2);")
        self.frame.setFrameShape(QtWidgets.QFrame.StyledPanel)
        self.frame.setFrameShadow(QtWidgets.QFrame.Raised)
        self.frame.setObjectName("frame")
### Frame2 is the radial equation box
        self.frame2 = QtWidgets.QFrame(self.centralwidget)
        self.frame2.setGeometry(QtCore.QRect(45, 525, 680, 141))
        font = QtGui.QFont()
        font.setPointSize(10)
        self.frame2.setFont(font)
        self.frame2.setStyleSheet("background-color: rgb(33, 235, 235);")
        self.frame2.setFrameShape(QtWidgets.QFrame.StyledPanel)
        self.frame2.setFrameShadow(QtWidgets.QFrame.Raised)
        self.frame2.setObjectName("frame2")
### Frame 3 is the wavefunction equation box
        self.frame3 = QtWidgets.QFrame(self.centralwidget)
        self.frame3.setGeometry(QtCore.QRect(40, 725, 980, 150))
        font = QtGui.QFont()
        font.setPointSize(10)
        self.frame3.setFont(font)
        self.frame3.setStyleSheet("background-color: rgb(255, 10, 206);")
        self.frame3.setFrameShape(QtWidgets.QFrame.StyledPanel)
        self.frame3.setFrameShadow(QtWidgets.QFrame.Raised)
        self.frame3.setObjectName("frame3")
#### Add a 4th frame for sigma = Zr/ao
        self.frame4 = QtWidgets.QFrame(self.centralwidget)
        self.frame4.setGeometry(QtCore.QRect(40, 925, 261, 141))
        font = QtGui.QFont()
        font.setPointSize(10)
        self.frame4.setFont(font)
        self.frame4.setStyleSheet("background-color: rgb(0, 187, 249);")
        self.frame4.setFrameShape(QtWidgets.QFrame.StyledPanel)
        self.frame4.setFrameShadow(QtWidgets.QFrame.Raised)
        self.frame4.setObjectName("frame4")


        self.Commentlab = QtWidgets.QLabel(self.frame)
        self.Commentlab.setGeometry(QtCore.QRect(20, 20, 201, 100))
        font = QtGui.QFont()
        font.setPointSize(10)
        self.Commentlab.setFont(font)
        self.Commentlab.setAlignment(QtCore.Qt.AlignLeading|QtCore.Qt.AlignLeft|QtCore.Qt.AlignTop)
        self.Commentlab.setWordWrap(True)
        self.Commentlab.setObjectName("Commentlab")

        self.Equationlab = QtWidgets.QLabel(self.frame2)
        self.Equationlab.setGeometry(QtCore.QRect(20, 20, 601, 100))
        font = QtGui.QFont()
        font.setPointSize(10)
        self.Equationlab.setFont(font)
        self.Equationlab.setAlignment(QtCore.Qt.AlignLeading | QtCore.Qt.AlignLeft | QtCore.Qt.AlignTop)
        self.Equationlab.setWordWrap(True)
        self.Equationlab.setObjectName("Equationlab")

        self.WaveEquationlab = QtWidgets.QLabel(self.frame3)
        self.WaveEquationlab.setGeometry(QtCore.QRect(20, 20, 701, 100))
        font = QtGui.QFont()
        font.setPointSize(10)
        self.WaveEquationlab.setFont(font)
        self.WaveEquationlab.setAlignment(QtCore.Qt.AlignLeading | QtCore.Qt.AlignLeft | QtCore.Qt.AlignTop)
        self.WaveEquationlab.setWordWrap(True)
        self.WaveEquationlab.setObjectName("WaveEquationlab")

        self.Cboxlabel = QtWidgets.QLabel(self.centralwidget)
        self.Cboxlabel.setGeometry(QtCore.QRect(470, 50, 171, 21))
        self.Cboxlabel.setObjectName("Cboxlabel")
        self.EQboxlabel = QtWidgets.QLabel(self.centralwidget)
        self.EQboxlabel.setGeometry(QtCore.QRect(45, 480, 200, 21))
        self.EQboxlabel.setObjectName("EQboxlabel")
        self.WEQboxlabel = QtWidgets.QLabel(self.centralwidget)
        self.WEQboxlabel.setGeometry(QtCore.QRect(45, 680, 200, 21))
        self.WEQboxlabel.setObjectName("WEQboxlabel")
        MainWindow.setCentralWidget(self.centralwidget)
        self.menubar = QtWidgets.QMenuBar(MainWindow)
        self.menubar.setGeometry(QtCore.QRect(0, 0, 1080, 22))
        self.menubar.setObjectName("menubar")
        MainWindow.setMenuBar(self.menubar)
        self.statusbar = QtWidgets.QStatusBar(MainWindow)
        self.statusbar.setObjectName("statusbar")
        MainWindow.setStatusBar(self.statusbar)

        self.RadialBut.clicked.connect(self.Radclicked)
        self.Orbbutton.clicked.connect(self.pltorbclicked)
        self.ContourBut.clicked.connect(self.pltContourclicked)
        self.Clearbutton.clicked.connect(self.Clearclicked)

        self.retranslateUi(MainWindow)
        QtCore.QMetaObject.connectSlotsByName(MainWindow)

    def retranslateUi(self, MainWindow):
        _translate = QtCore.QCoreApplication.translate
        MainWindow.setWindowTitle(_translate("MainWindow", "MainWindow"))
        self.titlelabel.setText(_translate("MainWindow", "Which orbital are you interested in?"))
        self.nlabel.setText(_translate("MainWindow", "n = "))
        self.llabel.setText(_translate("MainWindow", "l = "))
        self.mlabel.setText(_translate("MainWindow", "m = "))
        self.nBox.setItemText(0, _translate("MainWindow", "1"))
        self.nBox.setItemText(1, _translate("MainWindow", "2"))
        self.nBox.setItemText(2, _translate("MainWindow", "3"))
        self.nBox.setItemText(3, _translate("MainWindow", "4"))
        self.nBox.setItemText(4, _translate("MainWindow", "5"))
        self.nBox.setItemText(5, _translate("MainWindow", "6"))
        self.lBox.setItemText(0, _translate("MainWindow", "0"))
        self.lBox.setItemText(1, _translate("MainWindow", "1"))
        self.lBox.setItemText(2, _translate("MainWindow", "2"))
        self.lBox.setItemText(3, _translate("MainWindow", "3"))
        self.lBox.setItemText(4, _translate("MainWindow", "4"))
        self.lBox.setItemText(5, _translate("MainWindow", "5"))
        self.mBox.setItemText(0, _translate("MainWindow", "-5"))
        self.mBox.setItemText(1, _translate("MainWindow", "-4"))
        self.mBox.setItemText(2, _translate("MainWindow", "-3"))
        self.mBox.setItemText(3, _translate("MainWindow", "-2"))
        self.mBox.setItemText(4, _translate("MainWindow", "-1"))
        self.mBox.setItemText(5, _translate("MainWindow", "0"))
        self.mBox.setItemText(6, _translate("MainWindow", "1"))
        self.mBox.setItemText(7, _translate("MainWindow", "2"))
        self.mBox.setItemText(8, _translate("MainWindow", "3"))
        self.mBox.setItemText(9, _translate("MainWindow", "4"))
        self.mBox.setItemText(10, _translate("MainWindow", "5"))
        self.RadialBut.setText(_translate("MainWindow", "Plot Radial Part"))
        self.Orbbutton.setText(_translate("MainWindow", "Plot 3D Orbital"))
        self.ContourBut.setText(_translate("MainWindow", "Plot Contours"))
        self.Clearbutton.setText(_translate("MainWindow", "Clear Equation"))
        self.frame.setStatusTip(_translate("MainWindow", "Any errors will show up here"))
        self.Commentlab.setText(_translate("MainWindow", "No Errors"))
        self.Cboxlabel.setText(_translate("MainWindow", "Comment Box"))
        self.Equationlab.setText(_translate("MainWindow", "Equations will show here"))
        self.EQboxlabel.setText(_translate("MainWindow", "Radial Equation Box"))
        self.WEQboxlabel.setText(_translate("MainWindow", "Wave Equation Box"))

    def Clearclicked(self):
        self.svg.close()
        self.svg2.close()

    def Radclicked(self):
        n = int(self.nBox.currentText())
        l = int(self.lBox.currentText())
        m = int(self.mBox.currentText())

#### Setting l values to spdf notation
        self.Equationlab.setText("")
        if l == 0:
            econfig='s'
        elif l==1:
            econfig='p'
        elif l==2:
            econfig = 'd'
        elif l==3:
            econfig = 'f'
        elif l==4:
            econfig = 'g'
        elif l==5:
            econfig = 'h'
        else:
            econfig = 'na'
###End section
###Radial EQ box label
        EQboxtext = r'Radial part of the %s%s m=%s' % (n, econfig, m)
        self.EQboxlabel.setText(EQboxtext)
        WEQboxtext = r'Wavefunction of the %s%s m=%s' % (n, econfig, m)
        self.EQboxlabel.setText(EQboxtext)
###Converting text to image for display
        plt.rc('mathtext', fontset='cm')
        #plt.rcParams['text.usetex'] = True
        def tex2svg(formula, fontsize=12, dpi=300):
            """Render TeX formula to SVG.
            Args:
                formula (str): TeX formula.
                fontsize (int, optional): Font size.
                dpi (int, optional): DPI.
            Returns:
                str: SVG render.
            """

            fig = plt.figure(figsize=(0.01, 0.01))
            fig.text(0, 0, r'{}'.format(formula), fontsize=fontsize)

            output = BytesIO()
            fig.savefig(output, dpi=dpi, transparent=True, format='svg',
                        bbox_inches='tight', pad_inches=0.0)
            plt.close(fig)

            output.seek(0)
            return output.read()
### end txt to svg
### Radial Equations
        if n == 1:
            FORMULA = r'$R(1s) =2 (\frac{Z}{\alpha _0})^\frac{3}{2}e^\frac{-\sigma}{2}$'
        elif n==2 and l==0:
            FORMULA = r'$R(2s) = \frac{1}{2\sqrt{2}} (\frac{Z}{\alpha _0})^\frac{3}{2}(2-\sigma)e^\frac{-\sigma}{2} $'
        elif n==2 and l==1:
            FORMULA =r'$R(2p) = \frac{1}{2\sqrt{6}} (\frac{Z}{\alpha _0})^\frac{3}{2}e^\frac{-\sigma}{2}$'
        elif n == 3 and l == 0:
            FORMULA = r'$R(3s) = \frac{1}{9\sqrt{3}} (\frac{Z}{\alpha _0})^\frac{3}{2}(6-6\sigma + \sigma^2)e^\frac{-\sigma}{2} $'
        elif n == 3 and l == 1:
            FORMULA =r'$R(3p) = \frac{1}{9\sqrt{6}} (\frac{Z}{\alpha _0})^\frac{3}{2}(4-\sigma)\sigma)e^\frac{-\sigma}{2}$'
        elif n == 3 and l == 2:
            FORMULA =r'$R(3d) = \frac{1}{9\sqrt{30}} (\frac{Z}{\alpha _0})^\frac{3}{2}\sigma^2e^\frac{-\sigma}{2}$'
        else:
            FORMULA = r'Equations only display up to n=3'

###Display Radial Equation
        app = QApplication(sys.argv)
        self.svg = QSvgWidget(self.centralwidget)
        self.frame2.setStyleSheet("background-color: rgb(33, 235, 235);")
        self.svg.setWindowTitle(EQboxtext)
        self.svg.load(tex2svg(FORMULA))
        self.svg.setGeometry(QtCore.QRect(45, 525, 680, 141))
        self.svg.show()
### End Radial Equations
#####################Full wavefunction equation
        if n == 1:
            WaveEq = r'$\Psi(1s) =\frac{1}{\sqrt{\pi}}(\frac{Z}{\alpha _0})^\frac{3}{2}e^\frac{-\sigma}{2}$'
        elif n == 2 and l == 0:
            WaveEq = r'$\Psi(2s) = \frac{1}{4\sqrt{2\pi}} (\frac{Z}{\alpha _0})^\frac{3}{2}(2-\sigma)e^\frac{-\sigma}{2} $'
        elif n == 2 and l == 1 and m==0:
            WaveEq = r'$\Psi(2pz) = \frac{1}{4\sqrt{2\pi}} (\frac{Z}{\alpha _0})^\frac{3}{2}\sigma e^\frac{-\sigma}{2} \cos \ \theta$'
        elif n == 2 and l == 1 and m==-1:
            WaveEq = r'$\Psi(2px) = \frac{1}{4\sqrt{2\pi}} (\frac{Z}{\alpha _0})^\frac{3}{2}\sigma e^\frac{-\sigma}{2} \sin{\ \theta} \ \cos{\ \phi}$'
        elif n == 2 and l == 1 and m == 1:
            WaveEq = r'$\Psi(2py) = \frac{1}{4\sqrt{2\pi}} (\frac{Z}{\alpha _0})^\frac{3}{2}\sigma e^\frac{-\sigma}{2} \sin{\ \theta}\  \sin{\ \phi}$'
        elif n == 3 and l == 0:
            WaveEq = r'$\Psi(3s) = \frac{1}{81\sqrt{3\pi}} (\frac{Z}{\alpha _0})^\frac{3}{2}(27 - 18\sigma +2\sigma^2)e^\frac{-\sigma}{3} $'
        elif n == 3 and l == 1 and m==0:
            WaveEq = r'$\Psi(3pz) = \frac{\sqrt{2}}{81\sqrt{\pi}} (\frac{Z}{\alpha _0})^\frac{3}{2}\sigma(6-\sigma)e^\frac{-\sigma}{3}\cos{\ \theta}$'
        elif n == 3 and l == 1 and m == -1:
            WaveEq = r'$\Psi(3px) = \frac{\sqrt{2}}{81\sqrt{\pi}} (\frac{Z}{\alpha _0})^\frac{3}{2}\sigma(6-\sigma)e^\frac{-\sigma}{3} \sin{\ \theta}\ \cos{\ \phi}$'
        elif n == 3 and l == 1 and m == 1:
            WaveEq = r'$\Psi(3py) = \frac{1}{81\sqrt{3\pi}} (\frac{Z}{\alpha _0})^\frac{3}{2}\sigma(6-\sigma)e^\frac{-\sigma}{3} \sin{\ \theta} \ \sin{\ \phi}$'

            ##### then convert into a new svg2, display in new box..Eqns in McQuarrie pg 340
        elif n == 3 and l == 2 and m == 0 :
            WaveEq = r'$Psi(3dz^2) = \frac{1}{81\sqrt{6\pi}} (\frac{Z}{\alpha _0})^\frac{3}{2}\sigma^2e^\frac{-\sigma}{3}(3\cos^2{\ \theta-1})$'
        elif n == 3 and l == 2 and m == -1 :
            WaveEq = r'$Psi(3dxz) = \frac{\sqrt{2}}{81\sqrt{\pi}} (\frac{Z}{\alpha _0})^\frac{3}{2}\sigma^2e^\frac{-\sigma}{3}\sin{\ \theta}\ cos{\ \theta}\ cos{\phi}$'
        elif n == 3 and l == 2 and m == 1 :
            WaveEq = r'$Psi(3dyz) = \frac{\sqrt{2}}{81\sqrt{\pi}} (\frac{Z}{\alpha _0})^\frac{3}{2}\sigma^2e^\frac{-\sigma}{3}\sin{\ \theta}\ cos{\ \theta}\ sin{\ \phi}$'
        elif n == 3 and l == 2 and m == -2 :
            WaveEq = r'$Psi(3dx^2-y^2) = \frac{\sqrt{1}}{81\sqrt{2\pi}} (\frac{Z}{\alpha _0})^\frac{3}{2}\sigma^2e^\frac{-\sigma}{3}\sin^2{\ \theta}\ cos{\ 2\phi}$'
        elif n == 3 and l == 2 and m == 2 :
            WaveEq = r'$Psi(3dxy) = \frac{\sqrt{1}}{81\sqrt{2\pi}} (\frac{Z}{\alpha _0})^\frac{3}{2}\sigma^2e^\frac{-\sigma}{3}\sin^2{\ \theta}\ sin\ {2\phi}$'
        else:
            WaveEq = r'Equations only display up to n=3'

        SigmaEq = r'$\sigma = \frac{Zr}{a_0}$'
### Display wave equations
        self.svg2 = QSvgWidget(self.centralwidget)
        self.frame2.setStyleSheet("background-color: rgb(33, 235, 235);")
        self.svg2.setWindowTitle(WEQboxtext)
        self.svg2.load(tex2svg(WaveEq))
        self.svg2.setGeometry(QtCore.QRect(45, 725, 900, 141))
        self.svg2.show()

        self.svg3 = QSvgWidget(self.centralwidget)
        self.frame3.setStyleSheet("background-color: rgb(33, 235, 235);")
        self.svg3.setWindowTitle(WEQboxtext)
        self.svg3.load(tex2svg(SigmaEq))
        self.svg3.setGeometry(QtCore.QRect(45, 945, 261, 95))
        self.svg3.show()
### Error if quantum numbers are bonkers
        if abs(m)> l:
            print("Please check your n,l,m values remember the quantum # rules")
            self.Commentlab.setText("Please check your n,l,m values remember the quantum # rules")
        elif l >= n:
            print("Please check your n,l,m values remember the quantum # rules")
            self.Commentlab.setText("Please check your n,l,m values remember the quantum # rules")
        else:
            self.Commentlab.setText("Processing, check behind the main window for plots.")
### Defining how to plot the square of radial part of the wavefunction
            def Plot_R(n, l):
                style.use('fivethirtyeight')
                fig = plt.figure()
                ax = fig.add_subplot(1, 1, 1)
                newrange = np.round((2.5 * n ** 2 + 1.5 * n + 1) / 5) * 5
                r = np.linspace(0, newrange, 100)

                plt.ylabel("electron density")
                plt.xlabel("r / $a_o$")
                plt.plot(r, psi_R(n, l, r, ))
                plt.tight_layout()
                plt.show()

### Defining how to calc the normalized square of the radial part of the wavefunction
            def psi_R(n, l, r):
                psi_R = r ** 2 * R(n, l, r) ** 2
                return psi_R
### Define the Associated Legendre function (Used to Calc the radial part of the wavefunction)
            def AssociatedLegendre(l, m, x):
                lpmv(m, l, x)  # Annoyingly lpmv needs m to come first in the call order.
                return lpmv(m, l, x)
### Define the Associated Laguerre
            def AssociatedLaguerre(n, m, x):
                Anm = 0
                for a in range(0, n + 1):
                    Zeta = (factorial(m + n)) / ((factorial((m + n) - (n - a))) * factorial(n - a))
                    Anm = Anm + factorial(n + m) * (Zeta) / factorial(a) * (-x) ** a
                return Anm
### Define the spherical harmonics
            # Spherical Ylm function uses the Associated Legendre function
            def SphericalYlm(l, m, theta, phi):
                SphericalYlm = (-1 ** m) * np.sqrt(
                    ((2 * l + 1) / (4 * np.pi)) * ((factorial(l - abs(m))) / (factorial(l + abs(m))))) \
                               * AssociatedLegendre(l, m, np.cos(theta)) * np.exp(1j * m * phi)
                return SphericalYlm
### Defines the angular part of the wavefunction, however its not needed for the radial part TB deleted
            def Y(l, m, theta, phi):
                if m < 0:
                    Y = np.sqrt(2) * (-1) ** m * np.imag(SphericalYlm(l, abs(m), theta, phi))
                elif m == 0:
                    Y = (-1 ** m) * np.sqrt(
                        ((2 * l + 1) / (4 * np.pi)) * ((factorial(1 - abs(m))) / (factorial(1 + abs(m))))) \
                        * AssociatedLegendre(l, m, np.cos(theta))
                else:
                    Y = np.sqrt(2) * (-1) ** m * np.real(SphericalYlm(l, abs(m), theta, phi))

                return Y
### Defining the radial wavefunction. Note what is plotted is the normalized square of this function.
            def R(n, l, r):
                a = 1
                R = (np.sqrt((2 / (a * n)) ** 3 * factorial(n - l - 1) / (2 * n * factorial(n + l))) * np.exp(
                    -r / (a * n)) * (2 * r / (a * n)) ** l \
                     / factorial(n - l - 1 + 2 * l + 1) * AssociatedLaguerre(n - l - 1, 2 * l + 1, 2 * r / (a * n)))
                return R
            Plot_R(n,l)

##################### What happens when plot Orbital is clicked.
    def pltorbclicked(self):
        n = int(self.nBox.currentText())
        l = int(self.lBox.currentText())
        m = int(self.mBox.currentText())
        #   m = int(self.mBox.currentText())
        self.Equationlab.setText("")
        if l == 0:
            econfig='s'
        elif l==1:
            econfig='p'
        elif l==2:
            econfig = 'd'
        elif l==3:
            econfig = 'f'
        elif l==4:
            econfig = 'g'
        elif l==5:
            econfig = 'h'
        else:
            econfig = 'na'

        ###Radial EQ box label
        EQboxtext = r'Radial part of the %s%s m=%s' % (n, econfig, m)
        self.EQboxlabel.setText(EQboxtext)
        WEQboxtext = r'Wavefunction of the %s%s m=%s' % (n, econfig, m)
        self.EQboxlabel.setText(EQboxtext)
        ###Converting text to image for display
        plt.rc('mathtext', fontset='cm')

    # plt.rcParams['text.usetex'] = True
        def tex2svg(formula, fontsize=12, dpi=300):
            """Render TeX formula to SVG.
            Args:
                formula (str): TeX formula.
                fontsize (int, optional): Font size.
                dpi (int, optional): DPI.
            Returns:
                str: SVG render.
            """

            fig = plt.figure(figsize=(0.01, 0.01))
            fig.text(0, 0, r'{}'.format(formula), fontsize=fontsize)

            output = BytesIO()
            fig.savefig(output, dpi=dpi, transparent=True, format='svg',
                        bbox_inches='tight', pad_inches=0.0)
            plt.close(fig)

            output.seek(0)
            return output.read()

        ### end txt to svg
        ### Radial Equations
        if n == 1:
            FORMULA = r'$R(1s) =2 (\frac{Z}{\alpha _0})^\frac{3}{2}e^\frac{-\sigma}{2}$'
        elif n == 2 and l == 0:
            FORMULA = r'$R(2s) = \frac{1}{2\sqrt{2}} (\frac{Z}{\alpha _0})^\frac{3}{2}(2-\sigma)e^\frac{-\sigma}{2} $'
        elif n == 2 and l == 1:
            FORMULA = r'$R(2p) = \frac{1}{2\sqrt{6}} (\frac{Z}{\alpha _0})^\frac{3}{2}e^\frac{-\sigma}{2}$'
        elif n == 3 and l == 0:
            FORMULA = r'$R(3s) = \frac{1}{9\sqrt{3}} (\frac{Z}{\alpha _0})^\frac{3}{2}(6-6\sigma + \sigma^2)e^\frac{-\sigma}{2} $'
        elif n == 3 and l == 1:
            FORMULA = r'$R(3p) = \frac{1}{9\sqrt{6}} (\frac{Z}{\alpha _0})^\frac{3}{2}(4-\sigma)\sigma)e^\frac{-\sigma}{2}$'
        elif n == 3 and l == 2:
            FORMULA = r'$R(3d) = \frac{1}{9\sqrt{30}} (\frac{Z}{\alpha _0})^\frac{3}{2}\sigma^2e^\frac{-\sigma}{2}$'
        else:
            FORMULA = r'Equations only display up to n=3'

        ###Display Radial Equation
        app = QApplication(sys.argv)
        self.svg = QSvgWidget(self.centralwidget)
        self.frame2.setStyleSheet("background-color: rgb(33, 235, 235);")
        self.svg.setWindowTitle(EQboxtext)
        self.svg.load(tex2svg(FORMULA))
        self.svg.setGeometry(QtCore.QRect(45, 525, 680, 141))
        self.svg.show()
        ### End Radial Equations
        #####################Full wavefunction equation
        if n == 1:
            WaveEq = r'$\Psi(1s) =\frac{1}{\sqrt{\pi}}(\frac{Z}{\alpha _0})^\frac{3}{2}e^\frac{-\sigma}{2}$'
        elif n == 2 and l == 0:
            WaveEq = r'$\Psi(2s) = \frac{1}{4\sqrt{2\pi}} (\frac{Z}{\alpha _0})^\frac{3}{2}(2-\sigma)e^\frac{-\sigma}{2} $'
        elif n == 2 and l == 1 and m == 0:
            WaveEq = r'$\Psi(2pz) = \frac{1}{4\sqrt{2\pi}} (\frac{Z}{\alpha _0})^\frac{3}{2}\sigma e^\frac{-\sigma}{2} \cos \ \theta$'
        elif n == 2 and l == 1 and m == -1:
            WaveEq = r'$\Psi(2px) = \frac{1}{4\sqrt{2\pi}} (\frac{Z}{\alpha _0})^\frac{3}{2}\sigma e^\frac{-\sigma}{2} \sin{\ \theta} \ \cos{\ \phi}$'
        elif n == 2 and l == 1 and m == 1:
            WaveEq = r'$\Psi(2py) = \frac{1}{4\sqrt{2\pi}} (\frac{Z}{\alpha _0})^\frac{3}{2}\sigma e^\frac{-\sigma}{2} \sin{\ \theta}\  \sin{\ \phi}$'
        elif n == 3 and l == 0:
            WaveEq = r'$\Psi(3s) = \frac{1}{81\sqrt{3\pi}} (\frac{Z}{\alpha _0})^\frac{3}{2}(27 - 18\sigma +2\sigma^2)e^\frac{-\sigma}{3} $'
        elif n == 3 and l == 1 and m == 0:
            WaveEq = r'$\Psi(3pz) = \frac{\sqrt{2}}{81\sqrt{\pi}} (\frac{Z}{\alpha _0})^\frac{3}{2}\sigma(6-\sigma)e^\frac{-\sigma}{3}\cos{\ \theta}$'
        elif n == 3 and l == 1 and m == -1:
            WaveEq = r'$\Psi(3px) = \frac{\sqrt{2}}{81\sqrt{\pi}} (\frac{Z}{\alpha _0})^\frac{3}{2}\sigma(6-\sigma)e^\frac{-\sigma}{3} \sin{\ \theta}\ \cos{\ \phi}$'
        elif n == 3 and l == 1 and m == 1:
            WaveEq = r'$\Psi(3py) = \frac{1}{81\sqrt{3\pi}} (\frac{Z}{\alpha _0})^\frac{3}{2}\sigma(6-\sigma)e^\frac{-\sigma}{3} \sin{\ \theta} \ \sin{\ \phi}$'

            ##### then convert into a new svg2, display in new box..Eqns in McQuarrie pg 340
        elif n == 3 and l == 2 and m == 0:
            WaveEq = r'$Psi(3dz^2) = \frac{1}{81\sqrt{6\pi}} (\frac{Z}{\alpha _0})^\frac{3}{2}\sigma^2e^\frac{-\sigma}{3}(3\cos^2{\ \theta-1})$'
        elif n == 3 and l == 2 and m == -1:
            WaveEq = r'$Psi(3dxz) = \frac{\sqrt{2}}{81\sqrt{\pi}} (\frac{Z}{\alpha _0})^\frac{3}{2}\sigma^2e^\frac{-\sigma}{3}\sin{\ \theta}\ cos{\ \theta}\ cos{\phi}$'
        elif n == 3 and l == 2 and m == 1:
            WaveEq = r'$Psi(3dyz) = \frac{\sqrt{2}}{81\sqrt{\pi}} (\frac{Z}{\alpha _0})^\frac{3}{2}\sigma^2e^\frac{-\sigma}{3}\sin{\ \theta}\ cos{\ \theta}\ sin{\ \phi}$'
        elif n == 3 and l == 2 and m == -2:
            WaveEq = r'$Psi(3dx^2-y^2) = \frac{\sqrt{1}}{81\sqrt{2\pi}} (\frac{Z}{\alpha _0})^\frac{3}{2}\sigma^2e^\frac{-\sigma}{3}\sin^2{\ \theta}\ cos{\ 2\phi}$'
        elif n == 3 and l == 2 and m == 2:
            WaveEq = r'$Psi(3dxy) = \frac{\sqrt{1}}{81\sqrt{2\pi}} (\frac{Z}{\alpha _0})^\frac{3}{2}\sigma^2e^\frac{-\sigma}{3}\sin^2{\ \theta}\ sin\ {2\phi}$'
        else:
            WaveEq = r'Equations only display up to n=3'

        SigmaEq = r'$\sigma = \frac{Zr}{a_0}$'
        ### Display wave equations
        self.svg2 = QSvgWidget(self.centralwidget)
        self.frame2.setStyleSheet("background-color: rgb(33, 235, 235);")
        self.svg2.setWindowTitle(WEQboxtext)
        self.svg2.load(tex2svg(WaveEq))
        self.svg2.setGeometry(QtCore.QRect(45, 725, 900, 141))
        self.svg2.show()

        self.svg3 = QSvgWidget(self.centralwidget)
        self.frame3.setStyleSheet("background-color: rgb(33, 235, 235);")
        self.svg3.setWindowTitle(WEQboxtext)
        self.svg3.load(tex2svg(SigmaEq))
        self.svg3.setGeometry(QtCore.QRect(45, 945, 261, 95))
        self.svg3.show()

        if abs(m) > l:

            self.Commentlab.setText("Please check your n,l,m values remember the quantum # rules")
        elif l >= n:

            self.Commentlab.setText("Please check your n,l,m values remember the quantum # rules")
        else:
            self.Commentlab.setText("Processing, check behind the main window for plots.")

            def Plot_psi(n, l, m, ):

                border = 100
                accuracy = 220  # Higher number smoother surface, but longer processign time.
                raster = np.linspace(-border, border, accuracy)
                [x, y, z] = np.meshgrid(raster, raster, raster)
                r = np.sqrt(x ** 2 + y ** 2 + z ** 2)
                theta = np.arccos(z / r)
                phi = np.arctan2(y, x)
                Wave3 = psi_2(n, l, m, r, theta, phi)
                mlab.figure(1, fgcolor=(1, 1, 1))
                # We create a scalar field with the module of Phi as the scalar
                src = mlab.pipeline.scalar_field(Wave3)

                src.image_data.point_data.add_array(np.sign(psi(n, l, m, r, theta, phi)).T.ravel())

                src.image_data.point_data.get_array(1).name = 'phase'
                # Make sure that the dataset is up to date with the different arrays:
                src.update()

                # We select the 'scalar' attribute, ie the norm of Phi
                src2 = mlab.pipeline.set_active_attribute(src,
                                                          point_scalars='scalar')

                # Cut isosurfaces of the norm
                contour = mlab.pipeline.contour(src2)
                contour.filter.contours = [0.0015, ]

                # Now we select the 'angle' attribute, ie the phase of Phi
                contour2 = mlab.pipeline.set_active_attribute(contour,
                                                              point_scalars='phase')

                # And we display the surface. The colormap is the current attribute: the phase.
                mlab.pipeline.surface(contour2, colormap='plasma', opacity=0.5)

                mlab.colorbar(title='Phase', orientation='vertical', nb_labels=3)

                mlab.show()

            # Spherical Harmonic functions
            def AssociatedLegendre(l, m, x):
                lpmv(m, l, x)  # Annoyingly lpmv needs m to come first in the call order.

                return lpmv(m, l, x)

            def AssociatedLaguerre(n, m, x):
                Anm = 0

                for a in range(0, n + 1):
                    Zeta = (factorial(m + n)) / ((factorial((m + n) - (n - a))) * factorial(n - a))
                    Anm = Anm + factorial(n + m) * (Zeta) / factorial(a) * (-x) ** a
                return Anm

            # Spherical Ylm function uses the Associated Legendre function
            def SphericalYlm(l, m, theta, phi):
                SphericalYlm = (-1 ** m) * np.sqrt(
                    ((2 * l + 1) / (4 * np.pi)) * ((factorial(l - abs(m))) / (factorial(l + abs(m))))) \
                               * AssociatedLegendre(l, m, np.cos(theta)) * np.exp(1j * m * phi)
                return SphericalYlm

            def Y(l, m, theta, phi):
                if m < 0:
                    Y = np.sqrt(2) * (-1) ** m * np.imag(SphericalYlm(l, abs(m), theta, phi))
                elif m == 0:
                    Y = (-1 ** m) * np.sqrt(
                        ((2 * l + 1) / (4 * np.pi)) * ((factorial(1 - abs(m))) / (factorial(1 + abs(m))))) \
                        * AssociatedLegendre(l, m, np.cos(theta))
                else:
                    Y = np.sqrt(2) * (-1) ** m * np.real(SphericalYlm(l, abs(m), theta, phi))

                return Y

            def R(n, l, r):
                a = 1
                R = (np.sqrt((2 / (a * n)) ** 3 * factorial(n - l - 1) / (2 * n * factorial(n + l))) * np.exp(
                    -r / (a * n)) * (2 * r / (a * n)) ** l \
                     / factorial(n - l - 1 + 2 * l + 1) * AssociatedLaguerre(n - l - 1, 2 * l + 1, 2 * r / (a * n)))
                return R

            # Wavefunction
            def psi(n, l, m, r, theta, phi):
                psi = R(n, l, r) * Y(l, m, theta, phi)
                return psi

            def psi_R(n, l, r):
                psi_R = r ** 2 * R(n, l, r) ** 2
                return psi_R

            def psi_2(n, l, m, r, theta, phi):
                psi_2 = (r * R(n, l, r) * Y(l, m, theta, phi)) ** 2
                return psi_2

            def Plot_psi2plot1(n, l, m):
                size = np.round((2.5 * n ** 2 + 1.5 * n + 1) / 5) * 5
                X2 = np.linspace(-size, size, 50)
                Y2 = np.linspace(-size, size, 50)
                [X2, Y2] = np.meshgrid(X2, Y2)
                rad = np.sqrt(X2 ** 2 + Y2 ** 2)
                theta = 90
                phi = np.arctan2(Y2, X2)
                Z = psi_2(n, l, m, rad, theta, phi)
                fig, ax = plt.subplots(1, 1)
                cp = ax.contourf(X2, Y2, Z)
                fig.colorbar(cp, label='Probability density')  # Add a colorbar to a plot
                ax.set_title('Filled Contours Plot')
                ax.set_xlabel('r / $a_o$')
                ax.set_ylabel('r / $a_o$')
                plt.show()

            def Plot_psi2plot2(n, l, m):
                size = np.round((2.5 * n ** 2 + 1.5 * n + 1) / 5) * 5
                Y2 = np.linspace(-size, size, 50)
                Z2 = np.linspace(-size, size, 50)
                [Y2, Z2] = np.meshgrid(Y2, Z2)
                rad = np.sqrt(Y2 ** 2 + Z2 ** 2)
                theta = np.arccos(Z2 / rad)
                phi = 90
                X2 = psi_2(n, l, m, rad, theta, phi)
                fig, ax = plt.subplots(1, 1)
                cp = ax.contourf(Y2, Z2, X2)
                fig.colorbar(cp, label='Probability density')  # Add a colorbar to a plot
                ax.set_title('Filled Contours Plot')
                ax.set_xlabel('r / $a_o$')
                ax.set_ylabel('r / $a_o$')
                plt.show()

            def Plot_psi2plot3(n, l, m):
                size = np.round((2.5 * n ** 2 + 1.5 * n + 1) / 5) * 5
                X2 = np.linspace(-size, size, 50)
                Y2 = np.linspace(-size, size, 50)
                [X2, Y2] = np.meshgrid(X2, Y2)
                rad = np.sqrt(X2 ** 2 + Y2 ** 2)
                theta = 90
                phi = np.arctan2(Y2, X2)
                Z = psi_2(n, l, m, rad, theta, phi)
                fig = plt.figure()
                ax = plt.axes(projection='3d')
                ax.plot_surface(X2, Y2, Z, cmap='viridis', edgecolor='none')
                ax.set_title('Surface plot')
                ax.set_zlabel('probability density')
                ax.set_xlabel('r / $a_o$')
                ax.set_ylabel('r / $a_o$')
                plt.show()

            def intensity_func(n, l, m, r, theta, phi):

                return (psi_2(n, l, m, abs(r), theta, phi))

            def Plot_psi2plot4(n, l, m):
                size = np.round((2.5 * n ** 2 + 1.5 * n + 1) / 5) * 5
                Y2 = np.linspace(-size, size, 50)
                Z2 = np.linspace(-size, size, 50)
                [Y2, Z2] = np.meshgrid(Y2, Z2)
                rad = np.sqrt(Y2 ** 2 + Z2 ** 2)
                theta = np.arccos(Z2 / rad)
                phi = 90
                X2 = psi_2(n, l, m, rad, theta, phi)
                fig = plt.figure()
                ax = plt.axes(projection='3d')
                ax.set_xlabel('r / $a_o$')
                ax.set_ylabel('r / $a_o$')
                ax.set_zlabel('probability density')
                ax.plot_surface(Y2, Z2, X2, cmap='viridis', edgecolor='none')
                ax.set_title('Surface plot')
                plt.show()

            def Ppsi(n, l, m):
                probabilitydensity = 1e-6
                size = np.round((2.5 * n ** 2 + 1.5 * n + 1) / 5) * 5
                border = 80
                accuracy = 150
                raster = np.linspace(-border, border, accuracy)
                [x, y, z] = np.meshgrid(raster, raster, raster)
                r = np.sqrt(x ** 2 + y ** 2 + z ** 2)
                theta = np.arccos(z / r)
                phi = np.arctan2(y, x)
                Wave = psi_2(n, l, m, r, theta, phi).real
                return Wave

            Plot_psi(n,l,m)
           # Plot_psi2plot1(n,l,m)
           # Plot_psi2plot2(n,l,m)
           # Plot_psi2plot3(n,l,m)
           # Plot_psi2plot4(n,l,m)

    def pltContourclicked(self):
        n = int(self.nBox.currentText())
        l = int(self.lBox.currentText())
        m = int(self.mBox.currentText())
        #   m = int(self.mBox.currentText())
        self.Equationlab.setText("")
        if l == 0:
            econfig='s'
        elif l==1:
            econfig='p'
        elif l==2:
            econfig = 'd'
        elif l==3:
            econfig = 'f'
        elif l==4:
            econfig = 'g'
        elif l==5:
            econfig = 'h'
        else:
            econfig = 'na'

        ###Radial EQ box label
        EQboxtext = r'Radial part of the %s%s m=%s' % (n, econfig, m)
        self.EQboxlabel.setText(EQboxtext)
        WEQboxtext = r'Wavefunction of the %s%s m=%s' % (n, econfig, m)
        self.EQboxlabel.setText(EQboxtext)
        ###Converting text to image for display
        plt.rc('mathtext', fontset='cm')

        # plt.rcParams['text.usetex'] = True
        def tex2svg(formula, fontsize=12, dpi=300):
            """Render TeX formula to SVG.
            Args:
                formula (str): TeX formula.
                fontsize (int, optional): Font size.
                dpi (int, optional): DPI.
            Returns:
                str: SVG render.
            """

            fig = plt.figure(figsize=(0.01, 0.01))
            fig.text(0, 0, r'{}'.format(formula), fontsize=fontsize)

            output = BytesIO()
            fig.savefig(output, dpi=dpi, transparent=True, format='svg',
                        bbox_inches='tight', pad_inches=0.0)
            plt.close(fig)

            output.seek(0)
            return output.read()

        ### end txt to svg
        ### Radial Equations
        if n == 1:
            FORMULA = r'$R(1s) =2 (\frac{Z}{\alpha _0})^\frac{3}{2}e^\frac{-\sigma}{2}$'
        elif n == 2 and l == 0:
            FORMULA = r'$R(2s) = \frac{1}{2\sqrt{2}} (\frac{Z}{\alpha _0})^\frac{3}{2}(2-\sigma)e^\frac{-\sigma}{2} $'
        elif n == 2 and l == 1:
            FORMULA = r'$R(2p) = \frac{1}{2\sqrt{6}} (\frac{Z}{\alpha _0})^\frac{3}{2}e^\frac{-\sigma}{2}$'
        elif n == 3 and l == 0:
            FORMULA = r'$R(3s) = \frac{1}{9\sqrt{3}} (\frac{Z}{\alpha _0})^\frac{3}{2}(6-6\sigma + \sigma^2)e^\frac{-\sigma}{2} $'
        elif n == 3 and l == 1:
            FORMULA = r'$R(3p) = \frac{1}{9\sqrt{6}} (\frac{Z}{\alpha _0})^\frac{3}{2}(4-\sigma)\sigma)e^\frac{-\sigma}{2}$'
        elif n == 3 and l == 2:
            FORMULA = r'$R(3d) = \frac{1}{9\sqrt{30}} (\frac{Z}{\alpha _0})^\frac{3}{2}\sigma^2e^\frac{-\sigma}{2}$'
        else:
            FORMULA = r'Equations only display up to n=3'

        ###Display Radial Equation
        app = QApplication(sys.argv)
        self.svg = QSvgWidget(self.centralwidget)
        self.frame2.setStyleSheet("background-color: rgb(33, 235, 235);")
        self.svg.setWindowTitle(EQboxtext)
        self.svg.load(tex2svg(FORMULA))
        self.svg.setGeometry(QtCore.QRect(45, 525, 680, 141))
        self.svg.show()
        ### End Radial Equations
        #####################Full wavefunction equation
        if n == 1:
            WaveEq = r'$\Psi(1s) =\frac{1}{\sqrt{\pi}}(\frac{Z}{\alpha _0})^\frac{3}{2}e^\frac{-\sigma}{2}$'
        elif n == 2 and l == 0:
            WaveEq = r'$\Psi(2s) = \frac{1}{4\sqrt{2\pi}} (\frac{Z}{\alpha _0})^\frac{3}{2}(2-\sigma)e^\frac{-\sigma}{2} $'
        elif n == 2 and l == 1 and m == 0:
            WaveEq = r'$\Psi(2pz) = \frac{1}{4\sqrt{2\pi}} (\frac{Z}{\alpha _0})^\frac{3}{2}\sigma e^\frac{-\sigma}{2} \cos \ \theta$'
        elif n == 2 and l == 1 and m == -1:
            WaveEq = r'$\Psi(2px) = \frac{1}{4\sqrt{2\pi}} (\frac{Z}{\alpha _0})^\frac{3}{2}\sigma e^\frac{-\sigma}{2} \sin{\ \theta} \ \cos{\ \phi}$'
        elif n == 2 and l == 1 and m == 1:
            WaveEq = r'$\Psi(2py) = \frac{1}{4\sqrt{2\pi}} (\frac{Z}{\alpha _0})^\frac{3}{2}\sigma e^\frac{-\sigma}{2} \sin{\ \theta}\  \sin{\ \phi}$'
        elif n == 3 and l == 0:
            WaveEq = r'$\Psi(3s) = \frac{1}{81\sqrt{3\pi}} (\frac{Z}{\alpha _0})^\frac{3}{2}(27 - 18\sigma +2\sigma^2)e^\frac{-\sigma}{3} $'
        elif n == 3 and l == 1 and m == 0:
            WaveEq = r'$\Psi(3pz) = \frac{\sqrt{2}}{81\sqrt{\pi}} (\frac{Z}{\alpha _0})^\frac{3}{2}\sigma(6-\sigma)e^\frac{-\sigma}{3}\cos{\ \theta}$'
        elif n == 3 and l == 1 and m == -1:
            WaveEq = r'$\Psi(3px) = \frac{\sqrt{2}}{81\sqrt{\pi}} (\frac{Z}{\alpha _0})^\frac{3}{2}\sigma(6-\sigma)e^\frac{-\sigma}{3} \sin{\ \theta}\ \cos{\ \phi}$'
        elif n == 3 and l == 1 and m == 1:
            WaveEq = r'$\Psi(3py) = \frac{1}{81\sqrt{3\pi}} (\frac{Z}{\alpha _0})^\frac{3}{2}\sigma(6-\sigma)e^\frac{-\sigma}{3} \sin{\ \theta} \ \sin{\ \phi}$'

            ##### then convert into a new svg2, display in new box..Eqns in McQuarrie pg 340
        elif n == 3 and l == 2 and m == 0:
            WaveEq = r'$Psi(3dz^2) = \frac{1}{81\sqrt{6\pi}} (\frac{Z}{\alpha _0})^\frac{3}{2}\sigma^2e^\frac{-\sigma}{3}(3\cos^2{\ \theta-1})$'
        elif n == 3 and l == 2 and m == -1:
            WaveEq = r'$Psi(3dxz) = \frac{\sqrt{2}}{81\sqrt{\pi}} (\frac{Z}{\alpha _0})^\frac{3}{2}\sigma^2e^\frac{-\sigma}{3}\sin{\ \theta}\ cos{\ \theta}\ cos{\phi}$'
        elif n == 3 and l == 2 and m == 1:
            WaveEq = r'$Psi(3dyz) = \frac{\sqrt{2}}{81\sqrt{\pi}} (\frac{Z}{\alpha _0})^\frac{3}{2}\sigma^2e^\frac{-\sigma}{3}\sin{\ \theta}\ cos{\ \theta}\ sin{\ \phi}$'
        elif n == 3 and l == 2 and m == -2:
            WaveEq = r'$Psi(3dx^2-y^2) = \frac{\sqrt{1}}{81\sqrt{2\pi}} (\frac{Z}{\alpha _0})^\frac{3}{2}\sigma^2e^\frac{-\sigma}{3}\sin^2{\ \theta}\ cos{\ 2\phi}$'
        elif n == 3 and l == 2 and m == 2:
            WaveEq = r'$Psi(3dxy) = \frac{\sqrt{1}}{81\sqrt{2\pi}} (\frac{Z}{\alpha _0})^\frac{3}{2}\sigma^2e^\frac{-\sigma}{3}\sin^2{\ \theta}\ sin\ {2\phi}$'
        else:
            WaveEq = r'Equations only display up to n=3'

        SigmaEq = r'$\sigma = \frac{Zr}{a_0}$'
        ### Display wave equations
        self.svg2 = QSvgWidget(self.centralwidget)
        self.frame2.setStyleSheet("background-color: rgb(33, 235, 235);")
        self.svg2.setWindowTitle(WEQboxtext)
        self.svg2.load(tex2svg(WaveEq))
        self.svg2.setGeometry(QtCore.QRect(45, 725, 900, 141))
        self.svg2.show()

        self.svg3 = QSvgWidget(self.centralwidget)
        self.frame3.setStyleSheet("background-color: rgb(33, 235, 235);")
        self.svg3.setWindowTitle(WEQboxtext)
        self.svg3.load(tex2svg(SigmaEq))
        self.svg3.setGeometry(QtCore.QRect(45, 945, 261, 95))
        self.svg3.show()
        if abs(m) > l:

            self.Commentlab.setText("Please check your n,l,m values remember the quantum # rules")
        elif l >= n:

            self.Commentlab.setText("Please check your n,l,m values remember the quantum # rules")
        else:
            self.Commentlab.setText("Processing, check behind the main window for plots.")
            # Spherical Harmonic functions
            def AssociatedLegendre(l, m, x):
                lpmv(m, l, x)  # Annoyingly lpmv needs m to come first in the call order.

                return lpmv(m, l, x)

            def AssociatedLaguerre(n, m, x):
                Anm = 0

                for a in range(0, n + 1):
                    Zeta = (factorial(m + n)) / ((factorial((m + n) - (n - a))) * factorial(n - a))
                    Anm = Anm + factorial(n + m) * (Zeta) / factorial(a) * (-x) ** a
                return Anm

            # Spherical Ylm function uses the Associated Legendre function
            def SphericalYlm(l, m, theta, phi):
                SphericalYlm = (-1 ** m) * np.sqrt(
                    ((2 * l + 1) / (4 * np.pi)) * ((factorial(l - abs(m))) / (factorial(l + abs(m))))) \
                               * AssociatedLegendre(l, m, np.cos(theta)) * np.exp(1j * m * phi)
                return SphericalYlm

            def Y(l, m, theta, phi):
                if m < 0:
                    Y = np.sqrt(2) * (-1) ** m * np.imag(SphericalYlm(l, abs(m), theta, phi))
                elif m == 0:
                    Y = (-1 ** m) * np.sqrt(
                        ((2 * l + 1) / (4 * np.pi)) * ((factorial(1 - abs(m))) / (factorial(1 + abs(m))))) \
                        * AssociatedLegendre(l, m, np.cos(theta))
                else:
                    Y = np.sqrt(2) * (-1) ** m * np.real(SphericalYlm(l, abs(m), theta, phi))

                return Y

            def R(n, l, r):
                a = 1
                R = (np.sqrt((2 / (a * n)) ** 3 * factorial(n - l - 1) / (2 * n * factorial(n + l))) * np.exp(
                    -r / (a * n)) * (2 * r / (a * n)) ** l \
                     / factorial(n - l - 1 + 2 * l + 1) * AssociatedLaguerre(n - l - 1, 2 * l + 1, 2 * r / (a * n)))
                return R

            # Wavefunction
            def psi(n, l, m, r, theta, phi):
                psi = R(n, l, r) * Y(l, m, theta, phi)
                return psi

            def psi_R(n, l, r):
                psi_R = r ** 2 * R(n, l, r) ** 2
                return psi_R

            def psi_2(n, l, m, r, theta, phi):
                psi_2 = (r * R(n, l, r) * Y(l, m, theta, phi)) ** 2
                return psi_2

            def Plot_psi2plot1(n, l, m):
                size = np.round((2.5 * n ** 2 + 1.5 * n + 1) / 5) * 5
                X2 = np.linspace(-size, size, 50)
                Y2 = np.linspace(-size, size, 50)
                [X2, Y2] = np.meshgrid(X2, Y2)
                rad = np.sqrt(X2 ** 2 + Y2 ** 2)
                theta = 90
                phi = np.arctan2(Y2, X2)
                Z = psi_2(n, l, m, rad, theta, phi)
                fig, ax = plt.subplots(1, 1)
                cp = ax.contourf(X2, Y2, Z)
                fig.colorbar(cp, label='Probability density')  # Add a colorbar to a plot
                ax.set_title('Filled Contours Plot')
                ax.set_xlabel('x-axis r / $a_o$')
                ax.set_ylabel('y-axis r / $a_o$')
                plt.show()

            def Plot_psi2plot2(n, l, m):
                size = np.round((2.5 * n ** 2 + 1.5 * n + 1) / 5) * 5
                Y2 = np.linspace(-size, size, 50)
                Z2 = np.linspace(-size, size, 50)
                [Y2, Z2] = np.meshgrid(Y2, Z2)
                rad = np.sqrt(Y2 ** 2 + Z2 ** 2)
                theta = np.arccos(Z2 / rad)
                phi = 90
                X2 = psi_2(n, l, m, rad, theta, phi)
                fig, ax = plt.subplots(1, 1)
                cp = ax.contourf(Y2, Z2, X2)
                fig.colorbar(cp, label='Probability density')  # Add a colorbar to a plot
                ax.set_title('Filled Contours Plot')
                ax.set_xlabel('x-axis r / $a_o$')
                ax.set_ylabel('z-axis r / $a_o$')
                plt.show()

            def Plot_psi2plot3(n, l, m):
                size = np.round((2.5 * n ** 2 + 1.5 * n + 1) / 5) * 5
                X2 = np.linspace(-size, size, 50)
                Y2 = np.linspace(-size, size, 50)
                [X2, Y2] = np.meshgrid(X2, Y2)
                rad = np.sqrt(X2 ** 2 + Y2 ** 2)
                theta = 90
                phi = np.arctan2(Y2, X2)
                Z = psi_2(n, l, m, rad, theta, phi)
                fig = plt.figure()
                ax = plt.axes(projection='3d')
                ax.plot_surface(X2, Y2, Z, cmap='viridis', edgecolor='none')
                ax.set_title('Surface plot')
                ax.set_zlabel('probability density')
                ax.set_xlabel('x-axis r / $a_o$')
                ax.set_ylabel('y-axis r / $a_o$')
                plt.show()

            def intensity_func(n, l, m, r, theta, phi):

                return (psi_2(n, l, m, abs(r), theta, phi))

            def Plot_psi2plot4(n, l, m):
                size = np.round((2.5 * n ** 2 + 1.5 * n + 1) / 5) * 5
                Y2 = np.linspace(-size, size, 50)
                Z2 = np.linspace(-size, size, 50)
                [Y2, Z2] = np.meshgrid(Y2, Z2)
                rad = np.sqrt(Y2 ** 2 + Z2 ** 2)
                theta = np.arccos(Z2 / rad)
                phi = 90
                X2 = psi_2(n, l, m, rad, theta, phi)
                fig = plt.figure()
                ax = plt.axes(projection='3d')
                ax.set_xlabel('x-axis r / $a_o$')
                ax.set_ylabel('z-axis r / $a_o$')
                ax.set_zlabel('probability density')
                ax.plot_surface(Y2, Z2, X2, cmap='viridis', edgecolor='none')
                ax.set_title('Surface plot')
                plt.show()

            def Ppsi(n, l, m):
                probabilitydensity = 1e-6
                size = np.round((2.5 * n ** 2 + 1.5 * n + 1) / 5) * 5
                border = 80
                accuracy = 150
                raster = np.linspace(-border, border, accuracy)
                [x, y, z] = np.meshgrid(raster, raster, raster)
                r = np.sqrt(x ** 2 + y ** 2 + z ** 2)
                theta = np.arccos(z / r)
                phi = np.arctan2(y, x)
                Wave = psi_2(n, l, m, r, theta, phi).real
                return Wave

           # Plot_psi(n,l,m)
            Plot_psi2plot1(n,l,m)
            Plot_psi2plot2(n,l,m)
            Plot_psi2plot3(n,l,m)
            Plot_psi2plot4(n,l,m)

if __name__ == "__main__":
    import sys
    app = QtWidgets.QApplication(sys.argv)
    MainWindow = QtWidgets.QMainWindow()
    ui = Ui_MainWindow()
    ui.setupUi(MainWindow)
    MainWindow.show()
    sys.exit(app.exec_())
    time.sleep(14)
