from . import patchview
from PyQt5 import QtWidgets
import sys

__version__ = "0.2.5"
print("initilizating PatchView. use patchview.pvGUI() to launch the GUI")
pvApp = QtWidgets.QApplication(sys.argv)

def pvGUI():
    main = patchview.MainWindow(pvApp)
    main.show()
    pvApp.exec()
    return
