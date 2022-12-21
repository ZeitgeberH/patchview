from . import patchview
from . import __version__
from PyQt5 import QtWidgets
import sys

<<<<<<< HEAD
print(f"initilizating PatchView (version {__version__.version}).\nuse patchview.pvGUI() to launch the GUI")
=======
__version__ = "0.2.5"
print("initilizating PatchView. use patchview.pvGUI() to launch the GUI")
>>>>>>> main
pvApp = QtWidgets.QApplication(sys.argv)

def pvGUI():
    main = patchview.MainWindow(pvApp)
    main.show()
    pvApp.exec()
    return
