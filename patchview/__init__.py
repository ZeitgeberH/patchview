from . import patchview
from . import __version__
from PyQt5 import QtWidgets
import sys

print(f"initilizating PatchView (version {__version__.version}).\nuse patchview.pvGUI() to launch the GUI")
pvApp = QtWidgets.QApplication(sys.argv)

def pvGUI():
    main = patchview.MainWindow(pvApp)
    main.show()
    pvApp.exec()
    return
