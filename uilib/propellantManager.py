from PyQt5.QtCore import QObject, pyqtSignal

import motorlib

from .defaults import defaultPropellants
from .fileIO import loadFile, saveFile, fileTypes, getConfigPath
from .widgets.propellantMenu import PropellantMenu

class PropellantManager(QObject):

    updated = pyqtSignal()

    def __init__(self):
        super().__init__()
        self.propellants = []
        self.loadPropellants()

        self.propMenu = PropellantMenu(self)
        self.propMenu.closed.connect(self.updated.emit)

    def loadPropellants(self):
        try:
            propList = loadFile(getConfigPath() + 'propellants.yaml', fileTypes.PROPELLANTS)
            for propDict in propList:
                newProp = motorlib.hybrid_propellant.HybridPropellant()
                newProp.setProperties(propDict)
                self.propellants.append(newProp)
        except FileNotFoundError:
            self.propellants = [motorlib.hybrid_propellant.HybridPropellant(prop) for prop in defaultPropellants()]
            self.savePropellants()

    def savePropellants(self):
        propellants = [prop.getProperties() for prop in self.propellants]
        try:
            saveFile(getConfigPath() + 'propellants.yaml', propellants, fileTypes.PROPELLANTS)
        except:
            print('Unable to save propellants!')

    def getNames(self):
        return [prop.getProperty('name') for prop in self.propellants]

    def getPropellantByName(self, name):
        return self.propellants[self.getNames().index(name)]

    def showMenu(self):
        self.propMenu.setupPropList()
        self.propMenu.show()

    def setPreferences(self, pref):
        self.propMenu.ui.propEditor.setPreferences(pref)
