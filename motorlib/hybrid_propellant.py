"""HybridPropellant submodule that contains the hybrid fuel class."""

from .properties import PropertyCollection, FloatProperty, StringProperty, TabularProperty
from .simResult import SimAlert, SimAlertLevel, SimAlertType

from thermo.mixture import Mixture


class HybridPropellantTab(PropertyCollection):
    """Contains the combustion properties of a hybrid fuel."""
    def __init__(self, tabDict=None):
        super().__init__()
        self.props['minPressure'] = FloatProperty('Minimum Pressure', 'Pa', 0, 7e7)
        self.props['maxPressure'] = FloatProperty('Maximum Pressure', 'Pa', 0, 7e7)
        if tabDict is not None:
            self.setProperties(tabDict)


class HybridPropellant(PropertyCollection):
    """Contains the physical and thermodynamic properties of a hybrid fuel."""
    def __init__(self, propDict=None):
        super().__init__()
        self.props['name'] = StringProperty('Name')
        self.props['density'] = FloatProperty('Density', 'kg/m^3', 0, 10000)
        self.props['idealOxiFuelRatio'] = FloatProperty('Ratio', '', 0, 100)
        self.props['flameTemp'] = FloatProperty('Adiabatic Flame Temp', 'K', 0, 10000)
        self.props['tabs'] = TabularProperty('Properties', HybridPropellantTab)
        if propDict is not None:
            self.setProperties(propDict)

    def getCStar(self, pressure):
        """Returns the propellant's characteristic velocity."""
        temp = self.getCombustionTemp()
        gamma = self.getCombustionGamma(pressure)
        molarMass = self.getMolarMass(pressure)

        gasConst = 8314
        num = (gamma * gasConst / molarMass * temp)**0.5
        denom = gamma * ((2 / (gamma + 1))**((gamma + 1) / (gamma - 1)))**0.5
        return num / denom

    def getMolarMass(self, pressure):
        mixture = Mixture(['n2', 'h2o', 'co2'], zs=[.599, .224, .178], T=self.getProperty('flameTemp'), P=pressure)
        return mixture.MW

    def getCombustionTemp(self):
        return self.getProperty('flameTemp')

    def getGasConstant(self, pressure):
        mixture = Mixture(['n2', 'h2o', 'co2'], zs=[.599, .224, .178], T=self.getProperty('flameTemp'), P=pressure)
        return mixture.R_specific

    def getCombustionGamma(self, pressure):
        mixture = Mixture(['n2', 'h2o', 'co2'], zs=[.599, .224, .178], T=self.getProperty('flameTemp'), P=pressure)
        return mixture.Cpg / mixture.Cvg

    def getMinimumValidPressure(self):
        """Returns the lowest pressure value with associated combustion properties"""
        return min([tab['minPressure'] for tab in self.getProperty('tabs')])

    def getMaximumValidPressure(self):
        """Returns the highest pressure value with associated combustion properties"""
        return max([tab['maxPressure'] for tab in self.getProperty('tabs')])

    def getErrors(self):
        """Checks that all tabs have smaller start pressures than their end pressures, and verifies that no ranges
        overlap."""
        errors = []
        for tabId, tab in enumerate(self.getProperty('tabs')):
            if tab['maxPressure'] == tab['minPressure']:
                errText = 'Tab #' + str(tabId + 1) + ' has the same minimum and maximum pressures.'
                errors.append(SimAlert(SimAlertLevel.ERROR, SimAlertType.VALUE, errText, 'Propellant'))
            if tab['maxPressure'] < tab['minPressure']:
                errText = 'Tab #' + str(tabId + 1) + ' has reversed pressure limits.'
                errors.append(SimAlert(SimAlertLevel.ERROR, SimAlertType.VALUE, errText, 'Propellant'))
            for otherTabId, otherTab in enumerate(self.getProperty('tabs')):
                if tabId != otherTabId:
                    if otherTab['minPressure'] < tab['maxPressure'] < otherTab['maxPressure']:
                        err = 'Tabs #' + str(tabId + 1) + ' and #' + str(otherTabId + 1) + ' have overlapping ranges.'
                        errors.append(SimAlert(SimAlertLevel.ERROR, SimAlertType.VALUE, err, 'Propellant'))
        return errors

    def getPressureErrors(self, pressure):
        """Returns if the propellant has any errors associated with the supplied pressure such as not having set
        combustion properties"""
        errors = []
        for tab in self.getProperty('tabs'):
            if tab['minPressure'] < pressure < tab['maxPressure']:
                return errors
        aText = "Chamber pressure deviated from propellant's entered ranges. Results may not be accurate."
        errors.append(SimAlert(SimAlertLevel.WARNING, SimAlertType.VALUE, aText, 'Propellant'))
        return errors

    def addTab(self, tab):
        """Adds a set of combustion properties to the propellant"""
        self.props['tabs'].addTab(tab)
