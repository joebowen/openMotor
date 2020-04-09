import xml.etree.ElementTree as ET

import motorlib

from ..converter import Importer

# BS type -> oM class for all grains we can import
SUPPORTED_GRAINS = {
    '1': motorlib.grains.BatesGrain,
    '2': motorlib.grains.DGrain,
    '3': motorlib.grains.MoonBurner,
    '5': motorlib.grains.CGrain,
    '6': motorlib.grains.XCore,
    '7': motorlib.grains.Finocyl
}

# BS type -> label for grains we know about but can't import
UNSUPPORTED_GRAINS = {
    '4': 'Star',
    '8': 'Tablet',
    '9': 'Pie Segment'
}

def inToM(value):
    """Converts a string containing a value in inches to a float of meters"""
    return motorlib.units.convert(float(value), 'in', 'm')

def importPropellant(node):
    propellant = motorlib.hybrid_propellant.HybridPropellant()
    propTab = motorlib.hybrid_propellant.HybridPropellantTab()
    propellant.setProperty('name', node.attrib['Name'])
    ballN = float(node.attrib['BallisticN'])
    ballA = float(node.attrib['BallisticA']) * 1/(6895**ballN)
    propTab.setProperty('n', ballN)
    # Conversion only does in/s to m/s, the rest is handled above
    ballA = motorlib.units.convert(ballA, 'in/(s*psi^n)', 'm/(s*Pa^n)')
    propTab.setProperty('a', ballA)
    density = motorlib.units.convert(float(node.attrib['Density']), 'lb/in^3', 'kg/m^3')
    propellant.setProperty('density', density)
    propTab.setProperty('k', float(node.attrib['SpecificHeatRatio']))
    impMolarMass = node.attrib['MolarMass']
    # If the user has entered 0, override it to match the default propellant.
    if impMolarMass == '0':
        propTab.setProperty('m', 23.67)
    else:
        propTab.setProperty('m', float(impMolarMass))
    # Burnsim doesn't provide this property. Set it to match the default propellant.
    propTab.setProperty('t', 3500)
    propTab.setProperty('minPressure', 0)
    propTab.setProperty('maxPressure', 6.895e+06)
    propellant.setProperty('tabs', [propTab.getProperties()])
    return propellant

class BurnSimImporter(Importer):
    def __init__(self, manager):
        super().__init__(manager, 'BurnSim Motor', 'Loads motor files for BurnSim 3.0', {'.bsx': 'BurnSim Files'})

    def doConversion(self, path):
        motor = motorlib.hybrid_motor.HybridMotor()
        motor.config.setProperties(self.manager.preferences.general.getProperties())
        tree = ET.parse(path)
        root = tree.getroot()
        errors = ''
        propSet = False
        if root.find('Grain') is None:
            errors += "Motor contains no grains, but the propellant and nozzle can still be imported.\n"
        for child in root:
            if child.tag == 'Nozzle':
                motor.nozzle.setProperty('throat', inToM(child.attrib['ThroatDia']))
                motor.nozzle.setProperty('exit', inToM(child.attrib['ExitDia']))
                motor.nozzle.setProperty('efficiency', float(child.attrib['NozzleEfficiency']) / 100)
                motor.nozzle.setProperty('divAngle', 15)
                motor.nozzle.setProperty('convAngle', 45)
                errors += 'Nozzle angles not specified, assumed to be 15° and 45°.\n'
            if child.tag == 'Grain':
                if child.attrib['Type'] in SUPPORTED_GRAINS:
                    motor.grains.append(SUPPORTED_GRAINS[child.attrib['Type']]())
                    motor.grains[-1].setProperty('diameter', inToM(child.attrib['Diameter']))
                    motor.grains[-1].setProperty('length', inToM(child.attrib['Length']))

                    grainType = child.attrib['Type']

                    if child.attrib['EndsInhibited'] == '1':
                        motor.grains[-1].setProperty('inhibitedEnds', 'Top')
                    elif child.attrib['EndsInhibited'] == '2':
                        motor.grains[-1].setProperty('inhibitedEnds', 'Both')

                    if grainType in ('1', '3', '7'): # Grains with core diameter
                        motor.grains[-1].setProperty('coreDiameter', inToM(child.attrib['CoreDiameter']))

                    if grainType == '2': # D grain specific properties
                        motor.grains[-1].setProperty('slotOffset', inToM(child.attrib['EdgeOffset']))

                    elif grainType == '3': # Moonburner specific properties
                        motor.grains[-1].setProperty('coreOffset', inToM(child.attrib['CoreOffset']))

                    elif grainType == '5': # C grain specific properties
                        motor.grains[-1].setProperty('slotWidth', inToM(child.attrib['SlotWidth']))
                        radius = motor.grains[-1].getProperty('diameter') / 2
                        motor.grains[-1].setProperty('slotOffset', radius - inToM(child.attrib['SlotDepth']))

                    elif grainType == '6': # X core specific properties
                        motor.grains[-1].setProperty('slotWidth', inToM(child.attrib['SlotWidth']))
                        motor.grains[-1].setProperty('slotLength', inToM(child.attrib['CoreDiameter']) / 2)

                    elif grainType == '7': # Finocyl specific properties
                        motor.grains[-1].setProperty('finWidth', inToM(child.attrib['FinWidth']))
                        motor.grains[-1].setProperty('finLength', inToM(child.attrib['FinLength']))
                        motor.grains[-1].setProperty('numFins', int(child.attrib['FinCount']))

                    if not propSet: # Use propellant numbers from the forward grain
                        motor.propellant = importPropellant(child.find('Propellant'))
                        propSet = True

                else:
                    if child.attrib['Type'] in UNSUPPORTED_GRAINS:
                        errors += "File contains a "
                        errors += UNSUPPORTED_GRAINS[child.attrib['Type']]
                        errors += " grain, which can't be imported.\n"
                    else:
                        errors += "File contains an unknown grain of type " + child.attrib['Type'] + '.\n'

            if child.tag == 'TestData':
                errors += "\nFile contains test data, which is not imported."

            if child.tag == "Propellant":
                if not propSet:
                    motor.propellant = importPropellant(child)
                    propSet = True

        if errors != '':
            self.manager.app.outputMessage(errors + '\nThe rest of the motor will be imported.')

        self.manager.startFromMotor(motor)
