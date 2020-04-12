"""Contains the hybrid motor class and a supporting configuration property collection."""
from .grains import grainTypes
from .nozzle import Nozzle
from .nitrous_oxide import N2OTank
from .hybrid_propellant import HybridPropellant
from . import geometry
from . import units
from .simResult import SimulationResult, SimAlert, SimAlertLevel, SimAlertType
from .grains import EndBurningGrain
from .properties import PropertyCollection, FloatProperty, IntProperty

import math


class HybridMotorConfig(PropertyCollection):
    """Contains the settings required for simulation, including environmental conditions and details about
    how to run the simulation."""
    def __init__(self):
        super().__init__()
        # N2O Tank
        self.props['tankVolume'] = FloatProperty('N2O Tank Volume', 'm^3', 0, 100)
        self.props['tankPressure'] = FloatProperty('N2O Tank Pressure', 'Pa', 0, 7e7)
        self.props['tankUllage'] = FloatProperty('N2O Tank Ullage Percentage', '%', 0, 100)
        self.props['numOfInjectors'] = FloatProperty('Number of Injectors', '', 1, 10)
        self.props['injectorDiameter'] = FloatProperty('Injector Diameter', 'mm', .01, 5)
        self.props['injectorK2Coefficient'] = FloatProperty('Injector K2 Coefficient', '', 0, 4)
        # Limits
        self.props['maxPressure'] = FloatProperty('Maximum Allowed Pressure', 'Pa', 0, 7e7)
        self.props['maxMassFlux'] = FloatProperty('Maximum Allowed Mass Flux', 'kg/(m^2*s)', 0, 1e4)
        self.props['minPortThroat'] = FloatProperty('Minimum Allowed Port/Throat Ratio', '', 1, 4)
        # Simulation
        self.props['burnoutWebThres'] = FloatProperty('Web Burnout Threshold', 'm', 2.54e-5, 3.175e-3)
        self.props['burnoutThrustThres'] = FloatProperty('Thrust Burnout Threshold', '%', 0.01, 10)
        self.props['timestep'] = FloatProperty('Simulation Timestep', 's', 0.0001, 0.1)
        self.props['ambPressure'] = FloatProperty('Ambient Pressure', 'Pa', 0.0001, 102000)
        self.props['ambDensity'] = FloatProperty('Ambient Density', 'kg/m^3', 0.0001, 102000)
        self.props['igniterPressure'] = FloatProperty('Igniter Pressure', 'Pa', 0, 1e7)
        self.props['mapDim'] = IntProperty('Grain Map Dimension', '', 250, 2000)
        # Motor Design
        self.props['postCombustionVol'] = FloatProperty('Post Combustion Chamber Volume', 'm^3', 0, 1)
        self.props['initialCombustionVol'] = FloatProperty('Initial Combustion Chamber Volume', 'm^3', 0, 100)


class HybridMotor():
    """The motor class stores a number of grains, a nozzle instance, a propellant, and a configuration that it uses
    to run simulations. Simulations return a simRes object that includes any warnings or errors associated with the
    simulation and the data. The propellant field may be None if the motor has no propellant set, and the grains list
    is allowed to be empty. The nozzle field and config must be filled, even if they are empty property collections."""
    def __init__(self, propDict=None):
        self.grains = []
        self.propellant = None
        self.nozzle = Nozzle()
        self.config = HybridMotorConfig()

        if propDict is not None:
            self.applyDict(propDict)

        self.massInChamber = self.config.getProperty("initialCombustionVol") * self.config.getProperty("ambDensity")

        self.tank = N2OTank(
            delta_time=self.config.getProperty("timestep"),
            initial_tank_pressure=units.convert(self.config.getProperty("tankPressure"), 'Pa', 'bar'),
            initial_ullage_percentage=self.config.getProperty("tankUllage") / 100,
            tank_volume=self.config.getProperty("tankVolume"),
            num_of_injectors=self.config.getProperty("numOfInjectors"),
            injector_diameter=units.convert(self.config.getProperty("injectorDiameter"), 'mm', 'm'),
            injector_k2_coefficient=self.config.getProperty("injectorK2Coefficient")
        )

    def getDict(self):
        """Returns a serializable representation of the motor. The dictionary has keys 'nozzle', 'propellant',
        'grains', and 'config', which hold to the properties of their corresponding fields. Grains is a list
        of dicts, each containing a type and properties. Propellant may be None if the motor has no propellant
        set."""
        motorData = dict()
        motorData['nozzle'] = self.nozzle.getProperties()
        if self.propellant is not None:
            motorData['propellant'] = self.propellant.getProperties()
        else:
            motorData['propellant'] = None
        motorData['grains'] = [{'type': grain.geomName, 'properties': grain.getProperties()} for grain in self.grains]
        motorData['config'] = self.config.getProperties()
        return motorData

    def applyDict(self, dictionary):
        """Makes the motor copy properties from the dictionary that is passed in, which must be formatted like
        the result passed out by 'getDict'"""
        self.nozzle.setProperties(dictionary['nozzle'])
        if dictionary['propellant'] is not None:
            self.propellant = HybridPropellant(dictionary['propellant'])
        else:
            self.propellant = None
        self.grains = []
        for entry in dictionary['grains']:
            self.grains.append(grainTypes[entry['type']]())
            self.grains[-1].setProperties(entry['properties'])
        self.config.setProperties(dictionary['config'])

    def calcMassOutNozzle(self, dThroat, lastPressure):
        # TODO: I'm not sure where I found this.  I'm not sure this is the right way to calculate this.
        #  See solve_for_average_total_mass_flow_rate() in hybrid_rocket_motor_sim

        nozz = self.nozzle.getThroatArea(dThroat)

        temp = self.propellant.getCombustionTemp()
        gasConstant = self.propellant.getGasConstant(lastPressure)
        K = self.propellant.getCombustionGamma(lastPressure)
        molarMass = self.propellant.getMolarMass(lastPressure)

        print(f'temp: {temp}')
        print(f'gasConstant: {gasConstant}')
        print(f'K: {K}')
        print(f'molarMass: {molarMass}')

        massOutNozzle = lastPressure * nozz * math.sqrt(K / ((gasConstant / molarMass) * temp)) * ((2 / (K + 1)) ** ((K + 1) / (2 * (K - 1))))

        return massOutNozzle

    def calcIdealPressure(self, regDepth, dThroat, lastPressure, addedMass):
        if addedMass:   # There won't be any massOutNozzle if the N2O isn't flowing
            massOutNozzle = self.calcMassOutNozzle(dThroat, lastPressure)

            print(f'massOutNozzle: {massOutNozzle}')
            print(f'addedMass: {addedMass}')

            self.massInChamber += addedMass - massOutNozzle

        gWithReg = zip(self.grains, regDepth)
        perGrain = [gr.getCoreVolumeRegression(reg) for gr, reg in gWithReg]
        portVol = sum(perGrain)

        postCombustionVol = self.config.getProperty('postCombustionVol')

        volInChamber = portVol + postCombustionVol

        gasDensity = self.massInChamber / volInChamber

        temp = self.propellant.getCombustionTemp()
        gasConstant = self.propellant.getGasConstant(lastPressure)

        pressure = gasDensity * gasConstant * temp

        return pressure

    def calcIdealThrustCoeff(self, chamberPres, dThroat, exitPres=None):
        """Calculates C_f, the ideal thrust coefficient for the motor's nozzle and propellant, and the given chamber
        pressure. If nozzle exit pressure isn't provided, it will be calculated."""
        if chamberPres == 0:
            return 0

        temp = self.propellant.getCombustionTemp(chamberPres)
        gamma = self.propellant.getCombustionGamma(temp, chamberPres)
        if exitPres is None:
            exitPres = self.nozzle.getExitPressure(gamma, chamberPres)
        ambPres = self.config.getProperty("ambPressure")
        exitArea = self.nozzle.getExitArea()
        throatArea = self.nozzle.getThroatArea(dThroat)

        term1 = (2 * (gamma**2)) / (gamma - 1)
        term2 = (2 / (gamma + 1))**((gamma + 1) / (gamma - 1))
        term3 = 1 - ((exitPres / chamberPres) ** ((gamma - 1) / gamma))

        momentumThrust = (term1 * term2 * term3) ** 0.5
        pressureThrust = ((exitPres - ambPres) * exitArea) / (throatArea * chamberPres)

        return momentumThrust + pressureThrust

    def calcForce(self, chamberPres, dThroat, exitPres=None):
        """Calculates the force of the motor at a given regression depth per grain. Calculates exit pressure by
        default, but can also use a value passed in. This method uses a combination of the techniques described
        in these resources to adjust the thrust coefficient: https://apps.dtic.mil/dtic/tr/fulltext/u2/a099791.pdf
        and http://rasaero.com/dloads/Departures%20from%20Ideal%20Performance.pdf."""
        thrustCoeffIdeal = self.calcIdealThrustCoeff(chamberPres, dThroat, exitPres)
        divLoss = self.nozzle.getDivergenceLosses()
        throatLoss = self.nozzle.getThroatLosses(dThroat)
        skinLoss = self.nozzle.getSkinLosses()
        efficiency = self.nozzle.getProperty('efficiency')
        thrustCoeffAdj = divLoss * throatLoss * efficiency * (skinLoss * thrustCoeffIdeal + (1 - skinLoss))
        thrust = thrustCoeffAdj * self.nozzle.getThroatArea(dThroat) * chamberPres
        return max(thrust, 0)

    def runSimulation(self, callback=None):
        """Runs a simulation of the motor and returns a simRes instance with the results. Constraints are checked,
        including the number of grains, if the motor has a propellant set, and if the grains have geometry errors. If
        all of these tests are passed, the motor's operation is simulated by calculating Kn, using this value to get
        pressure, and using pressure to determine thrust and other statistics. The next timestep is then prepared by
        using the pressure to determine how the motor will regress in the given timestep at the current pressure.
        This process is repeated and regression tracked until all grains have burned out, when the results and any
        warnings are returned."""
        burnoutWebThres = self.config.getProperty('burnoutWebThres')
        burnoutThrustThres = self.config.getProperty('burnoutThrustThres')
        dTime = self.config.getProperty('timestep')

        simRes = SimulationResult(self)

        # Check for geometry errors
        if len(self.grains) == 0:
            aText = 'Motor must have at least one propellant grain'
            simRes.addAlert(SimAlert(SimAlertLevel.ERROR, SimAlertType.CONSTRAINT, aText, 'Motor'))
        for gid, grain in enumerate(self.grains):
            if isinstance(grain, EndBurningGrain) and gid != 0: # Endburners have to be at the foward end
                aText = 'End burning grains must be the forward-most grain in the motor'
                simRes.addAlert(SimAlert(SimAlertLevel.ERROR, SimAlertType.CONSTRAINT, aText, 'Grain ' + str(gid + 1)))
            for alert in grain.getGeometryErrors():
                alert.location = 'Grain ' + str(gid + 1)
                simRes.addAlert(alert)
        for alert in self.nozzle.getGeometryErrors():
            simRes.addAlert(alert)

        # Make sure the motor has a propellant set
        if self.propellant is None:
            alert = SimAlert(SimAlertLevel.ERROR, SimAlertType.CONSTRAINT, 'Motor must have a propellant set', 'Motor')
            simRes.addAlert(alert)
        else:
            for alert in self.propellant.getErrors():
                simRes.addAlert(alert)

        # If any errors occurred, stop simulation and return an empty sim with errors
        if len(simRes.getAlertsByLevel(SimAlertLevel.ERROR)) > 0:
            return simRes

        # Pull the required numbers from the propellant
        density = self.propellant.getProperty('density')

        # Generate coremaps for perforated grains
        for grain in self.grains:
            grain.simulationSetup(self.config)

        # Setup initial values
        perGrainReg = [0 for grain in self.grains]

        # At t = 0, the motor has ignited
        simRes.channels['time'].addData(0)
        igniterPres = self.config.getProperty('igniterPressure')

        simRes.channels['pressure'].addData(igniterPres)
        simRes.channels['force'].addData(0)
        simRes.channels['mass'].addData([grain.getVolumeAtRegression(0) * density for grain in self.grains])
        simRes.channels['massFlow'].addData([0 for grain in self.grains])
        simRes.channels['massFlux'].addData([0 for grain in self.grains])
        simRes.channels['unusedN2O'].addData(0)
        simRes.channels['regression'].addData([0 for grain in self.grains])
        simRes.channels['web'].addData([grain.getWebLeft(0) for grain in self.grains])
        simRes.channels['exitPressure'].addData(0)
        simRes.channels['dThroat'].addData(self.nozzle.getProperty('throat'))

        # Check port/throat ratio and add a warning if it is large enough
        aftPort = self.grains[-1].getPortArea(0)
        if aftPort is not None:
            minAllowed = self.config.getProperty('minPortThroat')
            ratio = aftPort / geometry.circleArea(self.nozzle.props['throat'].getValue())
            if ratio < minAllowed:
                desc = 'Initial port/throat ratio of ' + str(round(ratio, 3)) + ' was less than ' + str(minAllowed)
                simRes.addAlert(SimAlert(SimAlertLevel.WARNING, SimAlertType.CONSTRAINT, desc, 'N/A'))

        # Perform timesteps
        while simRes.shouldContinueSim(burnoutThrustThres):
            n2oMassFlow = self.tank.get_mass_flow_per_iteration(units.convert(simRes.channels['pressure'].getLast(), 'Pa', 'bar'))

            print(f'n2oMassFlow: {n2oMassFlow}')

            unusedN2OMassFlow = n2oMassFlow

            # Calculate regression
            massFlow = 0
            perGrainMass = [0 for grain in self.grains]
            perGrainMassFlow = [0 for grain in self.grains]
            perGrainMassFlux = [0 for grain in self.grains]
            perGrainWeb = [0 for grain in self.grains]
            for gid, grain in enumerate(self.grains):
                if grain.getWebLeft(perGrainReg[gid]) > burnoutWebThres:
                    idealOxiFuelRatio = self.propellant.getProperty('idealOxiFuelRatio')

                    prevAvgTotalMassFlux = simRes.channels['massFlux'].getLast()[gid]

                    # Find the mass flux through the grain based on the mass flow fed into from grains above it
                    perGrainMassFlux[gid] = grain.getAvgTotalMassFlux(massFlow, density, unusedN2OMassFlow, prevAvgTotalMassFlux)

                    # Calculate the combusted N2O for the next grain to have available
                    usedN2OMassFlow = grain.getUsedN2O(unusedN2OMassFlow, density, prevAvgTotalMassFlux, idealOxiFuelRatio)
                    unusedN2OMassFlow -= usedN2OMassFlow

                    # Calculate regression at the current pressure
                    reg = grain.getAvgRegRate(prevAvgTotalMassFlux) * dTime

                    # Find the mass of the grain after regression
                    perGrainMass[gid] = grain.getVolumeAtRegression(perGrainReg[gid]) * density

                    # Add the change in grain mass to the mass flow
                    massFlow += (simRes.channels['mass'].getLast()[gid] - perGrainMass[gid]) / dTime + usedN2OMassFlow

                    # Apply the regression
                    perGrainReg[gid] += reg
                    perGrainWeb[gid] = grain.getWebLeft(perGrainReg[gid])

                perGrainMassFlow[gid] = massFlow

            simRes.channels['unusedN2O'].addData(unusedN2OMassFlow)
            simRes.channels['regression'].addData(perGrainReg[:])
            simRes.channels['web'].addData(perGrainWeb)

            simRes.channels['mass'].addData(perGrainMass)
            simRes.channels['massFlow'].addData(perGrainMassFlow)
            simRes.channels['massFlux'].addData(perGrainMassFlux)

            # Calculate Pressure
            dThroat = simRes.channels['dThroat'].getLast()
            lastPressure = simRes.channels['pressure'].getLast()

            print(f'pressure.data: {simRes.channels["pressure"].data}')
            print(f'lastPressure: {lastPressure}')

            print(f'dThroat: {dThroat}')
            print(f'n2oMassFlow: {n2oMassFlow}')
            print(f'perGrainReg: {perGrainReg}')

            pressure = self.calcIdealPressure(perGrainReg, dThroat, lastPressure, n2oMassFlow)
            simRes.channels['pressure'].addData(pressure)

            print(f'pressure: {pressure}')

            # Calculate Exit Pressure
            temp = self.propellant.getCombustionTemp()
            gamma = self.propellant.getCombustionGamma(pressure)
            exitPressure = self.nozzle.getExitPressure(gamma, pressure)
            simRes.channels['exitPressure'].addData(exitPressure)

            # Calculate force
            force = self.calcForce(simRes.channels['pressure'].getLast(), dThroat, exitPressure)
            simRes.channels['force'].addData(force)

            simRes.channels['time'].addData(simRes.channels['time'].getLast() + dTime)

            # Calculate any slag deposition or erosion of the throat
            if pressure == 0:
                slagRate = 0
            else:
                slagRate = (1 / pressure) * self.nozzle.getProperty('slagCoeff')
            erosionRate = pressure * self.nozzle.getProperty('erosionCoeff')
            change = dTime * ((-2 * slagRate) + (2 * erosionRate))
            simRes.channels['dThroat'].addData(dThroat + change)

            if callback is not None:
                # Uses the grain with the largest percentage of its web left
                progress = max([g.getWebLeft(r) / g.getWebLeft(0) for g, r in zip(self.grains, perGrainReg)])
                if callback(1 - progress): # If the callback returns true, it is time to cancel
                    return simRes

        simRes.success = True

        if simRes.getPeakMassFlux() > self.config.getProperty('maxMassFlux'):
            desc = 'Peak mass flux exceeded configured limit'
            alert = SimAlert(SimAlertLevel.WARNING, SimAlertType.CONSTRAINT, desc, 'Motor')
            simRes.addAlert(alert)

        if simRes.getMaxPressure() > self.config.getProperty('maxPressure'):
            desc = 'Max pressure exceeded configured limit'
            alert = SimAlert(SimAlertLevel.WARNING, SimAlertType.CONSTRAINT, desc, 'Motor')
            simRes.addAlert(alert)

        # Note that this only adds all errors found on the first datapoint where there were errors to avoid repeating
        # errors. It should be revisited if getPressureErrors ever returns multiple types of errors
        for pressure in simRes.channels['pressure'].getData():
            if pressure > 0:
                err = self.propellant.getPressureErrors(pressure)
                if len(err) > 0:
                    simRes.addAlert(err[0])
                    break

        return simRes
