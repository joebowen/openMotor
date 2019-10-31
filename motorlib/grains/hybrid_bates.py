"""Hybrid BATES submodule"""

import numpy as np
import math
import skfmm
from skimage import measure

from ..grain import PerforatedGrain
from .. import geometry
from ..simResult import SimAlert, SimAlertLevel, SimAlertType
from ..properties import FloatProperty

class HybridBatesGrain(PerforatedGrain):
    """The Hybrid BATES grain has a simple cylindrical core. This type is not an FMM grain for performance reasons, as
    the calculations are easy enough to do manually.  This is intended to be used with a hybrid motor"""
    geomName = "Hybrid BATES"
    def __init__(self):
        super().__init__()
        self.props['coreDiameter'] = FloatProperty('Core Diameter', 'm', 0, 1)

    def simulationSetup(self, config):
        self.wallWeb = (self.props['diameter'].getValue() - self.props['coreDiameter'].getValue()) / 2

    def getCorePerimeter(self, regDist):
        return geometry.circlePerimeter(self.props['coreDiameter'].getValue() + (2 * regDist))

    def getCoreVolumeRegression(self, regDist):
        coreLength = self.props['length'].getValue()
        return geometry.circleArea(self.props['coreDiameter'].getValue() + (2 * regDist)) * coreLength

    def getFaceArea(self, regDist):
        outer = geometry.circleArea(self.props['diameter'].getValue())
        inner = geometry.circleArea(self.props['coreDiameter'].getValue() + (2 * regDist))
        return outer - inner

    def getDetailsString(self, preferences):
        lengthUnit = preferences.units.getProperty('m')
        out = 'Length: ' + self.props['length'].dispFormat(lengthUnit)
        out += ', Core: ' + self.props['coreDiameter'].dispFormat(lengthUnit)
        return out

    def getGeometryErrors(self):
        errors = super().getGeometryErrors()
        if self.props['coreDiameter'].getValue() == 0:
            errors.append(SimAlert(SimAlertLevel.ERROR, SimAlertType.GEOMETRY, 'Core diameter must not be 0'))
        if self.props['coreDiameter'].getValue() >= self.props['diameter'].getValue():
            aText = 'Core diameter must be less than grain diameter'
            errors.append(SimAlert(SimAlertLevel.ERROR, SimAlertType.GEOMETRY, aText))
        return errors

    def getAvgTotalMassFlux(self, massIn, density, n2oMassFlow, prevAvgTotalMassFlux):
        coreDiameter = self.props['coreDiameter'].getValue()

        avgPortSurfaceArea = coreDiameter * math.pi * self.props['length'].getValue()

        avgTotalMassFlowRate = self.getAvgTotalMassFlowRate(n2oMassFlow, density, avgPortSurfaceArea, prevAvgTotalMassFlux)

        avgTotalMassFlux = (massIn + avgTotalMassFlowRate) / avgPortSurfaceArea

        return avgTotalMassFlux

    def getAvgTotalMassFlowRate(self, n2oMassFlow, density, avgPortSurfaceArea, prevAvgTotalMassFlux):
        avgTotalMassFlowRate = n2oMassFlow

        tmp = avgTotalMassFlowRate + 10

        iterPrecision = self.props['iterPrecision'].getValue()

        while abs(avgTotalMassFlowRate - tmp) > iterPrecision:
            tmp = avgTotalMassFlowRate
            avgTotalMassFlowRate = self.getAvgFuelMassFlowRate(density, avgPortSurfaceArea, prevAvgTotalMassFlux) + n2oMassFlow

        return avgTotalMassFlowRate

    def getAvgFuelMassFlowRate(self, density, avgPortSurfaceArea, prevAvgTotalMassFlux):
        avgFuelMassFlowRate = density * self.getAvgRegRate(prevAvgTotalMassFlux) * avgPortSurfaceArea

        return avgFuelMassFlowRate

    def getAvgRegRate(self, prevAvgTotalMassFlux):
        portLength = self.props['length'].getValue()
        a = self.props['a'].getValue()
        n = self.props['n'].getValue()
        m = self.props['m'].getValue()

        avgRegRate = a * (prevAvgTotalMassFlux ** n) * (portLength ** m)

        return avgRegRate

    def getUsedN2O(self, n2oMassFlow, density, prevAvgTotalMassFlux, idealOxiFuelRatio):
        # NOTE(joebowen): I'm not sure this is a valid assumption for a hybrid motor given that post-combustion
        # chambers can be used to provide more complete burning.  If not all the combustion happens at the grain,
        # is it possible that the localized grain O/F ratio can't be used to determine 'wasted' N2O?
        coreDiameter = self.props['coreDiameter'].getValue()

        avgPortSurfaceArea = coreDiameter * math.pi * self.props['length'].getValue()

        avgTotalMassFlowRate = self.getAvgFuelMassFlowRate(density, avgPortSurfaceArea, prevAvgTotalMassFlux)

        oxiFuelRatio = n2oMassFlow / avgTotalMassFlowRate

        # TODO: something in here to compare the ideal to the actual ratios and figure out how much oxi is left...

        return n2oMassFlow  # TODO: Just returning n2oMassFlow until I revisit this

    # These two functions have a lot of code reuse, but it is worth it because making BATES an fmmGrain would make it
    # signficantly way slower
    def getFaceImage(self, mapDim):
        mapX, mapY = np.meshgrid(np.linspace(-1, 1, mapDim), np.linspace(-1, 1, mapDim))
        mask = mapX**2 + mapY**2 > 1
        coreMap = np.ones_like(mapX)

        # Normalize core diameter
        coreRadius = (self.props['coreDiameter'].getValue() / (0.5 * self.props['diameter'].getValue())) / 2

        # Open up core
        coreMap[mapX**2 + mapY**2 < coreRadius**2] = 0
        maskedMap = np.ma.MaskedArray(coreMap, mask)

        return maskedMap

    def getRegressionData(self, mapDim, numContours=15, coreBlack=True):
        masked = self.getFaceImage(mapDim)
        regressionMap = None
        contours = []
        contourLengths = {}

        try:
            cellSize = 1 / mapDim
            regressionMap = skfmm.distance(masked, dx=cellSize) * 2
            regmax = np.amax(regressionMap)
            regressionMap = regressionMap[:, :].copy()
            if coreBlack:
                regressionMap[np.where(masked == 0)] = regmax # Make the core black

            for dist in np.linspace(0, regmax, numContours):
                contours.append([])
                contourLengths[dist] = 0
                layerContours = measure.find_contours(regressionMap, dist, fully_connected='high')
                for contour in layerContours:
                    contours[-1].append(contour)
                    contourLengths[dist] += geometry.length(contour, mapDim)

        except ValueError as exc: # If there aren't any contours, do nothing
            print(exc)

        return (masked, regressionMap, contours, contourLengths)
