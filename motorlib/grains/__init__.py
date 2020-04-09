from .endBurner import *
from .bates import *
from .finocyl import *
from .moonBurner import *
from .star import *
from .xCore import *
from .cGrain import *
from .dGrain import *
from .rodTube import *
from .custom import *
from .hybrid_bates import *

# Generate grain geometry name -> constructor lookup table
grainTypes = {}

apcpGrainClasses = [BatesGrain, EndBurningGrain, Finocyl, MoonBurner, StarGrain, XCore, CGrain, DGrain, RodTubeGrain, CustomGrain]
hybridGrainClasses = [HybridBatesGrain]

for grainType in hybridGrainClasses:
    grainTypes[grainType.geomName] = grainType
