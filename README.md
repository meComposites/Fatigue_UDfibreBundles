# Fatigue_UDfibreBundles

Fatigue code for unidirectional composites under tension-tension loading

This code is related to the following publication:
"A computationally-efficient micromechanical model for the fatigue life of unidirectional composites under tension-tension loading"

DOI: https://doi.org/10.1016/j.ijfatigue.2018.05.017


Main.m is the main code of the fatigue model. This is where all the material properties and fatigue inputs are defined. 
The definition of the number of fatigue cycles to be modeled is implemented in a jump cycle strategy. For this the following is required:

The size of  cycle intervals  is defined in  DeltaCycles.
Ex: 3 intervals with a cycle interval size of 1 10 and 100
DeltaCycles = [1 10 1000 ]. 

The number of icrements done for each cycle interval is defined in CycleIntervals. 
Ex: 100 increments with a cycle interval of 1, and 3000 increments with a cycle interval of 10 and 1000 increments with a cycle interval of 1000.
CycleIntervals = [100 3000 1000]

The code automatically plots:

1 - The static strength distribution for each bundle size i (with 2^i fibres)
2 - a deterministic SN curve for a pre-defined bundle size.
2 - a stochastic SN curve for a pre-defined bundle size.


