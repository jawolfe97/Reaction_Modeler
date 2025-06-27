# Reaction_Modeler
Python code to simulate and visualize reaction progression over time

**Example 1: A <-> B**
___________________________________________________________________________
```r
# === USER INPUT ===
reaction_input = "A <-> B : 0.1//0.1"
concentration_input = "A: 100, B: 20"
 
# === UNITS ===
time_unit = "s"
conc_unit = "mM"

# === Simuplation Parameters ===
simulation_length = 100 #time_unit
simulation_timestep = 0.01 #time_unit


