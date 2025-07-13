# Reaction_Modeler
Python code to simulate and visualize reaction progression over time. 

Note: Reaction inputs are lists of single steps or equilibria with each species designated as a single letter code. A colon separates out the reaction or species from its rate or concentration. 
Coding examples displayed below. Great for introductory demonstrations of reaction kinetics and equilibria. Likely insufficient integration accuracy to demonstrate known complex oscillating reactions like Iodine Clock.

**Example 1: A <-> B**
___________________________________________________________________________
Code Snippet for Programming
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
```
___________________________________________________________________________
Printout for Simulation
```r
=== Initial Concentrations ===
A: 100.0 mM
B: 20.0 mM

=== List of All Reactions ===
R1: A -> B    (rate = k1 * [A], k1 = 0.1)
R2: B -> A    (rate = k2 * [B], k2 = 0.1)
```
___________________________________________________________________________
![Description](https://github.com/jawolfe97/Reaction_Modeler/blob/main/AtoB.svg)


**Example 2: Complex Rxn**
In this example, it's not obvious at first, but C is the limiting reagent. 
The initial spike in B can't be sustained and since nothing interconverts with C, the first reaction slowly dies out and A becomes the only remaining species.
___________________________________________________________________________
```r
# === USER INPUT ===
reaction_input = "2A + C -> B : 0.01, 2D -> A : 0.01, B <-> D: 0.01//0.01, E -> 2D: 0.1"
concentration_input = "A: 50, B: 50, C: 50, D: 50, E: 50"

# === UNITS ===
time_unit = "s"
conc_unit = "mM"

# === Simuplation Parameters ===
simulation_length = 100 #time_unit
simulation_timestep = 0.01 #time_unit
```
___________________________________________________________________________
![Description](https://github.com/jawolfe97/Reaction_Modeler/blob/main/ComplexRxn_100s.svg)
___________________________________________________________________________
```r
# === USER INPUT ===
reaction_input = "2A + C -> B : 0.01, 2D -> A : 0.01, B <-> D: 0.01//0.01, E -> 2D: 0.1"
concentration_input = "A: 50, B: 50, C: 50, D: 50, E: 50"

# === UNITS ===
time_unit = "s"
conc_unit = "mM"

# === Simuplation Parameters ===
simulation_length = 1000 #time_unit
simulation_timestep = 0.01 #time_unit
```
___________________________________________________________________________
![Description](https://github.com/jawolfe97/Reaction_Modeler/blob/main/ComplexRxn_1000s.svg)

## Acknowledgments

Parts of this code were generated or assisted by OpenAI's ChatGPT (GPT-4), based on user prompts. All generated code was reviewed and modified for correctness and performance.

