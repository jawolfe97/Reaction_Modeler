# === IMPORTS ===
import json
import re
from collections import defaultdict
import numpy as np
import matplotlib.pyplot as plt

# === FUNCTION DEFINITIONS ===

def parse_species_with_stoich(species_list):
    """Parse stoichiometry from expressions like '2A' or 'A'."""
    parsed = []
    for item in species_list:
        item = item.strip()
        match = re.match(r'^(\d+)?\s*([A-Za-z_][A-Za-z0-9_]*)$', item)
        if not match:
            raise ValueError(f"Invalid species format: '{item}'")
        stoich = int(match.group(1)) if match.group(1) else 1
        species = match.group(2)
        parsed.append((species, stoich))
    return parsed

def parse_initial_concentrations(conc_string):
    """Parse initial concentrations like 'A: 10, B: 12' into a dictionary."""
    conc_dict = {}
    for item in conc_string.split(','):
        item = item.strip()
        if not item:
            continue
        match = re.match(r'^([A-Za-z_][A-Za-z0-9_]*):\s*([\d.]+)$', item)
        if not match:
            raise ValueError(f"Invalid concentration format: '{item}'")
        species, value = match.groups()
        conc_dict[species] = float(value)
    return conc_dict

def parse_reactions(reaction_string, conc_string=""):
    reactions = [r.strip() for r in reaction_string.split(',')]
    initial_concs = parse_initial_concentrations(conc_string)
    species_set = set(initial_concs.keys())
    species_conc = defaultdict(lambda: 0.0, initial_concs)

    reaction_list = []
    readable_reactions = []
    reaction_id = 1

    for r in reactions:
        match = re.match(r'(.+?)(<->|->|<-)(.+?)(?::\s*([\d./]+))?$', r.strip())
        if not match:
            raise ValueError(f"Invalid reaction format: {r}")

        lhs, arrow, rhs, rate_info = match.groups()
        lhs_parsed = parse_species_with_stoich(lhs.split('+'))
        rhs_parsed = parse_species_with_stoich(rhs.split('+'))
        species_set.update(s for s, _ in lhs_parsed + rhs_parsed)

        def add_reaction(rid, reactants, products, rate):
            rate_law_terms = []
            for s, n in reactants:
                term = f"[{s}]" if n == 1 else f"[{s}]**{n}"
                rate_law_terms.append(term)
            rate_law = f"k{rid} * " + ' * '.join(rate_law_terms) if rate_law_terms else f"k{rid}"

            reaction_entry = {
                    "id": f"R{rid}",
                    "reactants": [{"species": s, "stoichiometry": n} for s, n in reactants],
                    "products": [{"species": s, "stoichiometry": n} for s, n in products],
                    "rate_law": rate_law,
                    "parameters": {f"k{rid}": float(rate)}
            }
            reaction_list.append(reaction_entry)

            lhs_str = ' + '.join(f"{n if n > 1 else ''}{s}" for s, n in reactants)
            rhs_str = ' + '.join(f"{n if n > 1 else ''}{s}" for s, n in products)
            readable_reactions.append(
                    f"R{rid}: {lhs_str} -> {rhs_str}    (rate = {rate_law}, k{rid} = {float(rate)})"
            )

        if arrow == '->':
            if rate_info is None:
                raise ValueError(f"Missing rate for irreversible reaction: {r}")
            add_reaction(reaction_id, lhs_parsed, rhs_parsed, rate_info)
            reaction_id += 1

        elif arrow == '<-':
            if rate_info is None:
                raise ValueError(f"Missing rate for irreversible reverse reaction: {r}")
            add_reaction(reaction_id, rhs_parsed, lhs_parsed, rate_info)
            reaction_id += 1

        elif arrow == '<->':
            if rate_info is None or '//' not in rate_info:
                raise ValueError(f"Missing or malformed rate for reversible reaction: {r}")
            kf, kr = rate_info.split('//')
            add_reaction(reaction_id, lhs_parsed, rhs_parsed, kf)
            reaction_id += 1
            add_reaction(reaction_id, rhs_parsed, lhs_parsed, kr)
            reaction_id += 1

    species_list = [{"id": s, "initial_concentration": species_conc[s]} for s in sorted(species_set)]

    graph = {
        "species": species_list,
        "reactions": reaction_list
    }

    return graph, readable_reactions

# === USER INPUT ===
reaction_input = "A <-> B : 0.1//0.1"
concentration_input = "A: 100, B: 20"
 
# === UNITS ===
time_unit = "s"
conc_unit = "mM"

# === Simuplation Parameters ===
simulation_length = 100 #time_unit
simulation_timestep = 0.01 #time_unit

# === PARSE INPUT ===
parsed_graph, readable_rxns = parse_reactions(reaction_input, concentration_input)

# === DISPLAY ===
print("=== JSON Representation ===")
print(json.dumps(parsed_graph, indent=2))

print("\n=== Initial Concentrations ===")
for s in parsed_graph["species"]:
    print(f"{s['id']}: {s['initial_concentration']} {conc_unit}")

print("\n=== List of All Reactions ===")
for rxn in readable_rxns:
    print(rxn)

# === SETUP ===
graph = parsed_graph
species_ids = [s['id'] for s in graph['species']]
initial_concs = np.array([s['initial_concentration'] for s in graph['species']], dtype=float)
species_index = {s: i for i, s in enumerate(species_ids)}

dt = simulation_timestep        # timestep (in time_unit)
t_max = simulation_length       # max time (in time_unit)
n_steps = int(t_max / dt)
t_values = np.linspace(0, t_max, n_steps)

y_values = np.zeros((len(species_ids), n_steps))
y_values[:, 0] = initial_concs

# === EULER INTEGRATION ===
for step in range(1, n_steps):
    y_prev = y_values[:, step - 1]
    dydt = np.zeros_like(y_prev)

    for reaction in graph['reactions']:
        k_name, k_val = next(iter(reaction['parameters'].items()))
        rate = k_val
        for reactant in reaction['reactants']:
            idx = species_index[reactant['species']]
            rate *= y_prev[idx] ** reactant['stoichiometry']
        for reactant in reaction['reactants']:
            idx = species_index[reactant['species']]
            dydt[idx] -= rate * reactant['stoichiometry']
        for product in reaction['products']:
            idx = species_index[product['species']]
            dydt[idx] += rate * product['stoichiometry']

    y_new = y_prev + dt * dydt
    y_values[:, step] = np.clip(y_new, 0, None)

# === PLOTTING ===
plt.figure(figsize=(10, 6))
for i, species in enumerate(species_ids):
    plt.plot(t_values, y_values[i], label=f"[{species}]")

plt.xlabel(f"Time ({time_unit})")
plt.ylabel(f"Concentration ({conc_unit})")
plt.title("Species Concentrations Over Time")
plt.legend()
plt.grid(True)
plt.tight_layout()
plt.show()
