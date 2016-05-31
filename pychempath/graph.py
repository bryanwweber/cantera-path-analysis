"""
Functions to produce a graph of the path analysis
"""


def get_species_reactions(species, gas):
    reactions = gas.reactions()
    want = []
    for r in reactions:
        if species in r.products or species in r.reactants:
            want.append(r)

    return want
