import pytest
import cantera as ct
import pychempath.graph as pcp


@pytest.fixture
def setup_cantera():
    return ct.Solution('gri30.xml')


def test_get_species_reactions(setup_cantera):
    species = 'CH4'
    want = pcp.get_species_reactions(species, setup_cantera)
    for r in want:
        assert (species in r.products or species in r.reactants)
