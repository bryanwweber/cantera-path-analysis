from cantera import Solution
from pandas import DataFrame
import scipy.sparse as sps
import numpy as np
import tables


class PathAnalysis(object):
    def __init__(self, save, output, percent, fuel, chem=None):
        save_file_name = save
        output_file_name = output
        self.conversion_percent = percent
        self.fuel = fuel

        with tables.open_file(save_file_name, mode='r') as save_file:
            table = save_file.root.reactor
            self.time = table.cols.time[:]
            self.mass_fractions = table.cols.mass_fractions[:]
            self.temperature = table.cols.temperature[:]
            self.pressure = table.cols.pressure[:]
            if chem is None:
                self.cti = save_file.root.cti_file.read().decode('utf-8')

        if chem is None:
            gas = Solution(source=self.cti)
        else:
            gas = Solution(chem)
        self.n_reactions = gas.n_reactions
        self.n_species = gas.n_species
        self.reaction_equations = gas.reaction_equations()
        self.species_names = gas.species_names
        self.stoich_diff = gas.product_stoich_coeffs() - gas.reactant_stoich_coeffs()
        self.molecular_weights = gas.molecular_weights

        self.conversion_indices = self.get_conversion_indices()
        self.rate_of_production = self.generate_rate_of_production(gas)
        self.integ_rop = self.generate_table()
        self.integ_rop.to_excel(output_file_name)

    def get_mole_fractions(self):
        avg_molar_mass = 1/np.sum(self.mass_fractions/self.molecular_weights, axis=1)
        return (self.mass_fractions.T*avg_molar_mass).T/self.molecular_weights

    def get_conversion_indices(self):
        mole_fractions = self.get_mole_fractions()
        mask = mole_fractions[0, :] > 0
        molar_conversion = (1.0 - mole_fractions[:, mask]/mole_fractions[0, mask])*100
        initial_species_names = np.array(self.species_names)[mask]
        fuel_index = np.where(initial_species_names == self.fuel)[0][0]
        return len(np.where(molar_conversion[:, fuel_index] <= self.conversion_percent)[0])

    def generate_rate_of_production(self, gas):
        rate_of_production = sps.lil_matrix(
            (self.conversion_indices*self.n_species, self.n_reactions), dtype=np.float64
        )

        for j in range(self.conversion_indices):
            gas.TPY = self.temperature[j], self.pressure[j], self.mass_fractions[j, :]
            rate_of_production[j*self.n_species:(j+1)*self.n_species, :] = gas.net_rates_of_progress*self.stoich_diff  # noqa

        return rate_of_production

    def generate_table(self):
        integ_rop = DataFrame(index=self.reaction_equations, columns=self.species_names,
                              dtype=np.float64)

        for k, spec in enumerate(self.species_names):
            integ_rop[spec] = np.trapz(y=self.rate_of_production[k::self.n_species].todense(),
                                       x=self.time[:self.conversion_indices], axis=0).T

        total_prod = (integ_rop[integ_rop > 0].fillna(0).values.sum(axis=0))/100
        total_dest = -(integ_rop[integ_rop < 0].fillna(0).values.sum(axis=0))/100

        integ_rop[integ_rop > 0] /= total_prod
        integ_rop[integ_rop < 0] /= total_dest
        return integ_rop
