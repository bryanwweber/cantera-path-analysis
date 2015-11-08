from cantera import Solution
from pandas import DataFrame
import numpy as np
import tables


class PathAnalysis(object):
    def __init__(self, chem, save, output, percent, fuel):
        self.chem_file_name = chem
        self.save_file_name = save
        self.output_file_name = output
        self.conversion_percent = percent
        self.fuel = fuel
        self.read_save_file()
        self.generate_table()
        self.write_output()

    def read_save_file(self):
        with tables.open_file(self.save_file_name, mode='r') as save_file:
            table = save_file.root.reactor
            self.time = table.cols.time[:]
            self.mass_fractions = table.cols.mass_fractions[:]
            self.temperature = table.cols.temperature[:]
            self.pressure = table.cols.pressure[:]

    def generate_table(self):
        gas = Solution(self.chem_file_name)
        self.n_reactions = gas.n_reactions
        self.reaction_equations = gas.reaction_equations()
        self.species_names = gas.species_names

        mole_fractions = np.zeros((len(self.time), gas.n_species))
        stoich_diff = gas.product_stoich_coeffs() - gas.reactant_stoich_coeffs()

        for i in range(len(self.time)):
            gas.TPY = self.temperature[i], self.pressure[i], self.mass_fractions[i, :]
            mole_fractions[i, :] = gas.X

        mask = mole_fractions[0, :] > 0
        molar_conversion = (1 - mole_fractions[:, mask]/mole_fractions[0, mask])*100
        initial_species_names = np.array(gas.species_names)[mask]
        fuel_index = np.where(initial_species_names == self.fuel)[0][0]
        conversion_indices = np.where(molar_conversion[:, fuel_index] <= self.conversion_percent)[0]
        rates_of_production = np.zeros((gas.n_species, gas.n_reactions, len(conversion_indices)))

        for i in range(len(self.time[conversion_indices])):
            gas.TPY = self.temperature[i], self.pressure[i], self.mass_fractions[i, :]
            rates_of_production[:, :, i] = gas.net_rates_of_progress*stoich_diff

        integrated_rop = np.trapz(
            rates_of_production, x=self.time[conversion_indices], axis=2
            )

        integrated_prod = np.where(integrated_rop > 0, integrated_rop, np.zeros(integrated_rop.shape))
        integrated_dest = np.where(integrated_rop < 0, integrated_rop, np.zeros(integrated_rop.shape))
        total_prod = np.sum(integrated_prod, axis=1)
        total_dest = -np.sum(integrated_dest, axis=1)

        with np.errstate(divide='ignore', invalid='ignore'):
            self.percent_prod = (integrated_prod.T/total_prod*100).T
            self.percent_prod[self.percent_prod == np.inf] = 0
            self.percent_prod = np.nan_to_num(self.percent_prod)
            self.percent_dest = (integrated_dest.T/total_dest*100).T
            self.percent_dest[self.percent_dest == np.inf] = 0
            self.percent_dest = np.nan_to_num(self.percent_dest)

        self.percent_combined = self.percent_prod + self.percent_dest

    def write_output(self):
        df = DataFrame(
            np.hstack((np.array((range(1, self.n_reactions+1),)).T, self.percent_combined.T)),
            index=self.reaction_equations, columns=(['Reaction No.'] + self.species_names)
            )

        df.to_excel(self.output_file_name)
