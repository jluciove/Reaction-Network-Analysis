import re

class Rxns:
    species = []
    reactants = []
    products = []
    stoichiometry = []
    stoich_dict = {}
    rxn_family = ''
    rxn_direction = ''
    rxn_family_index = 0
    mixprop_rxn = ''
    # Class constructor
    def __init__(self, rxn, index, rxn_fam):
        self.rxn = rxn
        self.species = []
        self.stoichiometry = []
        self.stoich_dict = {}
        self.rxn_family_index = index
        self.rxn_family = rxn_fam
        self.rxn_direction = ''
        self.species_reactants = []
        self.species_products = []
        self.mixprop_rxn = ''
        self.parse_rxn()

    # Parse reactions method  
    def parse_rxn(self):
        species = []
        stoichiometry = []
        # First check if reversible
        split_rxn = self.rxn.split('<->')
        if len(split_rxn) == 2:
            self.rxn_direction = '<->'
            # Split reactants
            reactants = split_rxn[0]
            reactants = reactants.replace(" ", "")
            reactants = reactants.split('+')
            for reactant in reactants:
                result = re.match(r"(?P<stoich>\d+)(?P<species>\w+)", reactant)
                if result:
                    self.stoichiometry.append(-1*int(result.group('stoich')))
                    self.species.append(result.group('species'))
                    self.species_reactants.append(result.group('species'))
                else:
                    self.stoichiometry.append(-1)
                    self.species.append(reactant)
                    self.species_reactants.append(reactant)
            # Split products
            products = split_rxn[1]
            products = products.replace(" ", "")
            products = products.split('+')
            for product in products:
                result = re.match(r"(?P<stoich>\d+)(?P<species>\w+)", product)
                if result:
                    self.stoichiometry.append(int(result.group('stoich')))
                    self.species.append(result.group('species'))
                    self.species_products.append(result.group('species'))
                else:
                    self.stoichiometry.append(1)
                    self.species.append(product)
                    self.species_products.append(product)
        else:
            # Resction is irreversible
            self.rxn_direction = '->'
            split_rxn = self.rxn.split('->')
            # Split reactants
            reactants = split_rxn[0]
            reactants = reactants.replace(" ", "")
            reactants = reactants.split('+')
            for reactant in reactants:
                result = re.match(r"(?P<stoich>\d+)(?P<species>\w+)", reactant)
                if result:
                    self.stoichiometry.append(-1*int(result.group('stoich')))
                    self.species.append(result.group('species'))
                    self.species_reactants.append(result.group('species'))
                else:
                    self.stoichiometry.append(-1)
                    self.species.append(reactant)
                    self.species_reactants.append(reactant)
            # Split products
            products = split_rxn[1]
            products = products.replace(" ", "")
            products = products.split('+')
            for product in products:
                result = re.match(r"(?P<stoich>\d+)(?P<species>\w+)", product)
                if result:
                    self.stoichiometry.append(int(result.group('stoich')))
                    self.species.append(result.group('species'))
                    self.species_products.append(result.group('species'))
                else:
                    self.stoichiometry.append(1)
                    self.species.append(product)
                    self.species_products.append(product)

        self.stoich_dict = { sp : st for sp, st in zip(self.species, self.stoichiometry) }

    # Stoich info return method
    def stoich_info(self, species_name):
        if species_name in self.stoich_dict.keys():
            return self.stoich_dict[species_name]
        else:
            return 0