from rxns import Rxns
from species import Species
import numpy as np

# Model class
# - Contains species and reaction information
# - Species and reactions are separate classes and 
#   are loaded during run time

class Models:
    species = []
    rxns = []
    stoich_matrix = []
    prop_names_dict = {}
    first_rxn_index = 0
    rxn_no = 0
    rxnfamily_index = []
    rxnfamily_names = []
    mixprop_definitions = {}
    mixprop_values = {}
    properties = {}
    map_species_to_mixprops = {}
    mixprop_rxns = []
    unique_mixprop_rxns = []
    # Class constructor
    def __init__(self):
        self.species = []
        self.rxns = []
        self.prop_names_dict = {}
        self.stoich_matrix = []
        self.first_rxn_index = 0
        self.rxn_no = 0
        self.rxnfamily_index = []
        self.rxnfamily_names = []
        self.mixprop_definitions = {}
        self.properties = {}
        self.mixprop_rxns = []
        self.map_species_to_mixprops = {}
        self.unique_mixprop_rxns = []

    # Load rxns
    def load_rxns(self, rxns_filepath):
        rxn_file = open(rxns_filepath)
        rxns_file = rxn_file.readlines()
        rxn_index = 0
        for rxn in rxns_file[self.first_rxn_index-1:(self.first_rxn_index + self.rxn_no)-1]:
            rxn = rxn.strip('\n')
            rxn = rxn.replace(" ", "")
            # Load into rxn_obj
            rxn_o = Rxns(rxn, self.rxnfamily_index[rxn_index], self.rxnfamily_names[self.rxnfamily_index[rxn_index]-1])
            self.rxns.append(rxn_o)
            self.species = np.concatenate([self.species,rxn_o.species])

            rxn_index = rxn_index + 1

        self.species = list(set(self.species))
        self.species.sort()
        rxn_file.close()

    # Load property names
    def load_propnames(self, propnames_filepath):
        # ---- Read Mix Prop Name File ---- 
        propnames_file = open(propnames_filepath)
        prop_names = propnames_file.readlines()
        lindex = 0
        for name in prop_names:
            name = name.strip('\n')
            self.prop_names_dict[name] = lindex 
            lindex += 1
        propnames_file.close()

    # Load definitions of mixture properties
    def load_mixprop_definitions(self, mixpropdef_filepath):
        mixpropdef_file = open(mixpropdef_filepath)
        mix_prop_def = mixpropdef_file.readlines()
        for definition in mix_prop_def[1:]:
            definition = definition.strip('\n')
            definition = definition.split(',')
            temp_def = []
            logical_def = definition[1].split('&')
            for logical in logical_def:
                logical = logical.replace(" ", "")
                temp_def.append(logical)
            self.mixprop_definitions[definition[0]] = temp_def
        mixpropdef_file.close()

    # Load reaction family information
    def load_rxnfamilies(self, rxnfam_filepath):
        # ---- Read Rxn Fam File ---- 
        rxnfam_file = open(rxnfam_filepath)
        rxnfam = rxnfam_file.readlines()
        locali = 0
        for rf in rxnfam:
            rf = rf.strip('\n')
            if(locali == 0):
                splitrf = rf.split(' ')
                self.first_rxn_index = int(splitrf[0])
                self.rxn_no = int(splitrf[2])
                self.rxnfamily_index.append(int(splitrf[1]))
            else:
                splitrf = rf.split(' ')
                if(locali >= self.rxn_no):
                    self.rxnfamily_names.append(splitrf[0])
                else:
                    self.rxnfamily_index.append(int(splitrf[1]))
            locali = locali + 1
        rxnfam_file.close()
    
    # Build stoichiometry matrix
    def build_stoich_matrix(self):
        for s in self.species:
            s_array = []
            for r in self.rxns:
                s_value = r.stoich_info(s)
                s_array.append(s_value)
            self.stoich_matrix.append(s_array)
        return self.stoich_matrix

    # Load species properties
    def load_properties(self, properties_filepath):
        properties_file = open(properties_filepath)
        properties = properties_file.readlines()
        # Skipping first row and number of reaction families
        for props in properties[1+len(self.rxnfamily_names):]:
            lindex = 0
            lspecies = ''
            lprops = []
            props = props.strip('\n')
            for s_props in props.split(' '):
                if lindex == 0:
                    lspecies = s_props
                else:
                    lprops.append(s_props)
                lindex += 1
            self.properties[lspecies] = lprops
                
    # Calculate mixture properties
    def mapping_species_to_mixprops(self,):
        for s in self.species:
            for key, value in self.mixprop_definitions.items():
                criteria_met = []
                for v in value:
                    foundEQ = v.split('==')
                    foundNEQ = v.split('!=')
                    foundLT = v.split('<')
                    foundGT = v.split('>')
                    foundGTE = v.split('>=')
                    if(len(foundEQ)==2):
                        if s in self.properties:
                            if(float(self.properties[s][self.prop_names_dict[foundEQ[0]]]) == float(foundEQ[1])):
                                criteria_met.append(True)
                            else:
                                criteria_met.append(False)
                    elif (len(foundNEQ)==2):
                        if s in self.properties:
                            if(float(self.properties[s][self.prop_names_dict[foundNEQ[0]]]) != float(foundNEQ[1])):
                                criteria_met.append(True)
                            else:
                                criteria_met.append(False)
                    elif (len(foundLT)==2):
                        if s in self.properties:
                            if(float(self.properties[s][self.prop_names_dict[foundLT[0]]]) < float(foundLT[1])):
                                criteria_met.append(True)
                            else:
                                criteria_met.append(False)
                    elif (len(foundGTE)==2):
                        if s in self.properties:
                            if(float(self.properties[s][self.prop_names_dict[foundGTE[0]]]) >= float(foundGTE[1])):
                                criteria_met.append(True)
                            else:
                                criteria_met.append(False)
                if(all(criteria_met) and len(criteria_met) != 0):
                    if s not in self.map_species_to_mixprops:
                        self.map_species_to_mixprops[s] = [key]
                    else:
                        self.map_species_to_mixprops[s].append(key)

    # Calculate alternative reaction network based on species properties
    def calc_mixprop_rxns(self,):
        for rx in self.rxns:
            mixprop_rxn = ''
            lindex = 0
            findex = len(rx.species_reactants)-1
            for reactant in rx.species_reactants:
                if lindex == findex:
                    if reactant in self.map_species_to_mixprops:
                        mixprop_rxn += self.map_species_to_mixprops[reactant][0]
                    else:
                        mixprop_rxn += reactant
                else:
                    if reactant in self.map_species_to_mixprops:
                        mixprop_rxn += self.map_species_to_mixprops[reactant][0]+'+'
                    else:
                        mixprop_rxn += reactant+'+' 
                lindex += 1
            mixprop_rxn += rx.rxn_direction
            lindex = 0
            findex = len(rx.species_products)-1
            for product in rx.species_products:
                if lindex == findex:
                    if product in self.map_species_to_mixprops:
                        mixprop_rxn += self.map_species_to_mixprops[product][0]
                    else:
                        mixprop_rxn += product
                else:
                    if product in self.map_species_to_mixprops:
                        mixprop_rxn += self.map_species_to_mixprops[product][0]+'+'
                    else:
                        mixprop_rxn += product+'+'
                lindex += 1
            rx.mixprop_rxn = mixprop_rxn
            self.mixprop_rxns.append(mixprop_rxn)
        self.unique_mixprop_rxns = np.unique(self.mixprop_rxns)
