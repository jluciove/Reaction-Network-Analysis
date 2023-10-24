# ------------------------------------------------
# ------- Reaction Network Model Analysis --------
# by Juan C. Lucio-Vega, 2021

# -------------------------------------------------------------------------------
# Code is used to analyze reaction networks from a legacy kinetic modeling tool
# -------------------------------------------------------------------------------

# Libraries
from models import Models

#  ---- File paths ---- #
propnames_filepath = "csv/propNames.csv"
mixpropdef_filepath = "csv/mixpropdefs.csv"
rxns_filepath = "model/rxn.eqn"
rxnfam_filepath = "model/rxnfamily.dat"
rates_filepath = "model/output/ratefile.dat"
dydt_filepath = "model/output/dydtfile.dat"
properties_filepath = "model/output/proplist.dat"

# ---- Initialize Models ---- #
model = Models()
model.load_rxnfamilies(rxnfam_filepath)
model.load_properties(properties_filepath)
model.load_rxns(rxns_filepath)
model.load_propnames(propnames_filepath)
model.load_mixprop_definitions(mixpropdef_filepath)
model.build_stoich_matrix()
model.mapping_species_to_mixprops()
model.calc_mixprop_rxns()

# ---- Print Information ---- #
# A stoichiometry matrix is an M-by-R matrix, where M equals the total number of species in a model, 
# and R equals the total number of reactions in a model. 
print(model.stoich_matrix)

# print(model.prop_names_dict)
# print(model.mixprop_definitions)
print('Number of species properties: {}'.format(len(model.prop_names_dict)))
print('Number of reactions: {}'.format(model.rxn_no))
print('Number of species: {}'.format(len(model.species)))
print('Number of reaction families: {}'.format(len(model.rxnfamily_names)))
# print('Reaction: {}, Reaction Family: {}'.format(model.rxns[5985].rxn,model.rxns[5985].rxn_family))
# print(model.properties[model.species[0]])

# print(model.map_species_to_mixprops)
print('Print double counted species: {}'.format(model.map_species_to_mixprops['Species000230']))
print('Print lumped rxn network set: \n {}'.format(model.unique_mixprop_rxns))
