# Reaction Network Analysis

Python code developed to analyze results from a legacy kinetic modeling tool. The code parses information from a reaction network and property files. The goal of the analysis is to understand where molecular flow is moving during kinetic parameter estimation.

Analysis:
 - Stoichmetric Matrix
 - Species
 - Bulk Property Test - when converting molecular information to bulk properties, are we double counting molecules? 
 - Convert molecular-level model to lumped model based on bulk property information (The molecular model should match the lumped model)

Important Files:
  * reaction_network_analysis.py - Main file
  * models.py - Parse data from files and transfer to class data structure
  * rxns.py - Parse reaction info and store in class
