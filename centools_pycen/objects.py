"""
                                    ---------------------
                                     C E N  -  T O O L S

                                           Objects
                                    ---------------------
"""

########################################################################################################################
########################################################################################################################
# Importing packages and path info

import sys, os, pickle, pandas
import curation

path = os.getcwd() +"/../data/"

########################################################################################################################
########################################################################################################################
# PATHS #

curated_path = path + "curated_data/"
object_path = path + "objects/"

########################################################################################################################
########################################################################################################################
# ESSENTIALITY #

sanger_essentiality, broad_essentiality = curation.essentiality()

# MUTATION INFO #

sanger_mutation, broad_mutation = curation.mutation()

#  EXPRESSION - FPKM INFO #

sanger_expression, broad_expression = curation.expression()

# CNV INFO #

sanger_cnv, broad_cnv = curation.cnv()

########################################################################################################################
########################################################################################################################
# PROJECT CLASSES #

# SANGER #

class Sanger(object):
    def __init__(self, model_name):
        self.model_name = model_name
        self.sanger = None
        self.EXPRESSION, self.MUTATION, self.CNV, self.FUSION, self.DRUG, self.ESSENTIALITY =\
            False,False,False,False,False,False
        self.mutations, self.oncogenic_mutations, self.hotspot_mutations = None, None, None
        self.mutated_genes, self.oncogenic_mutated_genes = None, None
        self.expression = None
        self.cnv = None
        self.fusion = None
        self.growth_property = None
        self.msi = None
        self.ploidy = None
        self.drug = None
        self.tissue, self.site = None, None
        self.cancer, self.cancer_subtype = None, None
        self.broad = None
        self.essentiality = None


# BROAD #

class Broad(object):
    def __init__(self, broad):
        self.broad = broad
        self.model_name = None
        self.EXPRESSION, self.MUTATION, self.CNV, self.FUSION, self.DRUG, self.ESSENTIALITY =\
            False,False,False,False,False,False
        self.mutations, self.oncogenic_mutations, self.hotspot_mutations = None, None, None
        self.mutated_genes, self.oncogenic_mutated_genes = None, None
        self.expression = None
        self.cnv = None
        self.drug = None
        self.fusion = None
        self.site, self.tissue, self.growth_property = None, None, None
        self.cancer, self.cancer_subtype = None, None
        self.sanger = None
        self.essentiality = None


########################################################################################################################
########################################################################################################################
# DE-SERIALISATION OF THE OBJECTS #

def deserialisation_project():
    sanger_obj = pickle.load(open(object_path + "SANGER_OBJECT.pkl", "rb"))
    broad_obj = pickle.load(open(object_path + "BROAD_OBJECT.pkl", "rb"))
    return sanger_obj, broad_obj


########################################################################################################################
########################################################################################################################



