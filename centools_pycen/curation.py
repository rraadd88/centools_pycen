"""
                                    ---------------------
                                     C E N  -  T O O L S

                                      Data  Preparation
                                    ---------------------
"""
########################################################################################################################
########################################################################################################################
# Importing packages and path info

import sys, os, pickle, pandas

path = os.getcwd() +"/../data/"

########################################################################################################################
########################################################################################################################
# PATHS #

raw_path = path + "raw_data/"
curated_path = path + "curated_data/"

########################################################################################################################
########################################################################################################################
# ESSENTIALITY DATA FRAMEs #
# Read curated essentiality tables

def essentiality():
    sanger_essentiality = pandas.read_csv(curated_path + "Sanger_Essentiality.csv", index_col=0)
    broad_essentiality = pandas.read_csv(curated_path + "Broad_Essentiality.csv", index_col=0)
    return sanger_essentiality, broad_essentiality

########################################################################################################################
########################################################################################################################
# CELL LINE MAPPING AND DISEASE INFO #

# Read curated conversion tables
def conversion():
    ess_common_df = pandas.read_csv(curated_path + "Model_ID_Conversion.csv", index_col = 0)
    return ess_common_df

# Read manually curated file

def manual_disease():
    disease_map = pandas.read_csv(curated_path + "Manually_curated_tissue_cancer_map.csv", index_col = 0)
    return disease_map

########################################################################################################################
########################################################################################################################
# MUTATION INFO #
# Read curated mutation tables

def mutation():
    sanger_mutation = pandas.read_csv(curated_path + "Sanger_Mutation.csv", index_col=0)
    broad_mutation = pandas.read_csv(curated_path + "Broad_Mutation.csv", index_col=0)
    return sanger_mutation, broad_mutation

########################################################################################################################
########################################################################################################################
#  EXPRESSION - FPKM/TPM INFO #
# Read the curated expression files

def expression():
    sanger_expression = pandas.read_csv(curated_path + "Sanger_Expression.csv", index_col=0)
    broad_expression = pandas.read_csv(curated_path + "Broad_Expression.csv", index_col=0)
    return sanger_expression, broad_expression

########################################################################################################################
########################################################################################################################
# FUSION INFO #
# Read the curated data

def fusion():
    sanger_fusion = pandas.read_csv(curated_path + "Sanger_Fusion.csv", index_col=0)
    broad_fusion = pandas.read_csv(curated_path + "Broad_Fusion.csv", index_col=0)
    return sanger_fusion, broad_fusion


########################################################################################################################
########################################################################################################################
# CNV INFO #
# Read the curated data

def cnv():
    sanger_cnv = pandas.read_csv(curated_path + "Sanger_CNV.csv", index_col=0)
    broad_cnv = pandas.read_csv(curated_path + "Broad_CNV.csv", index_col=0)
    return sanger_cnv, broad_cnv

########################################################################################################################
########################################################################################################################
# DRUG INFO #
# Read the curated data

def drug():
    sanger_drug = pandas.read_csv(curated_path + "Sanger_Drug.csv", index_col=0)
    broad_drug = pandas.read_csv(curated_path + "Broad_Drug.csv", index_col=0)
    return sanger_drug, broad_drug


########################################################################################################################
########################################################################################################################

