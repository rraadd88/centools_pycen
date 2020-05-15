########################################################################################################################
########################################################################################################################
"""
                                       -----------------
                                       C E N - T O O L S
                                       -----------------
"""

########################################################################################################################
########################################################################################################################
# PACKAGES #

import sys, os, scipy, pandas, numpy, pickle, networkx, argparse, warnings
from scipy import stats
from scipy.stats import mannwhitneyu
warnings.simplefilter(action='ignore', category=FutureWarning)


print("Packages are ready to use!")
print("\n****************************************************************************\n")

########################################################################################################################
########################################################################################################################
# INPUT #

print("Taking inputs...\n")

def check_file_exists(file):
    if not os.path.exists(file): parser.error("The file %s does not exist in %s!" %(file,os.getcwd()+"/"))
    else: return pandas.read_csv(os.getcwd()+ "/" + file, index_col = 0)


def take_input():
    parser = argparse.ArgumentParser(prog="PyCEN",
                                     usage="%(prog)s [inputs]",
                                     formatter_class=argparse.RawDescriptionHelpFormatter,
                                     description=("""
                                     -----------------
                                     C E N - T O O L S
                                         
                                         P y C E N
                                     -----------------

Terms and conditions CEN-tools is an integrative tool that identifies the underlying context that is associated with the 
essentiality of a gene from large-scale genome-scale CRISPR screens.

Copyright © 2020 EMBL- European Bioinformatics Institute

This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public 
License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later 
version. This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the 
implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more 
details. Neither the institution name nor the name CEN-tools can be used to endorse or promote products derived from 
this software without prior written permission. For written permission, please contact petsalaki@ebi.ac.uk.

Products derived from this software may not be called CEN-tools nor may CEN-tools appear in their names without prior 
written permission of the developers. You should have received a copy of the GNU General Public License along with 
this program. If not, see gnu.org/licenses/.

                                     """))

    for grp in parser._action_groups:
        if grp.title == "optional arguments": grp.title = "Inputs"
        elif "positional arguments": grp.title = "Mandatory Inputs"

    # QUERY GENE(S)

    # The query can be done with one gene or a list of genes in a csv file!

    gene_info = parser.add_mutually_exclusive_group(required = True)
    gene_info.add_argument("-gene", dest = "GENES", type = lambda x: [x],
                           help = "The interested gene whose CEN will be created!")
    gene_info.add_argument("-genefile", dest = "GENES",type = lambda x: list(check_file_exists(x).index),
                           help="The interested gene file whose CEN will be created!")

    # PROJECT

    # The project to used for essentiality, cell line separation, context unless otherwise specified!
    parser.add_argument("-project", dest = "PROJECT", choices = ["SANGER", "BROAD"], default = "BROAD",
                        help = "The interested project whose data will be used! (BROAD or SANGER)")

    # ESSENTIALITY

    # The user can use the pre-defined screens (SANGER or BROAD),
    # or can load a different essentiality data!
    # If the user will add essentiality information, the file should be provided!

    parser.add_argument("-add_essentiality", action = "store_true", dest = "ADD_ESSENTIALITY",
                        help = "Boolean value to represent if the user will provide essentiality data.")

    parser.add_argument("-essfile", dest = "ESSENTIALITY", type = lambda x: check_file_exists(x),
                        required = "-add_essentiality" in sys.argv,
                        help = "The user-defined essentiality data (rows as gene names and columns as cell lines)")

    # CONTEXT

    # The user can use the pre-defined context files or can load different files!
    parser.add_argument("-add_context", action = "store_true", dest = "ADD_CONTEXT",
                        help = "Boolean value to represent if the user will provide a new context data"
                               "rather than using expression, amplification, depletion, or mutation.")

    # If the ADD_CONTEXT is provided as TRUE then the file should be given with its data structure!

    parser.add_argument("-contextfile", dest = "CONTEXT", type = lambda x: check_file_exists(x),
                        required="-add_context" in sys.argv,
                        help="The context data table (rows as gene names and columns as cell lines)")

    parser.add_argument("-contextstr", dest = "CONTEXT_STR", choices = ["continuous", "discrete"],
                        required="-add_context" in sys.argv,
                        help= "The structure of the data provided as a new context.")

    # The user should select a context!
    # If the add_context is TRUE than it should be written --data or
    # otherwise any other pre-defined context should be given!

    context_data_info = parser.add_mutually_exclusive_group()

    context_data_info.add_argument("-exp", action='store_true', dest= "EXPRESSION",
                                   help = "Use if the context is pre-defined expression.")
    context_data_info.add_argument("-compess", action='store_true', dest= "COMPARE_ESSENTIALITY",
                                   help = "Use if the context is to compare the"
                                          "essentialities of user define cell line lists.")
    context_data_info.add_argument("-mut", action='store_true', dest= "MUTATION",
                                   help = "Use if the context is pre-defined mutation.")
    context_data_info.add_argument("-amp", action='store_true', dest= "AMPLIFICATION",
                                   help = "Use if the context is pre-defined amplification (CNV).")
    context_data_info.add_argument("-dep", action='store_true', dest= "DEPLETION",
                                   help="Use if the context is pre-defined depletion (CNV).")

    # DEPENDANT GENE(S)

    # If the user want to search for the dependence of a specific gene/gene list to query gene/gene_list
    # They can provide it.
    # Otherwise, all the genes in the essentiality data will be searched!

    parser.add_argument("-add_dependant", action="store_true", dest = "ADD_DEPENDANT",
                        help = "Boolean value. If the user wants to provide a dependant gene/gene list")

    dependant_gene_info= parser.add_mutually_exclusive_group()
    dependant_gene_info.add_argument("-cogene", dest="CoGENES",
                                     help=" Specified dependant gene!", type = lambda x: [x])
    dependant_gene_info.add_argument("-cogenefile", dest="CoGENES",
                                     help=" Specified dependant gene!",
                                     type = lambda x: list(check_file_exists(x).index))

    # GROUPING

    # CEN-tools is offering two types of grouping for the cancer cell lines.
    # The user can select "tissue" or "cancer".

    parser.add_argument("-group", dest = "SITE", choices=["tissue", "cancer"], default="tissue",
                        help= "According to what cell lines will be grouped. (default: %(default)s)")

    # If the user wants to search for a specific tissue/cancer type, they can provide it by name of file.
    # Otherwise, all the tissues/cancer types will be searched!

    group_info= parser.add_mutually_exclusive_group()
    group_info.add_argument("-groupname", dest="GROUPS",
                            help="Specified tissue/cancer type!", type = lambda x: [x])
    group_info.add_argument("-groupfile", dest="GROUPS", type = lambda x: list(check_file_exists(x).index),
                            help="Specified tissues/cancer types!")

    # If the user wants to separate cell lines manually, they need to specify it!
    parser.add_argument("-separate", action="store_true", dest="SEPARATION",
                        help="Boolean value. If the user wants to provide separated cell lines lists!")

    parser.add_argument("-interested_group",dest="INTERESTED",help="Specified interested group of cell lines!",
                        type = lambda x: list(check_file_exists(x).index))
    parser.add_argument("-other_group", dest="OTHER", help="Specified other group of cell lines!",
                        type = lambda x: list(check_file_exists(x).index))

    # THRESHOLDS

    # The user can provide the limitation of the statistical tests!

    parser.add_argument("-limcor", default=0.5, type=float, dest= "CORLIMIT",
                        help= "The minimum threshold for absolute correlation!")
    parser.add_argument("-limp", default = 0.05, type = float, dest = "PLIMIT",
                        help = "The maximum threshold for p-value!")

    # OUTPUT

    # Output folder can be specified! Otherwise the script write into current directory!

    parser.add_argument("-o", dest = "OUTPUT_PATH", default = os.getcwd() + "/",
                        help="The path for output. If not specified the current directory will be used!")

    parser.add_argument("-ofile", dest = "OUTPUT_FILE", default = "CEN",
                        help="The output file name, if not specified CEN will be used!")

    # MUTATION
    mut_info = parser.add_mutually_exclusive_group()
    mut_info.add_argument("-hotspot", action="store_true", dest= "HOTSPOT",
                          help = "Cell lines will be separated according to presence "
                                 "of mutations on genes having hotspot mutations!")
    mut_info.add_argument("-oncogenic", action="store_true", dest= "ONCOGENIC",
                          help = "Cell lines will be separated according to presence "
                                 "of mutations on oncogenic genes!")

    parser.print_help()

    parsed_input = parser.parse_args()
    input_dict = vars(parsed_input)

    return input_dict

args = take_input()

########################################################################################################################
########################################################################################################################
# PATHS #

# The user should be on pyCEN directory:

if os.getcwd().split("/")[-1] == "centools_pycen": path = os.getcwd()
else: sys.exit("The script should be run inside centools_pycen directory!")

print("The path is %s\n" %path)

# If RESULT_PATH was provided , it will be used!
# If not current folder will be tried!

if "RESULT_PATH" in args.keys():
    # Create the output folder
    os.system("mkdir %s/%s" %(args["RESULT_PATH"], "OUTPUT"))
    result_path = args["RESULT_PATH"] + "/OUTPUT/"

else:
    os.system("mkdir %s/%s" % (path, "OUTPUT/"))
    result_path = path + "/OUTPUT/"

print("The result will be written in %s\n" %result_path)

# Add path to system for modules!
if path + "centools_pycen" not in sys.path: sys.path.append(path + "/centools_pycen")

print("The path of centools_pycen modules are added to the python path!")

########################################################################################################################
########################################################################################################################
# CLASSES #

print("Importing the objects...\n")

if args["PROJECT"]:
    import objects
    from objects import Sanger, Broad
    sanger_obj, broad_obj = objects.deserialisation_project()

else: sanger_obj, broad_obj = None, None

if sanger_obj is None and broad_obj is None: sys.exit("Cannot import objects!")

########################################################################################################################
########################################################################################################################
# ESSENTIALITY #

essentiality_data, obj, project = False, None, None

if args["PROJECT"]:
    try:
        from curation import essentiality
        project = args["PROJECT"]
        obj = sanger_obj if project == "SANGER" else broad_obj
        essentiality_df = essentiality()[0] if project == "SANGER" else essentiality()[1]
        essentiality_data = True

        print("The essentiality data from %s project is loaded!\n" % project)

    except ImportError: sys.exit("The error occurred when importing data!")

else: obj = None

if args["ESSENTIALITY"] is not None:
    essentiality_df = args["ESSENTIALITY"]
    essentiality_data = True
    print("The essentiality data from the user is loaded!\n")

if project is None and args["ESSENTIALITY"] is None:
    sys.exit("No essentiality information was found!!")

if not essentiality_data: sys.exit("No essentiality information was found!!")

########################################################################################################################
########################################################################################################################
# CONTEXT #

context, context_df, context_str = None, None, None
EXP, MUT, AMP, DEP, CNT, COMPESS = False, False, False, False, False, False
HOTSPOT, ONCOGENE = False, False

if not args["ADD_CONTEXT"]:
    if project:
        arg_context = [c for c in ["EXPRESSION", "MUTATION", "AMPLIFICATION", "DEPLETION",
                                   "COMPARE_ESSENTIALITY"] if args[c]]
        if len(arg_context) == 1:
            context = arg_context[0]
            if context == "EXPRESSION":
                from curation import expression
                context_df = expression()[0] if project == "SANGER" else expression()[1]
                context_str = "continuous"
                EXP = True
            elif context == "MUTATION":
                from curation import mutation
                context_df = mutation()[0] if project == "SANGER" else mutation()[1]
                context_str = "discrete"
                MUT = True
                if "HOTSPOT" in args.keys(): HOTSPOT = True
                elif "ONCOGENIC" in args.keys(): ONCOGENE = True
                else: sys.exit("No Hotspot or oncogenic information was given!")
            elif context in ["AMPLIFICATION", "DEPLETION"]:
                from curation import cnv
                context_df = cnv()[0] if project == "SANGER" else cnv()[1]
                context_str = "discrete"
                if context == "AMPLIFICATION": AMP = True
                elif context == "DEPLETION": DEP = True
            elif context =="COMPARE_ESSENTIALITY":
                from curation import essentiality
                context_df = essentiality()[0] if project == "SANGER" else essentiality()[1]
                context_str = "discrete"
                COMPESS = True
        elif len(arg_context) == 0: sys.exit("One of the pre-defined context should be given!")
        else: sys.exit("Only one context should be given!")

    else: sys.exit("No context information was found!!")

else:
    if args["ADD_CONTEXT"]:
        if args["CONTEXT"] is not None and args["CONTEXT_STR"] is not None:
            context = "user_defined"
            context_df = args["CONTEXT"]
            context_str = args["CONTEXT_STR"]
            CNT = True
        else: sys.exit("There should be new context file with the data structure information!")
    else: sys.exit("Context information cannot be found!")

########################################################################################################################
########################################################################################################################
# GENE(S) #

# First check if the query gene(s) is in the essentiality data!
# Second check if the query gene(s) is in the context data!

query_genes = []
if args["GENES"]:
    ess_data_available, context_data_available = {}, {}

    for gene in args["GENES"]:
        if gene in list(essentiality_df.index): ess_data_available[gene] = True
        if context != "MUTATION":
            if gene in list(context_df.index): context_data_available[gene] = True
        else:
            if gene in set([mut.split(" p.")[0] for mut in context_df["Mutation"]]):
                context_data_available[gene] = True

    # Write not available genes
    not_ess_data_available = [gene for gene, availability in ess_data_available.items() if availability == False]
    not_context_data_available = [gene for gene, availability in context_data_available.items() if
                                  availability == False]
    not_data_available = set(not_ess_data_available).union(set(not_context_data_available))

    if len(not_ess_data_available) != 0: print("The essentiality data for below query gene(s) is not found:")
    for gene in not_ess_data_available: print(gene)
    print("\n")

    if len(not_context_data_available) != 0: print("The context data for below query gene(s) is not found:")
    for gene in not_context_data_available: print(gene)
    print("\n")

    # If there is no any available query genes, system exit should be given!
    if len(not_ess_data_available) == len(args["GENES"]): sys.exit("There is no available query gene to analyse!")
    if len(not_context_data_available) == len(args["GENES"]): sys.exit("There is no available query gene to analyse!")

    if len(not_data_available) != len(args["GENES"]):
        query_genes = list(set(args["GENES"]).difference(not_data_available))

else: sys.exit("At least one query gene should be given!")


# Only check if the co-gene(s) is in the essentiality data!

co_genes = []

if args["ADD_DEPENDANT"]:

    # Used for only essentiality data with this(ese) gene(s)
    if args["CoGENES"]:

        co_ess_data_available = {}

        for cogene in args["CoGENES"]:
            if cogene in list(essentiality_df.index): co_ess_data_available[cogene] = True

        # Write not available co-genes
        not_co_ess_data_available = [cogene for cogene, availability in co_ess_data_available.items()
                                     if availability == False]

        if len(not_co_ess_data_available) != 0: print("The essentiality data for below co-gene(s) is not found:")
        for cogene in not_co_ess_data_available: print(cogene)
        print("\n")

        # If there is no any available query genes, all the genes in the essentiality table should be used!

        if len(not_co_ess_data_available) != len(args["CoGENES"]):
            co_genes = list(set(args["CoGENES"]).difference(not_co_ess_data_available))

    else:
        sys.exit("If -add_dependant was added, then user needs to specify -cogene or -cofile!")

else:
    print("No dependant gene was specified! All genes in the essentiality table will be searched!")
    co_genes = set(essentiality_df.index)

########################################################################################################################
########################################################################################################################
# CELL LINE SEPARATION  #

site_given, sub_site_given = False, False
new_discrete_data = False
separation = False  # Input for separate_CLs's user_defined_separation
interested, other, group, sites = False, False, None, None
other_cls = []

if args["SEPARATION"] and args["SITE"]:
    args["SITE"] = None

if args["SEPARATION"]:
    if args["INTERESTED"] and args["OTHER"]:
        separation = True
        group = "USER"
        interested, other = True, True
    elif args["INTERESTED"]:
        separation = True
        group = "USER"
        interested, other = True, False
        other_cls = list(set(essentiality_df.columns).difference(set(args["INTERESTED"])))
    else: sys.exit("If -separate was added, then user needs to specify -interested and -other!")

elif args["SITE"]:
    site_given = True
    group = args["SITE"]
    if args["GROUPS"]:
        sub_site_given = True
        sites = [g for g in args["GROUPS"]]
    else: print("All %s types will be used!\n" %group)

elif context_str == "discrete":
    new_discrete_data = True
    group = "USER"

else: sys.exit("No cell line separation information was found!")


########################################################################################################################
########################################################################################################################
# LIMITATIONS #

cor_limit, significancy = None, None

if args["CORLIMIT"]: cor_limit = args["CORLIMIT"]
if args["PLIMIT"]: significancy = args["PLIMIT"]


########################################################################################################################
########################################################################################################################
# FUNCTIONS #

def separate_CLs(project_object, site_property, site = None, gene = None,
                 hotspot = False, oncogenic = False, amplified = False, depleted = False,
                 user_defined_separation = False, new_discrete_data = False, pre_defined = True):
    """
    Separation of Cell lines according to their tissue of origin
    :param project_object: BROAD or SANGER (the project that is used)
    :param site: Interested tissue/cancer
    :param site_property: If site is a "tissue" name or a "cancer" name.
    :param gene: The gene having mutation that we need to separate the cell lines in accordance.
    :param hotspot: This variable can be used only with the gene and mutation variable.
                    If True only the hotspot mutation or hotspot carrying genes are taking into account!
    :param oncogenic: This variable can be used only with the gene and mutation variable.
                    If True only the oncogenic mutation or hotspot carrying genes are taking into account!
    :param amplified: If the given gene have any amplification
    :param depleted: If the given gene have any depletion
    :param user_defined_separation: If user wants to specify separated lists of cell lines.
    :return: The cell line dictionary for interested site and cell lines, cell line dictionary for
             interested site (or all) but separated inside the group such as by mutation
    """

    grouped_CLs, two_grouped_CLs = {}, {}
    all_CLs, normal_CLs = [], []

    for i, l in project_object.items():

        # Which SITE the user specified!

        project_site = l.tissue if site_property == "tissue" else (
            l.cancer if site_property == "cancer" and l.cancer != None else None)

        if project_site != False and project_site != None and type(project_site) != float:
            if project_site.upper() != "NORMAL":
                if i not in all_CLs: all_CLs.append(i)
                if project_site.upper() not in grouped_CLs.keys():
                    grouped_CLs[project_site.upper()] = [i]
                else:
                    grouped_CLs[project_site.upper()].append(i)
            else:
                if i not in normal_CLs: normal_CLs.append(i)

    if pre_defined and user_defined_separation == False and new_discrete_data == False:
        if site != None:
            CLs = grouped_CLs[site.upper()]
            non_CLs = [CL for CL in all_CLs if CL not in CLs]

            if gene != None and amplified == False and depleted == False:

                mutated_CLs, non_mutated_CLs = [], []

                # We will focus on only CLs because there is site and mutated gene are specified!

                if hotspot == False and oncogenic == True:
                    mutated_CLs = [cl for cl in CLs if project_object[cl].oncogenic_mutated_genes != None and
                                   gene in project_object[cl].oncogenic_mutated_genes]
                    non_mutated_CLs = [cl for cl in CLs if cl not in mutated_CLs]

                    two_grouped_CLs = {"interested": mutated_CLs, "other": non_mutated_CLs}

                elif hotspot == True and oncogenic == False:

                    mutation_info_CLs = [cl for cl in CLs if project_object[cl].hotspot_mutations != None]

                    mutated_CLs = [cl for cl in mutation_info_CLs if gene in
                                   set([m.split(" p.")[0] for m in project_object[cl].hotspot_mutations])]

                    non_mutated_CLs = [cl for cl in CLs if cl not in mutated_CLs]

                    two_grouped_CLs = {"interested": mutated_CLs, "other": non_mutated_CLs}

            if amplified == True and gene != None:

                amplified_CLs, non_amplified_CLs = [], []

                # We will focus in CLs because there is site and amplification were specified

                amplified_CLs = [cl for cl in CLs if project_object[cl].CNV and gene in project_object[cl].cnv.index and
                                 project_object[cl].cnv.loc[gene] >= 1.5]

                non_amplified_CLs = [cl for cl in CLs if cl not in amplified_CLs]

                two_grouped_CLs = {"interested": amplified_CLs, "other": non_amplified_CLs}

            if depleted == True and gene != None:

                depleted_CLs, non_depleted_CLs = [], []

                # We will focus in CLs because there is site and depletion were specified

                depleted_CLs = [cl for cl in CLs if project_object[cl].CNV and gene in project_object[cl].cnv.index and
                                project_object[cl].cnv.loc[gene] <= -1.5]

                non_depleted_CLs = [cl for cl in CLs if cl not in depleted_CLs]

                two_grouped_CLs = {"interested": depleted_CLs, "other": non_depleted_CLs}

        else:

            # There is no tissue or cancer type was specified!

            if gene != None and amplified == False and depleted == False:

                mutated_CLs, non_mutated_CLs = [], []

                if hotspot == False and oncogenic == True:

                    mutated_CLs = [cl for cl in all_CLs if project_object[cl].oncogenic_mutated_genes != None and
                                   gene in project_object[cl].oncogenic_mutated_genes]

                    non_mutated_CLs = [cl for cl in all_CLs if cl not in mutated_CLs]

                    two_grouped_CLs = {"interested": mutated_CLs, "other": non_mutated_CLs}

                elif hotspot == True and oncogenic == False:

                    mutation_info_CLs = [cl for cl in all_CLs if project_object[cl].hotspot_mutations != None]

                    mutated_CLs = [cl for cl in mutation_info_CLs if gene in
                                   set([m.split(" p.")[0] for m in project_object[cl].hotspot_mutations])]

                    non_mutated_CLs = [cl for cl in all_CLs if cl not in mutated_CLs]

                    two_grouped_CLs = {"interested": mutated_CLs, "other": non_mutated_CLs}

                elif hotspot == True and oncogenic == False:

                    mutation_info_CLs = [cl for cl in all_CLs if project_object[cl].hotspot_mutations != None]

                    mutated_CLs = [cl for cl in mutation_info_CLs if mutation in project_object[cl].hotspot_mutations]

                    non_mutated_CLs = [cl for cl in all_CLs if cl not in mutated_CLs]

                    two_grouped_CLs = {"interested": mutated_CLs, "other": non_mutated_CLs}


            if amplified == True and gene != None:

                amplified_CLs, non_amplified_CLs = [], []

                amplified_CLs = [cl for cl in all_CLs if project_object[cl].CNV and gene in project_object[cl].cnv.index and
                                 project_object[cl].cnv.loc[gene] >= 1.5]
                non_amplified_CLs = [cl for cl in all_CLs if cl not in amplified_CLs]

                two_grouped_CLs = {"interested": amplified_CLs, "other": non_amplified_CLs}

            if depleted == True and gene != None:

                depleted_CLs, non_depleted_CLs = [], []

                depleted_CLs = [cl for cl in all_CLs if project_object[cl].CNV and gene in project_object[cl].cnv.index and
                                project_object[cl].cnv.loc[gene] <= -1.5]
                non_depleted_CLs = [cl for cl in all_CLs if cl not in depleted_CLs]
                two_grouped_CLs = {"interested": depleted_CLs, "other": non_depleted_CLs}

    elif user_defined_separation and new_discrete_data == False and pre_defined == False:
        global other
        if other: two_grouped_CLs = {"interested" : args["INTERESTED"], "other": args["OTHER"]}
        else:
            global other_cls
            two_grouped_CLs = {"interested" : args["INTERESTED"],
                               "other": other_cls}

    elif new_discrete_data and user_defined_separation == False and pre_defined == False:

        cl_df = pandas.DataFrame(new_discrete_data.loc[gene])
        two_grouped_CLs = {"interested": list(cl_df[cl_df[gene] == True].index),
                           "other": list(cl_df[cl_df[gene] == False].index)}

    return grouped_CLs, two_grouped_CLs


def difference_between_multiple_groups(grouped_CL_dict, df, gene_y, gene_x,
                                       group_by, searched_context, significancy = 0.05):
    """
    Kruskal-Wallis Test between multiple categorical and 1 continuous groups
    :param grouped_CL_dict: The dictionary of group of cell lines in different tissue/cancer.
    :param df: DataFrame having values for interested property.
    :param gene_y: Intersted Gene that are measured.
    :param gene_x: The gene change the classification of cell lines or same as gene_y.
    :param group_by : How we are separating the group according to "TISSUE" or "CANCER".
    :param searched_context : The context we are searching on.
    :param significancy: The maximum P value for the test. Default is 0.05.
    :return: group_value_df (Data frame for values and corresponding tissue/cancer,
            significant_groups (The groups having significantly higher value for the interested property.)
    """

    significant_groups = {}
    group_value_df = None

    # LIST OF VALUES

    if gene_x in list(df.index) and gene_y in list(df.index):
        value_dfs = []
        for group, cl_list in grouped_CL_dict.items():
            cl_values = list(df.loc[gene_y, [cl for cl in cl_list if cl in df.columns]])

            if len(cl_values) >= 3: value_dfs.append(pandas.DataFrame({group: cl_values}))

        group_value_df = pandas.melt(pandas.concat(value_dfs, axis=1))
        group_value_df = group_value_df.dropna()
        group_value_df.columns = ["Group", "Value"]

        # KRUSKAL-WALLIS TEST BETWEEN GROUPS

        if len(set(list(group_value_df.Value))) != 1:

            _, kruskal_p_val = stats.kruskal(*[list(i[1].Value) for i in group_value_df.groupby(["Group"])])

            if kruskal_p_val < significancy:

                # NOT POST-HOC BUT ONE BY ONE MANN WHITNEY U TEST

                for group in set(group_value_df.Group):
                    data = [list(group_value_df[group_value_df.Group == group].Value),
                            list(group_value_df[group_value_df.Group != group].Value)]

                    if sum(data[0]) != 0 and sum(data[1]) != 0:

                         _, mann_p_val = mannwhitneyu(data[0], data[1])

                         if mann_p_val < significancy: significant_groups[group] = {"MWU_P": mann_p_val, "DATA": data}

    return group_value_df, significant_groups


def difference_between_groups(grouped_CL_dict, df, gene_y, gene_x,
                              group_by, searched_context, significancy = 0.05):
    """
    Mann Whitney Test between 1 categorical and 1 continuous groups
    :param grouped_CL_dict: The dictionary of group of cell lines in different tissue/cancer/cancer subtype.
    :param df: DataFrame having values for interested property.
    :param gene_y: Intersted Gene that are measured.
    :param gene_x: The gene change the classification of cell lines or same as gene_y.
    :param group_by : How we are separating the group according to "TISSUE" or "CANCER", "CANCER SUBTYPE".
    :param searched_context : The context we are searching on.
    :param significancy: The maximum P value for the test. Default is 0.05.
    :return: significance (boolean value if the result significant or not), p (p value of the test), data
    """

    # SITE

    significance, mann_p_val, data = None, None, None

    # LIST OF VALUES

    if gene_x in list(df.index) and gene_y in list(df.index):

        cl_values = list(df.loc[gene_y, [CL for CL in grouped_CL_dict["interested"] if CL in df.columns]])
        other_cl_values = list(df.loc[gene_y, [CL for CL in grouped_CL_dict["other"] if CL in df.columns]])

        n = 0
        group_value_df = pandas.DataFrame(0, index=range(len(cl_values) + len(other_cl_values)),
                                          columns=["Group", "Value"])

        if len(cl_values) >= 3 and len(other_cl_values) >=3:
            for k in cl_values:
                group_value_df.loc[n, "Group"] = "interested"
                group_value_df.loc[n, "Value"] = k
                n += 1
            for k in other_cl_values:
                group_value_df.loc[n, "Group"] = "other"
                group_value_df.loc[n, "Value"] = k
                n += 1

        # MANN WHITNEY U TEST BETWEEN TWO GROUPS

        data = [list(group_value_df[group_value_df.Group == "interested"].Value),
                list(group_value_df[group_value_df.Group == "other"].Value)]

        if len(data[0]) >= 3 and len(data[1]) >= 3 and data[0] != data[1]:

            _, mann_p_val = mannwhitneyu(data[0], data[1])

            if mann_p_val < significancy: significance = True

    return significance, mann_p_val, data


def association_between_groups(grouped_CL_dict, gene_x, gene_y, df_x, df_y, site,
                               group_by, searched_context, correlation_limit = None,
                               x_by=None, y_by=None, significancy=0.05):
    """
    Correlation and Regression analysis between 2 continuous groups
    :param grouped_CL_dict: The dictionary of group of cell lines in different tissue/cancer
    :param gene_x : The gene whose information will be taken from df_x
    :param gene_y : The gene whose information will be taken from df_y
    :param df_x : LMPlot x axis property (if essentiality will be used, use in x axis)
    :param df_y : LMPlot y axis property (if expression will be used, use in y axis)
    :param x_by : Which property will be represented "EXPRESSION", "ESSENTIALITY", "CNVs" in df_x
    :param y_by : Which property will be represented "EXPRESSION", "ESSENTIALITY", "CNVs" in df_y
    :param site : The name of the interested tissue/cancer.
    :param group_by : How we are separating the group according to "TISSUE" or "MUTATION"
    :param searched_context : The context we are searching on
    :param significancy: The maximum P value for the test. Default is 0.05.
    :param correlation_limit: The minimum absolute correlation coefficient.
    :return: True or None for association, correlation coefficient, p value, value data frame
    """

    # SITE

    association, cor_coef, cor_p, df = None, None, None, None

    if site is not None:

        cls = [cl for cl in grouped_CL_dict[site.upper()] if cl in df_x.columns and cl in df_y.columns]

    else:
        cls = list(set(df_x.columns).intersection(set(df_y.columns)))

    if len(cls) >= 5:

        if gene_x in list(df_x.index) and gene_y in list(df_y.index):

            if type(df_y.loc[gene_y]) == pandas.Series and type(df_x.loc[gene_x]) == pandas.Series:

                df = pandas.concat([pandas.DataFrame(df_x.loc[gene_x, cls]),
                                    pandas.DataFrame(df_y.loc[gene_y, cls])], axis=1)

                df = df.dropna()

                if len(df.columns) == 2 and len(df.index) >= 5:
                    df.columns = [x_by, y_by]

                    # Correlation between properties

                    cor_coef, cor_p = stats.pearsonr(df[x_by], df[y_by])

                    if cor_p < significancy:
                        if correlation_limit != None and numpy.abs(cor_coef) >=correlation_limit: association = True
                        elif correlation_limit == None: association = True

    return association, cor_coef, cor_p, df

def store_network_object(network_obj):

    network_df = pandas.DataFrame(0, index=list(range(len(network_obj.edges()))),
                                  columns=["from", "to", "pvalue", "median_interested", "median_other",
                                           "effector", "affected", "interested_cl", "other_cl",
                                           "site", "association", "to_Who", "group", "project"])
    i = 0
    for edge in network_obj.edges:
        network_df.loc[i, "from"] = edge[0]
        network_df.loc[i, "to"] = edge[1]
        network_df.loc[i, "pvalue"] = network_obj.get_edge_data(edge[0], edge[1], edge[2])["p_value"]
        network_df.loc[i, "median_interested"] = network_obj.get_edge_data(edge[0], edge[1], edge[2])["median_interested"]
        network_df.loc[i, "median_other"] = network_obj.get_edge_data(edge[0], edge[1], edge[2])["median_other"]
        network_df.loc[i, "effector"] = network_obj.get_edge_data(edge[0], edge[1], edge[2])["effector"]
        network_df.loc[i, "affected"] = network_obj.get_edge_data(edge[0], edge[1], edge[2])["affected"]
        network_df.loc[i, "interested_cl"] = network_obj.get_edge_data(edge[0], edge[1], edge[2])["interested_cl"]
        network_df.loc[i, "other_cl"] = network_obj.get_edge_data(edge[0], edge[1], edge[2])["other_cl"]
        network_df.loc[i, "site"] = network_obj.get_edge_data(edge[0], edge[1], edge[2])["site"]
        network_df.loc[i, "association"] = network_obj.get_edge_data(edge[0], edge[1], edge[2])["association"]
        network_df.loc[i, "to_Who"] = network_obj.get_edge_data(edge[0], edge[1], edge[2])["to_who"]
        network_df.loc[i, 'group'] = network_obj.get_edge_data(edge[0], edge[1], edge[2])["group"]
        network_df.loc[i, "project"] = network_obj.get_edge_data(edge[0], edge[1], edge[2])["project"]
        i += 1
    return network_df


########################################################################################################################
########################################################################################################################
# FLOWS #

def flow_continuous_pre_defined_CLs():

    global query_genes, co_genes, obj, essentiality_df, context, context_df, group, significancy, cor_limit, \
        site_given, sub_site_given, new_discrete_data, separation, interested, other, group, sites, other_cls

    cen_network = networkx.MultiDiGraph()

    for query_gene in query_genes:
        if co_genes:
            if query_gene not in co_genes:
                co_genes.append(query_gene)
            affected_genes = co_genes
        else:
            affected_genes = list(essentiality_df.index)

        # Group selected - Not user-defined
        # Multiple group difference - continuous data / EDGE 1

        whole_group_CLs, _ = separate_CLs(project_object=obj, site_property=group)

        whole_group_value_df1, whole_significancy_dict1 = difference_between_multiple_groups(
            grouped_CL_dict= whole_group_CLs, df = context_df, gene_y = query_gene, gene_x = query_gene,
            group_by = group, searched_context = group, significancy=significancy)

        if not sub_site_given:

            for significant_group1, d1 in whole_significancy_dict1.items():

                # NETWORK EDGE 1 (LINKAGE)

                cen_network.add_edge(query_gene, significant_group1, p_value=d1["MWU_P"],
                                     median_interested = numpy.median(d1["DATA"][0]),
                                     median_other = numpy.median(d1["DATA"][1]), effector=group,
                                     affected=context.lower(), site = significant_group1,
                                     interested_cl = len(d1["DATA"][0]), other_cl = len(d1["DATA"][1]),
                                     to_who="PANCANCER", association=None, group = group, project = project)

        else:
            for sub_group in sites:
                sub_group = sub_group.upper()
                if sub_group in whole_significancy_dict1.keys():

                    # NETWORK EDGE 1 (LINKAGE)
                    d1 = whole_significancy_dict1[sub_group]

                    cen_network.add_edge(query_gene, sub_group, p_value=d1["MWU_P"],
                                         median_interested=numpy.median(d1["DATA"][0]),
                                         median_other=numpy.median(d1["DATA"][1]), effector=group,
                                         affected=context.lower(), site=sub_group,
                                         interested_cl=len(d1["DATA"][0]), other_cl=len(d1["DATA"][1]),
                                         to_who="PANCANCER", association=None, group=group, project=project)

        # Multiple group difference - Essentiality / EDGE 2

        whole_group_value_df2, whole_significancy_dict2 = difference_between_multiple_groups(
            grouped_CL_dict= whole_group_CLs, df = essentiality_df, gene_y = query_gene, gene_x = query_gene,
            group_by = group, searched_context = group, significancy=significancy)

        if not sub_site_given:
            for significant_group2, d2 in whole_significancy_dict2.items():

                # NETWORK EDGE 2 (LINKAGE)
                cen_network.add_edge(query_gene, significant_group2, p_value=d2["MWU_P"],
                                     median_interested = numpy.median(d2["DATA"][0]),
                                     median_other = numpy.median(d2["DATA"][1]), effector=group,
                                     affected="essentiality", site= significant_group2,
                                     interested_cl = len(d2["DATA"][0]), other_cl = len(d2["DATA"][1]),
                                     to_who="PANCANCER", association=None, group = group, project = project)

        else:
            for sub_group in sites:
                sub_group = sub_group.upper()
                if sub_group in whole_significancy_dict2.keys():

                    # NETWORK EDGE 2 (LINKAGE)
                    d2 = whole_significancy_dict2[sub_group]

                    cen_network.add_edge(query_gene, sub_group, p_value=d2["MWU_P"],
                                         median_interested=numpy.median(d2["DATA"][0]),
                                         median_other=numpy.median(d2["DATA"][1]), effector=group,
                                         affected="essentiality", site=sub_group, interested_cl=len(d2["DATA"][0]),
                                         other_cl=len(d2["DATA"][1]), to_who="PANCANCER", association=None,
                                         group=group, project=project)

        for co_gene in affected_genes:

            # Context Data to Essentiality Association - No site specific / EDGE 3

            data_association, data_cor_coef, data_cor_p, data_cor_df = association_between_groups(
                grouped_CL_dict= whole_group_CLs, gene_x = co_gene, gene_y = query_gene, df_x = essentiality_df,
                df_y = context_df, site = None, group_by = None, searched_context = context,
                correlation_limit = cor_limit, significancy=significancy, x_by = "ESSENTIALITY", y_by = context)

            if data_association is not None and data_cor_p < significancy:

                # NETWORK EDGE 3 (DATA-ESS ASSOCIATION)
                cen_network.add_edge(query_gene, co_gene, p_value=data_cor_p, median_interested=None,
                                     median_other = None, effector=context.lower(), affected="essentiality",
                                     association=data_cor_coef, interested_cl = len(list(data_cor_df.index)),
                                     other_cl = len(list(data_cor_df.index)), to_who="PANCANCER",
                                     site = "PANCANCER", group = group,project = project)

            if not sub_site_given:

                # Context Data to Essentiality Association - Site specific / EDGE 4

                for g in whole_group_CLs.keys():

                    data_inside_association, data_inside_cor_coef, data_inside_cor_p, \
                    data_inside_cor_df = association_between_groups(
                        grouped_CL_dict=whole_group_CLs, gene_x=co_gene, gene_y=query_gene, df_x=essentiality_df,
                        df_y=context_df, site=g, group_by=group, searched_context=context,
                        correlation_limit=cor_limit, significancy=significancy, x_by="ESSENTIALITY", y_by=context)

                    if data_inside_association != None and data_inside_cor_p < significancy:

                        # NETWORK EDGE 4 (DATA-ESS ASSOCIATION IN GROUP)
                        cen_network.add_edge(query_gene, co_gene, p_value=data_inside_cor_p, median_interested=None,
                                             median_other=None,interested_cl = len(list(data_inside_cor_df.index)),
                                             other_cl = len(list(data_inside_cor_df.index)), effector=context.lower(),
                                             affected="essentiality", association=data_inside_cor_coef,
                                             to_who=g, site=g, group=group, project=project)


            else:

                whole_group_CLs, _ = separate_CLs(project_object=obj, site_property=group)

                for sub_group in sites:
                    sub_group = sub_group.upper()

                    data_inside_association, data_inside_cor_coef, data_inside_cor_p, \
                    data_inside_cor_df = association_between_groups(
                        grouped_CL_dict=whole_group_CLs, gene_x=co_gene, gene_y=query_gene, df_x=essentiality_df,
                        df_y=context_df, site=sub_group.upper(), group_by=group, searched_context=context,
                        correlation_limit=cor_limit, significancy=significancy, x_by="ESSENTIALITY", y_by=context)

                    if data_inside_association is not None and data_inside_cor_p < significancy:
                        # NETWORK EDGE 4 (DATA-ESS ASSOCIATION IN GROUP)
                        cen_network.add_edge(query_gene, co_gene, p_value=data_inside_cor_p, median_interested=None,
                                             median_other=None, interested_cl=len(list(data_inside_cor_df.index)),
                                             other_cl=len(list(data_inside_cor_df.index)), effector=context.lower(),
                                             affected="essentiality", association=data_inside_cor_coef,
                                             to_who=sub_group.upper(), site=sub_group.upper(), group=group,
                                             project=project)

    return cen_network


def flow_continuous_user_defined_CLs():

    global query_genes, co_genes, obj, essentiality_df, context, context_df, group, significancy, cor_limit, \
        site_given, sub_site_given, new_discrete_data, separation, interested, other, group, sites, other_cls

    cen_network = networkx.MultiDiGraph()

    affected_genes = []
    for query_gene in query_genes:
        if co_genes:
            if query_gene not in co_genes:
                co_genes.append(query_gene)
            affected_genes = co_genes
        else: affected_genes = list(essentiality_df.index)

        for co_gene in affected_genes:

            # User-defined Groups

            _, user_group = separate_CLs(user_defined_separation = separation, pre_defined = site_given)

            # Context Data to Essentiality Association - No site specific / EDGE 1

            data_association, data_cor_coef, data_cor_p, data_cor_df = association_between_groups(
                grouped_CL_dict= user_group, gene_x = co_gene, gene_y = query_gene, df_x = essentiality_df,
                df_y = context_df, site = None, group_by = None, searched_context = context.lower(),
                correlation_limit = cor_limit, significancy=significancy, x_by = "ESSENTIALITY", y_by = context)

            if data_association is not None and data_cor_p < significancy:

                # NETWORK EDGE 2 (DATA-ESS ASSOCIATION)
                cen_network.add_edge(query_gene, co_gene, p_value=data_cor_p, median_interested=None, median_other = None,
                                     effector=context.lower(), affected="essentiality", association=data_cor_coef,
                                     interested_cl = len(list(data_cor_df.index)), other_cl = None, to_who="OTHER",
                                     site = group, group = group, project = project)

    return cen_network


def flow_discrete():

    global query_genes, co_genes, obj, essentiality_df, context, context_df, group, significancy, cor_limit, \
        site_given, sub_site_given, new_discrete_data, separation, interested, other, group, sites, other_cls,\
        MUT, AMP, DEP, CNT, COMPESS, HOTSPOT, ONCOGENE

    cen_network = networkx.MultiDiGraph()

    affected_genes = []
    for query_gene in query_genes:
        if co_genes:
            if query_gene not in co_genes:
                co_genes.append(query_gene)
            affected_genes = co_genes
        else: affected_genes = list(essentiality_df.index)

        for co_gene in affected_genes:

            # Not a specific group
            whole_grouped_CLs, whole_inside_grouped_CLs = separate_CLs(
                project_object=obj, site_property=group, site=None, gene=query_gene,
                hotspot=HOTSPOT, oncogenic=ONCOGENE, amplified = AMP, depleted = DEP,
                user_defined_separation = separation, new_discrete_data = CNT, pre_defined = site_given)

            whole_significance, whole_mann_p_val, whole_data = difference_between_groups(
                grouped_CL_dict=whole_inside_grouped_CLs, df=essentiality_df, gene_y=co_gene, gene_x=query_gene,
                group_by=context, searched_context=context, significancy=significancy)

            if whole_significance:

                # NETWORK EDGE - EDGE 1 - not site specific - PANCANCER
                cen_network.add_edge(
                    query_gene, co_gene, p_value=whole_mann_p_val, median_interested=numpy.median(whole_data[0]),
                    median_other = numpy.median(whole_data[1]), effector=context.lower(), affected="essentiality",
                    interested_cl = len(whole_data[0]), other_cl = len(whole_data[1]),association=None,
                    to_who="PANCANCER", site="PANCANCER", group = group, project = project)

            if not sub_site_given:

                # Group by group

                for g in whole_grouped_CLs.keys():

                    sub_grouped_CLs, sub_inside_grouped_CLs = separate_CLs(
                        project_object=obj, site_property=group, site=g, gene = query_gene,
                        hotspot=HOTSPOT, oncogenic=ONCOGENE, amplified=AMP, depleted=DEP,
                        user_defined_separation=separation, new_discrete_data=CNT, pre_defined=site_given)

                    sub_significance, sub_mann_p_val, sub_data = \
                        difference_between_groups(
                            grouped_CL_dict=sub_inside_grouped_CLs, df=essentiality_df, gene_y=co_gene,
                            gene_x=query_gene, group_by=context, searched_context=context, significancy=significancy)

                    if sub_significance:
                        if query_gene == co_gene:

                            # NETWORK EDGE - EDGE 2 - inside group - gene's effect itself on site
                            cen_network.add_edge(
                                query_gene, g, p_value=sub_mann_p_val, median_interested=numpy.median(sub_data[0]),
                                median_other = numpy.median(sub_data[1]), effector=context.lower(),
                                affected="essentiality", interested_cl = len(sub_data[0]), other_cl = len(sub_data[1]),
                                association=None, to_who=g, site=g, group = group, project = project)
                        else:

                            # NETWORK EDGE - EDGE 3 - inside group - gene's effect on another on site
                            cen_network.add_edge(
                                query_gene, co_gene, p_value=sub_mann_p_val,
                                median_interested=numpy.median(sub_data[0]), median_other=numpy.median(sub_data[1]),
                                effector=context.lower(), affected="essentiality", association=None,
                                interested_cl=len(sub_data[0]), other_cl=len(sub_data[1]),
                                to_who=g, site=g, group=group, project=project)

            else:

                for g in sites:

                    sub_grouped_CLs, sub_inside_grouped_CLs = separate_CLs(
                        project_object=obj, site_property=group, site=g, gene=query_gene,
                        hotspot=HOTSPOT, oncogenic=ONCOGENE, amplified=AMP, depleted=DEP,
                        user_defined_separation=separation, new_discrete_data=CNT, pre_defined=site_given)

                    sub_significance, sub_mann_p_val, sub_data = \
                        difference_between_groups(
                            grouped_CL_dict=sub_inside_grouped_CLs, df=essentiality_df, gene_y=co_gene,
                            gene_x=query_gene, group_by=context, searched_context=context, significancy=significancy)

                    if sub_significance:
                        if query_gene == co_gene:

                            # NETWORK EDGE - EDGE 2 - inside group - gene's effect itself on site

                            cen_network.add_edge(
                                query_gene, g, p_value=sub_mann_p_val, median_interested=numpy.median(sub_data[0]),
                                median_other=numpy.median(sub_data[1]), effector=context.lower(),
                                affected="essentiality", interested_cl=len(sub_data[0]), other_cl=len(sub_data[1]),
                                association=None, to_who=g, site=g, group=group, project=project)
                        else:

                            # NETWORK EDGE - EDGE 3 - inside group - gene's effect on another on site

                            cen_network.add_edge(
                                query_gene, co_gene, p_value=sub_mann_p_val,
                                median_interested=numpy.median(sub_data[0]), median_other=numpy.median(sub_data[1]),
                                effector=context.lower(), affected="essentiality", association=None,
                                interested_cl=len(sub_data[0]), other_cl=len(sub_data[1]),
                                to_who=g, site=g, group=group, project=project)
    return cen_network


########################################################################################################################
########################################################################################################################
# MAIN #

def main():

    global path, result_path, context_str, site_given, separation, context, group

    if context_str == "continuous" and site_given:

        cen_object = flow_continuous_pre_defined_CLs()

    if context_str == "continuous" and site_given == False:

        cen_object = flow_continuous_user_defined_CLs()

    if context_str == "discrete":

        cen_object = flow_discrete()


    print("The CEN object was created.\n")
    print("Making it as a csv file...\n")

    df = store_network_object(cen_object)

    df.to_csv(result_path + args["OUTPUT_FILE"] + ".csv")

    print("The CEN was written in %s\n" %result_path)


main()

