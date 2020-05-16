
# CEN-tools

**PyCEN** is the *Python* implementation of the **CEN-tools** web interface ([http://cen-tools.com/](http://cen-tools.com/)) to analyse contexts in pre-defined or user-defined data in detail.

## How to use ?
**1. Download the project.**

	> cd ~/
	> git clone https://github.com/rraadd88/centools_pycen.git

**2. Go to the folder.**

	> cd centools_pycen/

**3. Install Required Packages.**

	> pip3 install -r requirements.txt -U

**4. Download data.**     

Download two directories in the centools_pycen folder.     

1. [data](https://gitlab.ebi.ac.uk/petsalakilab/centools_pycen/tree/master/data)    
2. [test_data](https://gitlab.ebi.ac.uk/petsalakilab/centools_pycen/tree/master/test_data)    


**2. Prepare you input arguments:**

**GENES** : One of the arguments **-gene** or **-genefile** *should* be specified. In here the user has an opportunity to give only one interested gene (**-gene**) or a list of genes as a ".csv" file (**-genefile**).

**PROJECT** : CEN-tools uses two projects as *SANGER* and *BROAD*. The user can select one of them for the further analysis with **-project**. The data will be used for essentiality, context, and cell line information for separation unless user data is defined!

**ESSENTIALITY** :  Unless specified, the essentiality data of the selected project will be used. However, the user can also specify their own essentiality data with **-add_essentiality** option. If the user choose to give their own essentiality data, they need to provide the filepath with **-essfile**. 
	*Important note: If the pre-defined context is to be used, the new essentiality data should have the cell line names matching with the context data. Otherwise, the new context data should also match with the new essentiality data!* 

**CONTEXT** : Moreover, the user can specify a new data for *context* with **-add_context** or pre-defined contexts of the selected project will be used! If **-add_context** is given, then the user *should* provide both the data (**-contextfile**) and the data structure (**-contextstr**) as *continuous* or *discrete*.    The *contextfile* should be consists of rows as genes, columns as cell lines, and entities as Boolean values (True/False). *Such as for fusion data, if a cell line (column) has a fusion version of the gene in the corresponding row, entity can be True.* The *continuous* data will be used for correlation, the *discrete* data will be used for separation and then analysis will be done for searching associations. On the other hand, if the one of the pre-defined contexts will be used, then corresponding context should be provided with one of them shown below:
		 EXPRESSION : **-exp**
		 MUTATION : **-mut**
		 AMPLIFICATION : **-amp**
		 DEPLETION : **-dep**
		 COMPARE_ESSENTIALITY: **-compess** (This is a special context option for comparison of essentiality profiles of query gene(s) between user selected cell lines with **-separate**. -example below as *WRN* essentiality with MSI status of cell lines.)

**DEPENDANTS** : Using the default settings, PyCEN analyses all dependences between the query gene/genes; however, the user can also specify the dependant gene/genes to reduce the time for execution of the script. To specify an interested dependance, the user firstly should provide **-add_dependent**, then give one gene with **-cogene** or a list of gene (.csv) with **-cogenefile** .

**GROUPING** : The cell lines can be grouped according to their tissues or cancer types, new discrete context data, or user specified separated cell line files. The user can use *tissue* and *cancer* separation by using **-group** as *tissue* or *cancer*.  Moreover, they can specify a specific tissue/cancer by using **-groupname** or several tissue/cancer types by **-groupfile** with a ".csv" file. On the other hand, after using **Cell line selector**, or manually created lists of cell lines, the user can give their lists to the analysis. To do this, the user first needs to specify **-separate**, then submit two ".csv" files for interested **-interested_group** and other **other_group** groups for comparison. If the **-other_group** will not be given, the rest of the cell lines in the corresponding project will be used. If a new **discrete context data** will be given, then cell lines will be grouped according to the data.   
	*Important note* : Be sure the given tissue/cancer type name in **-separate** option is as the format in CEN-tools. All usable tissue and cancer names can be found in *test_data* folder.

**THRESHOLDS** : For the analysis, the user *can* specify the thresholds for the statistical tests with **-limcor** for absolute value of *minimum correlation coefficient* and **-limp** for *maximum p value* for significance. 

**OUTPUT** : The user can specify the output path with **-o**. Otherwise, the current directory will be used!

**OUTPUT FILE** : The user can specify the output file name without extension, otherwise "CEN" will be used as name. *Important note: When the user run PyCEN script more than one and not give a specific output file name, the newer result can overwrite the existing one.* 

**MUTATION** : If the user wants to use the pre-defined **-mut** context, then the separation type should be specified with one of the **-hotspot**, **-oncogenic**. These options separate the cell lines based on the presence of any mutation on the hotspot or oncogenic mutation carrying query genes.  
	
**4. Examples:**

*BRAF* CEN in the context of *hotspot* *mutation* of *BRAF* gene  in all *tissue* types from the *BROAD* project data with *BRAF_Mutation_Tissue_BROAD* result file name:	

`python3 pyCEN.py -gene BRAF -mut -add_dependant -cogene BRAF -group tissue -project BROAD -hotspot -ofile BRAF_Mutation_Tissue_BROAD`

*SOX10* CEN in the context of *expression* of *SOX10* gene in all *tissue* types with the *BROAD* project data with *SOX10_Expression_Tissue_BROAD* result file name: 

`python3 pyCEN.py -gene SOX10 -exp -add_dependant -cogene SOX10 -group tissue -project BROAD -ofile SOX10_Expression_Tissue_BROAD`

*MAPK1* CEN in the context of *oncogenic* *mutation* of *BRAF* gene in *skin  cancer* type with the *BROAD* project data with *BRAF_MAPK1_Mutation_Cancer_Skin_BROAD* result file name:

`python3 pyCEN.py -gene BRAF -mut -add_dependant -cogene MAPK1 -group cancer -groupname Skin\ cancer -project BROAD -oncogenic -ofile BRAF_MAPK1_Mutation_Cancer_Skin_BROAD`

*ERBB2* CEN in the context of *amplification* of *ERBB2* gene in all *tissue* types with the *BROAD* project data with *ERBB2_Amplification_Tissue_BROAD* result file name:

`python3 pyCEN.py -gene ERBB2 -amp -add_dependant -cogene ERBB2 -group tissue -project BROAD -ofile ERBB2_Amplification_Tissue_BROAD`

Other usages that require user input files:

*MITF* and *SOX9* CEN in the context of *expression* of *SOX10* gene in all *cancer* types with the *BROAD* project data with *SOX10_dependants_Expression_Cancer_BROAD* result file name:

`python3 pyCEN.py -gene SOX10 -exp -add_dependant -cogenefile ../test_data/sox10_dependants.csv -project BROAD -group tissue -ofile SOX10_dependants_Expression_Cancer_BROAD`

CEN of all kinases in the context of *BRAF* *hotspot* *mutation* in *skin* *tissue* in *BROAD* project with *BRAF_Kinases_Mutation_Skin_Tissue_BROAD*:
	
`python3 pyCEN.py -gene BRAF -add_dependant -cogenefile ../test_data/kinases.csv -mut -group tissue -groupname Skin -project BROAD -hotspot -ofile BRAF_Kinases_Mutation_Skin_Tissue_BROAD`

CEN of *WRN* in the context of MSI as a user separated cell line lists in *SANGER* project data with *WRN_MSI_SANGER* result file name:

`python3 pyCEN.py -gene WRN -add_dependant -cogene WRN -compess -separate -interested_group ../test_data/SANGER_MSI_CLs.csv -other_group ../test_data/SANGER_MSS_CLs.csv -project SANGER -ofile WRN_MSI_SANGER`
	
***In all ".csv" files, the first column should not include numbers!***
___
### Citation
If you find CEN-Tools useful please cite:  

Sharma S*, Dincer C*, Weidemüller P, Wright GJ, Petsalaki E., *CEN-tools: An integrative platform to identify the ‘contexts’ of essential genes.*
___
### Contact

Contact us [here](mailto:cen-tools@googlegroups.com) for any problems or feedback on CEN-tools.
___
### Terms and conditions
 
 **CEN-tools**  is an integrative tool that identifies the underlying context that is associated with the essentiality of a gene from large-scale genome-scale CRISPR screens.

Copyright © 2020 EMBL- European Bioinformatics Institute

This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

Neither the institution name nor the name CEN-tools can be used to endorse or promote products derived from this software without prior written permission. For written permission, please contact  [petsalaki@ebi.ac.uk](mailto:petsalaki@ebi.ac.uk).

Products derived from this software may not be called CEN-tools nor may CEN-tools appear in their names without prior written permission of the developers.

You should have received a copy of the GNU General Public License along with this program. If not, see  [gnu.org/licenses/](http://www.gnu.org/licenses/).

You can also view the licence  [here.](http://www.gnu.org/licenses/lgpl-3.0.txt).
___
### Further Disclaimer

For policies regarding the underlying data, please also refer to:

-   [DepMap: terms and conditions](https://depmap.org/portal/terms/)
-   [Project Score: terms and conditions](https://score.depmap.sanger.ac.uk/documentation)


  


