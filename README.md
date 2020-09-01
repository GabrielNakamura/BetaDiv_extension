# Data and scripts for reproduce the analysis contained in the manuscript:
["A multifaceted approach to analyzing taxonomic, functional, and phylogenetic β diversity"](https://esajournals.onlinelibrary.wiley.com/doi/abs/10.1002/ecy.3122)

## Authors

Gabriel Nakamura, Yzel Rondon Súarez, Wagner Vicentin and Leandro Duarte


## Repository structure

```shell
├── MS_BD_extend_Ecology.Rproj
├── R
│   ├── Analises_Ivinhema-16-01-2019.R # script with empirical analysis
│   ├── C_plotPhylo.R # script to plot phylogenetic tree
│   ├── C_readData_09-05-20.R # read data for empirical analysis
│   ├── barplot_LCBD.R # generate figure for plot LCBD values
│   └── unimodal_resp.R # generate figure to plot unimodal curve response
├── README.md
├── analysis_performance # folder containing all scripts to simulate the scenarios presented in manuscript
│   ├── All_scenario_lastVersion_18-1-2019.R
│   ├── All_scenarios_table.xlsx
│   ├── Metadata_folder-analysis_performance.rtf
│   ├── Performance_all.png
│   ├── Supp_AIC_analysis.R
│   ├── Supp_material_AICvalues.RData
│   ├── Table1_allScenarios_2018.docx
│   └── dissimilarity_properties # folder with scripts to test for dissimilarity properties in beta diversity metrics
│       ├── Properties_P3_6_05-03.pdf
│       ├── Toy_Wmatrices.xlsx
│       ├── dissimilarity_properties.R
│       ├── dissimilarity_properties.Rproj
│       ├── functions
│       │   ├── BetaDiv_extension_22-10-2019.R
│       │   ├── beta-div.R
│       │   ├── simul_dimensionality.R
│       │   └── simulate_metacommunity_adapt.R
│       └── graph_properties.svg
├── data # all data for empirical analysis
│   ├── Sub_portions_concatenate.xlsx
│   ├── atributos_Ivinhema.xlsx
│   ├── comu.Ivinhema.completa.csv
│   ├── dados ivinhema.xls
│   ├── metadados_atributos_Ivinhema.txt
│   └── processed
│       ├── Org_Amb.xlsx
│       ├── atributos_Ivinhema.txt
│       ├── comu_Ivinhema.txt
│       └── filogenia.ivinhema.riachos.txt
├── functions # all functions needed to run simulation and empirical analysis
│   ├── Beta.div_adapt.R
│   ├── LCBDsig_Parfunction_18-01-2019.R
│   ├── beta-div.R
│   └── simulate_metacommunity_adapt.R
└── shp_Ivinhema
    ├── drenagens_ivinhema.cpg
    ├── drenagens_ivinhema.dbf
    ├── drenagens_ivinhema.prj
    ├── drenagens_ivinhema.sbn
    ├── drenagens_ivinhema.sbx
    └── drenagens_ivinhema.shp
```



# Beta_div function
Function to compute phylogenetic or functional **Beta Diversity** (BD<sub>Fun</sub> and BD<sub>Phy</sub>). This is an extension of the method originally proposed by [Legendre and De Cáceres 2013](https://onlinelibrary.wiley.com/doi/full/10.1111/ele.12141).

### Details

Beta.div_adapt function computes functional and phylogenetic BD and its components (local and species contributions) as a variance presented in a metacommunity. The function also tests the null hyphotesis of no effect of phylogenetic and functional relationships in the distribution of species in metacommunity.

# arguments
### inputs:

Y= metacommunity composition matrix, communities in rows and species in columns;

dist_spp= matrix of pairwise distance among species;

nperm= number of permutation used to compute p values of BD<sub>F</sub> or BD<sub>P</sub> and its components;

method= method used to calculate Beta diversity and its components.  "raw" (default) or "distance.based";

### output:

List object containing the following components: BDextend.obs - BD metric extended for functional or phylogenetic dimensions (BD<sub>Fun</sub> or BD<sub>Phy</sub>);
LCBDextend.obs - local components of BDextend.obs;
SCBDextend.obs - species contribution component of BDextend.obs (applicable only when used with "raw" method);
ptaxa.BDextend and ptaxa.LCBDextend - p values calculating according to taxa shuffle null model for BDextend.obs and LCBDextend.obs
