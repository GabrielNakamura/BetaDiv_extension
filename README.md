# Beta_div function
Function to compute phylogenetic or functional **Beta Diversity** (BD<sub>Fun</sub> and BD<sub>Phy</sub>). This is an extension of the method originally proposed by [Legendre and De Cáceres 2013](https://onlinelibrary.wiley.com/doi/full/10.1111/ele.12141).

### Important

Before run the function the user must read in R workspace the original function used to calculate [BD](https://onlinelibrary.wiley.com/doi/full/10.1111/ele.12141) 

### Details

Beta.div_adapt function computes functional and phylogenetic BD and its components (local and species contributions) as a variance presented in a metacommunity. The function also tests the null hyphotesis of no effect of phylogenetic and functional relationships in the distribution of species in metacommunity.

# arguments
### inputs:

Y= metacommunity composition matrix, communities in rows and species in columns; 

dist_spp= matrix of pairwise distance among species;

nperm= number of permutation used to compute p values of BD<sub>F</sub> or BD<sub>P</sub> and its components;

method= method used to calculate Beta diversity and its components. "distance.based"  or "raw" (default);

### output:

List object containing the following components: BDextend.obs - BD metric extended for functional or phylogenetic dimensions (BD<sub>Fun</sub> or BD<sub>Phy</sub>);
LCBDextend.obs - local components of BDextend.obs;
ptaxa.BDextend and ptaxa.LCBDextend - p values calculating according to taxa shuffle null model for BDextend.obs and LCBDextend.obs

