# Beta_div function
Function to compute functional and phylogenetic **Beta Diversity** (BD) (BD<sub>F</sub> and BD<sub>P</sub>). This is an extension of the method originally proposed by [Legendre and De Cáceres 2013](https://onlinelibrary.wiley.com/doi/full/10.1111/ele.12141).

###Important
Before run the function the user must read the original function used to calculate [BD](Legendre and De Cáceres 2013](https://onlinelibrary.wiley.com/doi/full/10.1111/ele.12141) 

# arguments
###inputs:

Y= metacommunity composition matrix, communities in rows and species in columns; 

tree= phylogenetic tree, must be in newick format;

trait= vector of quantitative traits for all species in Y;

nperm= number of permutation used to compute p values of BD<sub>F</sub> and BD<sub>P</sub> and its components;

method= method used to calculate Beta diversity and its components. "distance.based" (default) or "raw";

###output:

List object containing the following components: BD.obs - functional and phylogenetic beta diversity metrics (BD<sub>F</sub> and BD<sub>P</sub>);
LCBD.obs - local components of BD<sub>F</sub> and BD<sub>P</sub>/;
p.BD and p.LCBD - p values calculating according to taxa shuffle null model.

