# Idea
The Idea is to check up to 250 of the top correlations for each clinical parameter based on the (Stichit outer cv filter, only coding genes) prefiltering of the DisGeNet and Go enrichment analysis input creation.
For each clinical parameter the  [[visualize_top_correlations.r |correlation visualizations]] are manually curated und the result of DisGeNet and Go enrichment are integrated.

# Workflow
**==(use obsidian data management structure for curation documentation)==**

1. Select genes of interest for each clinical parameter (maybe 5-10)
	- Check all top correlations (enrichment input set) for:
		-
	- DisGeNet
		- Check the strongest correlations that is associated to similar (based on the name) disease to pulmonary hypertension 
		- Check the most abundant genes that play a role in many of the similar sounding disease
	- GO enrichment
		- Check the correlations of genes hits with lowest p-values

2. Manually curate the correlations
	- solid separation between the disease and healthy conditions of the patient specific clinical parameters
	- dense cluster of the patient conditions in the dotplot that displays the linear correlations
	- linear correlation that are not driven mainly by outliers
	- correlations that have a solid flow
	- only a few zero ATAC capabilities, unless mostly a single condition group is exclusive zero
	- ==Curate the segments==
		- size of the open frame
		- other regulatory elements of the same gene in direct proximity?

3.  Mine literature information
	- via consensus
		- Build more informative prompts
			- include connection between gene and:
				- PH
				- pulmonary
				- hypertension
				- inflammatory disease
				- vascular system
				- ==SEE KEEP IN MIND SECTION BELOW==
	- via scholar
	- bring in again
		- GO enrichment cycles and check the other correlations again
		- DisGeNet disease results
	- Keep in mind where the data is coming from! (blood samples):
		- we sequenced blood cells and not the affected tissues
		- the samples contain sequenced immune cells 
			- interesting because their interaction with the epithelial cells to get into the affected tissues
			- blood is not restricted to the locality of PH and strong reflect other causalities like age or other disease of the patients
			- Inflammatory effects on the tissues interesting - keyword cytokine
			- connections to inflammatory effects within the vascular system
			- Interaction of immune cells with other cells
			- ANgiogenese ebenfalls iteressant und linekd zy Hypertension

4. Generate additional visualization
	- Boxplot of the condition groups along the clinical parameters (can the parameter be used to distinct patient groups in the dataset) - maybe include the whole clinical data of all 46 patient to have more informative parameter distributions
	- Boxplot of the ATAC accessibilities across the patient groups
	- improved dot-plots
		- exclude outliers
		- use square triangle and dot for the visualization together with the color coding for the condition groups
		- colorblind color coding
		- bigger dots, higher resolution, bigger and more informative text
	- IGV visualization of the gene loci together with the regulatory positions
		- check the position other the other reported regulatory segments

5. Gather further infos
	- GWAS, TF footprint