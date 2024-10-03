The [[CheckTopCorrelationsOfEnrichInput | checking of the top correlations]]  for each clinical parameter by manual curation of the [[visualize_top_correlations.r]] plots is very time intensive and difficult to find for each of the many potentially interesting correlations the potential connection to genes/pathways/mechanisms associated to pulmonary hypertension, thus the search is limited to checking a few very good performing correlations.

The contrary approach described here aims to establish associations by finding already established connections  of PH and genes/pathways/mechanisms as well as already established and upcoming treatment targets. This was done by extracting information of:

* PH reviews ([[@montaniPulmonaryArterialHypertension2013]], [[@sommerCurrentFutureTreatments2021]])
* Consensus mining ([[Consensus overview]])
* Recent publications of Soni ([[ResearchSoni]])

A collection of potential genes of interest was created by adding the not only the direct hits in literature but also an gene name mining in the ENSEMBL database. The gene name mining is simply a search of ENSEMBL database human gene hits, that carry the gene name or the mechanism in their entry and thus have a sort of connection to the PH literature hit.
Afterwards were all of the ENSEBL gene ids searched in the [[CheckTopCorrelationsOfEnrichInput | top hits of the correlations]] to check if they are represented there. If this is the case they were added to the list [[PotentialHitsOfInterest]] (associated hits are represented in the respective gene page (GenesDump directory)).

The extracted keywords in form of mechanisms or functional structures were search in the input of the [[goEnrichment.r | enrichment analysis]] as a more targeted search compared to ENSEMBL described above.
# KeywordExtraction
The goal is to get a better unerstanding of pulmanory hypertension to extract a set of keywords, which can be used for high trougput literature mining by assotiating candidate gene names with keywords.
## Overview
This approach is a semi manual light version of what could be done for high throughput literature research for candidates. It can function as a proof of concept for the usage of consenus as data mining tool that provides an output which can be further processed by another Language model to extract and filter the relevance and compare and filter the different outputs for single gene. 
But its has to be noted that consensus is only based on public access publication and thus can be only a part of a bigger data mining application which is based on different sources and approaches like a pubmed crawler.

### (From [[@montaniPulmonaryArterialHypertension2013]])
Blood tests Serological tests for HIV, hepatitis B or C serology should be performed to screen for associated diseases. The thyroid hormone measurement may reveal either hyperthyroid dysfunction or autoimmune thyroiditis, frequently encountered in PAH.
Many different treatments target prostacyclin and its receptors as well as endothelin and its receptor due to its role as vasoconstrictor. Also targeted with PDE-5 inhibitor
### (From [[@sommerCurrentFutureTreatments2021]])
The advent of pharmacological therapies targeting the prostacyclin, endothelin, and NO pathways has significantly improved outcomes.
Current research focuses  on targeting the underlying pathways of aberrant proliferation, migration, and apoptosis. Despitesuccess in preclinical models, using a plethora of novel approaches targeting cellular GPCRs,  ion channels, metabolism, epigenetics, growth factor receptors, transcription factors, and inflammation, successful transfer to human disease with positive outcomes  in clinical trials is limited.

... However, interest in the field has been revived by an improved understanding of immune regulation in PAH and the notion that perivascular inflammatory infiltrates (macrophages, B-cells, T-cells, and dendritic cells) often precede structural pulmonary vascu lar remodelling (Tamosiuniene et al., 2011)

# ==General Keywords==
```
Pulamnory hypertension, inflammation, epitel, immune , Angiogenese, vascular, cytokine, TGF-ß signaling, vasoconstrictive/vasodilatory, vascular remodelling, NO pathways, prostacyclin, endothelin, PDE-5, proliferative/mitogenic, apoptosis, Coagulation, adhesion, growth factors, macrophages, B-cells, T-cells, dendritic cells, neutrophils, lymphocytes
```

==Maybe check also connection to Inflammatory mediators==
```
FOCUS FOCUS FOCUS
### Inflammatory mediators
## IL-1
## IL-6 
## LTB4
## macrophage migration inhibitory factor
## leptin
## TNF-α
## CD20
## elafin
## neutrophil derived serine proteases elastase and proteinase-3
```