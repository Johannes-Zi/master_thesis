The goal is to get a better unerstanding of pulmanory hypertension to extract a set of keywords, which can be used for high trougput literature mining by assotiating candidate gene names with keywords.

# Concept that can be described in the Master thesis
## Overview
This approach is a semi manual light version of what could be done for high throughput literature research for candidates. It can function as a proof of concept for the usage of consenus as data mining tool that provides an output which can be further processed by another Language model to extract and filter the relevance and compare and filter the different outputs for single gene. 
But its has to be noted that consensus is only based on public access publication and thus can be only a part of a bigger data mining application which is based on different sources and approaches like a pubmed crawler.

### From the 2013 review
Blood tests Serological tests for HIV, hepatitis B or C serology should be performed to screen for associated diseases. The thyroid hormone measurement may reveal either hyperthyroid dysfunction or autoimmune thyroiditis, fre quently encountered in PAH.

Many different treatments target prostacyclin and its receptors as well as endothelin and its receptor due to its role as vasoconstrictor. Also targeted with PDE-5 inhibitor
## ==General Keywords==
```
Pulamnory hypertension, inflammation, epitel, immune , Angiogenese, vascular, cytokine, oxidative stress, TGF-ß signaling, vasodilator, vasoconstrictor, thyroid, prostacyclin, endothelin, PDE-5
```

# Research from Soni
## Overview
==In this step i screened systematically the top 1720 correlations (DisGeNet and Go enrich input) if they are part of the top gene dataset.==

We hypothesize that ADAR1-mediated RNA editing via MDA5-PKR IFN-β activation plays a critical role in abnormal SMC proliferation and autoimmune inflammation during vascular remodeling.
*(Kim, Y., et al. "PKR Activation Caused by ADAR1 Deficeincy Contributes to Vascular Remodeling in PAH." D95. SURVEYING THE ZOO:-OMICS FOR MECHANISTIC DISCOVERY IN PULMONARY VASCULAR DISEASE. American Thoracic Society, 2024. A7244-A7244.)*

We hypothesized that NLRP3 inflammasome activation mediates RVCM contractile dysfunction and that this is attenuated by E2 via estrogen receptor (ER) α.
*(Sobrano Fais, R., et al. "NLRP3-induced Contractile Dysfunction in Right Ventricle Cardiomyocytes Is Sexually Dimorphic and Modified by 17β-estradiol Via Estrogen Receptor α." B95. SURFING IN THE GREEN ROOM: DEFINING MOLECULAR MECHANISMS OF PULMONARY VASCULAR DISEASE. American Thoracic Society, 2024. A4695-A4695.)*

Single-cell RNA sequencing on sorted lineage-labeled GLI1+ cells revealed an _Acta2__high_ fraction of cells with pathways in cancer and MAPK (mitogen-activated protein kinase) signaling as potential players in reprogramming these cells during vascular remodeling. Our data highlight GLI1+ cells as an alternative cellular source of VSMCs (Vascular smooth muscle cell) in pulmonary hypertension and suggest that these cells and the associated signaling pathways represent an important therapeutic target for further studies.
*(Chu, Xuran, et al. "GLI1+ Cells Contribute to Vascular Remodeling in Pulmonary Hypertension." Circulation Research 134.11 (2024): e133-e149.)*

Genetic diagnostics and molecular approaches in pulmonary arterial hypertension
![[Pasted image 20240925142515.png]]
(BMPR2, ATP13A3, AQP1, ABCC8, KCNK3, SMAD9, SOX17, CAV1, TBX4, EIF2AK4, KDR, ENG, ACVRL1, GDF2)
*(Eichstaedt, Christina A., et al. "Genetische Diagnostik und molekulare Ansätze bei pulmonalarterieller Hypertonie." Pneumologie 77.11 (2023): 862-870.)*

The circulating levels of several extracellular matrix proteins deregulated in decompensated RV subgroups were assessed in two independent cohorts of individuals with pulmonary arterial hypertension, revealing that NID1, C1QTNF1 and CRTAC1 predicted the development of a maladaptive RV state, as defined by magnetic resonance imaging parameters, and were associated with worse clinical outcomes.
*(Khassafi, Fatemeh, et al. "Transcriptional profiling unveils molecular subgroups of adaptive and maladaptive right ventricular remodeling in pulmonary hypertension." Nature cardiovascular research 2.10 (2023): 917-936.)*

## ==Keywords==
==Ive checked for the following genes if they were found in the set of 1720 target genes (DisGeNet input)==
```
# REspective ENSEMBL IDs OR of IDs of related genes

ADAR1, ENSG00000160710
MDA5-PKR IFN-β, ENSG00000115267
SMC, ENSG00000231822/ENSG00000268364/ENSG00000270332
NLRP3, ENSG00000162711
estrogen,
GLI1+, ENSG00000111087
Acta2, ENSG00000107796
MAPK, ENSG00000250062/ ENSG00000259438/ ENSG00000261399/ ENSG00000157833/ ENSG00000282110/ ENSG00000225663/ ENSG00000141441/ ENSG00000168175/ 
VSMCs,
NID1, ENSG00000116962
C1QTNF1, ENSG00000173918
CRTAC1, ENSG00000095713
PKA, ENSG00000095713
NOCT, ENSG00000151014
IL6, ENSG00000136244
IFNa2 and IFNAR2, ENSG00000188379 and ENSG00000159110
MTC1, ENSG00000232119,

# From the treatment guidelines
BMPR2, ENSG00000204217
ATP13A3, ENSG00000133657
AQP1, ENSG00000240583
ABCC8, ENSG00000006071
KCNK3, ENSG00000171303
SMAD9, ENSG00000120693
SOX17, ENSG00000164736
CAV1, ENSG00000105974
TBX4, ENSG00000121075
EIF2AK4, ENSG00000128829
KDR, ENSG00000128052
ENG, ENSG00000106991
ACVRL1, ENSG00000139567
GDF2, ENSG00000263761

```
# CONSENSUS - Connection between blood cells and pulmanory hypertension
## Overview
### Introduction

Pulmonary hypertension (PH) is a complex and multifactorial disease characterized by elevated blood pressure in the pulmonary arteries, leading to right heart failure and other severe complications. Recent research has explored the potential role of blood cells, particularly immune cells, in the pathogenesis and progression of PH.

### Key Insights

- **Immune Cells and Hypertension**:
    
    - Immune cells, including lymphocytes, monocytes, and neutrophils, have been associated with increased blood pressure, suggesting a potential causal relationship between these cells and hypertension[1](https://consensus.app/results/?q=is%20there%20a%20conncetion%20between%20blood%20cells%20and%20pulmanory%20hypertension%3F&synthesize=on&copilot=on&lang=en#result-1) [3](https://consensus.app/results/?q=is%20there%20a%20conncetion%20between%20blood%20cells%20and%20pulmanory%20hypertension%3F&synthesize=on&copilot=on&lang=en#result-3) [4](https://consensus.app/results/?q=is%20there%20a%20conncetion%20between%20blood%20cells%20and%20pulmanory%20hypertension%3F&synthesize=on&copilot=on&lang=en#result-4).
    - Both innate and adaptive immune cells, such as macrophages, T cells, and dendritic cells, contribute to hypertension and related end-organ damage[3](https://consensus.app/results/?q=is%20there%20a%20conncetion%20between%20blood%20cells%20and%20pulmanory%20hypertension%3F&synthesize=on&copilot=on&lang=en#result-3) [4](https://consensus.app/results/?q=is%20there%20a%20conncetion%20between%20blood%20cells%20and%20pulmanory%20hypertension%3F&synthesize=on&copilot=on&lang=en#result-4).
      *(Adventitial dendritic cells play a critical role in the initiation of vasculitis by switching from nonstimulatory to T-cell activating, regulating T-cell recruitment and retention in arterial walls.*
      *Vascular Dendritic Cells in Giant Cell Arteritis; C. Weyand, W. Ma-Krupa, O. Pryshchep; 1. Dez. 2005)*
- **Inflammatory Cells in Pulmonary Hypertension**:
    
    - Inflammatory cell infiltrates, including various immune effector cells, are a prominent feature in the vascular remodeling seen in pulmonary arterial hypertension (PAH), indicating their role in disease progression[2](https://consensus.app/results/?q=is%20there%20a%20conncetion%20between%20blood%20cells%20and%20pulmanory%20hypertension%3F&synthesize=on&copilot=on&lang=en#result-2).
    - Genome-wide blood expression profiling in PAH patients has identified enriched mechanisms involving interferon, mTOR/p70S6K, stress kinase, and toll-like receptor signaling, which are consistent with inflammatory and immune responses[2](https://consensus.app/results/?q=is%20there%20a%20conncetion%20between%20blood%20cells%20and%20pulmanory%20hypertension%3F&synthesize=on&copilot=on&lang=en#result-2).
- **Specific Immune Cell Contributions**:
    
    - Neutrophils and their extracellular traps (NETs) are implicated in hypertension, suggesting that these cells may also play a role in pulmonary hypertension through similar inflammatory mechanisms[3](https://consensus.app/results/?q=is%20there%20a%20conncetion%20between%20blood%20cells%20and%20pulmanory%20hypertension%3F&synthesize=on&copilot=on&lang=en#result-3).
    - γδ T cells have been shown to mediate angiotensin II-induced hypertension and vascular injury, indicating that specific T cell subsets may be involved in the immune response contributing to hypertension[5](https://consensus.app/results/?q=is%20there%20a%20conncetion%20between%20blood%20cells%20and%20pulmanory%20hypertension%3F&synthesize=on&copilot=on&lang=en#result-5).

### Conclusion

The research indicates a significant connection between blood cells, particularly immune cells, and pulmonary hypertension. Immune cells such as lymphocytes, monocytes, neutrophils, and specific T cell subsets contribute to the inflammatory processes that underlie both systemic and pulmonary hypertension. These findings suggest that targeting immune cell-mediated inflammation could be a potential therapeutic strategy for managing pulmonary hypertension.

## ==Keywords==
```
Immune cells, lymphocytes, monocytes, neutrophils, blood pressure, macrophages, T cells, dendritic cells, interferon, mTOR/p70S6K, stress kinase, toll-like receptor signaling
```
# CONSENSUS - Current treatment targets
## Overview
```
	Is there a conncetion between geneXXX and the endothelin system or endothelin receptors?
```
**In Usage Drug Classes for PAH**:
These drug classes have significantly improved the management of PAH, enhancing progression-free survival and quality of life
- Endothelin Receptor Antagonists,
- Phosphodiesterase-5 Inhibitors,
- Soluble Guanylate Cyclase Stimulators,
- Prostacyclin Analogues,
- Prostacyclin Receptor Agonists
## ==Keywords==

```
Endothelin, Phosphodiesterase, Guanylate Cyclase, Prostacyclin
```
# CONSENSUS - Emerging Therapies:
## Overview
- **Novel Targets**: Research is focusing on targeting pathways involved in vascular remodeling, including GPCRs, ion channels, metabolism, epigenetics, growth factor receptors, and inflammation[3](https://consensus.app/results/?q=pulmonary%20hypertension%20treatments&synthesize=on&copilot=on&lang=en#result-3) [8](https://consensus.app/results/?q=pulmonary%20hypertension%20treatments&synthesize=on&copilot=on&lang=en#result-8).
- **Phytochemicals and Medicinal Plants**: These have shown potential due to their antiproliferative, antioxidant, antivascular remodeling, anti-inflammatory, vasodilatory, and apoptosis-inducing actions[5](https://consensus.app/results/?q=pulmonary%20hypertension%20treatments&synthesize=on&copilot=on&lang=en#result-5).
*(Sommer, Natascha, et al. "Current and future treatments of pulmonary arterial hypertension." British journal of pharmacology 178.1 (2021): 6-30.;*
*Baliga, Reshma S., Raymond J. MacAllister, and Adrian J. Hobbs. "New perspectives for the treatment of pulmonary hypertension." British journal of pharmacology 163.1 (2011): 125-140.)*
## ==Keywords==

```
vascular remodeling, GPCRs, ion channels, epigenetics, growth factor receptors, inflammation, apoptosis
```

# CONSENSUS - Blood cells and pulmanory arterial hypertension V2
## Overview
### Introduction

Pulmonary arterial hypertension (PAH) is a severe condition characterized by increased blood pressure in the pulmonary arteries, leading to right heart failure and premature death. Recent research has explored the role of various blood cells in the pathogenesis and progression of PAH, aiming to identify potential biomarkers and therapeutic targets.

### Key Insights

- **Inflammatory and Immune Cells:**
    
    - Inflammatory cell infiltrates, including macrophages, T cells, and dendritic cells, are prominent in PAH, suggesting that immune effector cells contribute to disease progression[1](https://consensus.app/results/?q=is%20there%20a%20conncetion%20between%20blood%20cells%20and%20pulmanory%20arterial%20hypertension%3F&synthesize=on&copilot=on&lang=en#result-1) [5](https://consensus.app/results/?q=is%20there%20a%20conncetion%20between%20blood%20cells%20and%20pulmanory%20arterial%20hypertension%3F&synthesize=on&copilot=on&lang=en#result-5) [9](https://consensus.app/results/?q=is%20there%20a%20conncetion%20between%20blood%20cells%20and%20pulmanory%20arterial%20hypertension%3F&synthesize=on&copilot=on&lang=en#result-9).
    - Abnormal T lymphocyte subsets, particularly increased CD8+ cytotoxic effector-memory cells and regulatory T cells, are observed in PAH patients, indicating a dysfunctional immune system[5](https://consensus.app/results/?q=is%20there%20a%20conncetion%20between%20blood%20cells%20and%20pulmanory%20arterial%20hypertension%3F&synthesize=on&copilot=on&lang=en#result-5) [9](https://consensus.app/results/?q=is%20there%20a%20conncetion%20between%20blood%20cells%20and%20pulmanory%20arterial%20hypertension%3F&synthesize=on&copilot=on&lang=en#result-9).
- **Red Blood Cells (RBCs):**
    
    - PAH patients exhibit decreased endothelial nitric oxide synthase (eNOS) activity and nitric oxide (NO) generation in RBCs. This reduction is linked to increased Rho-Kinase (ROCK) activity, which can be mitigated by ROCK inhibitors[3](https://consensus.app/results/?q=is%20there%20a%20conncetion%20between%20blood%20cells%20and%20pulmanory%20arterial%20hypertension%3F&synthesize=on&copilot=on&lang=en#result-3).
    *(Im Menschen wird eNOS vornehmlich in [Endothelzellen](https://de.wikipedia.org/wiki/Endothel "Endothel") gebildet, welche die innerste Zellschicht in Blut- und Lymphgefäßen darstellen. Dort spielen eNOS und Stickstoffmonoxid bei der Regulation des [Blutdruckes](https://de.wikipedia.org/wiki/Blutdruck "Blutdruck") und für die Funktion von Blutgefäßen eine zentrale Rolle. Verringerte Aktivität oder eine Fehlfunktion der eNOS begünstigen die Entstehung von Gefäßerkrankungen wie [Atherosklerose](https://de.wikipedia.org/wiki/Arteriosklerose "Arteriosklerose"). Besonders wichtig ist dabei das Phänomen der eNOS-Entkopplung. Entkoppelte eNOS produziert [Superoxid](https://de.wikipedia.org/wiki/Hyperoxide "Hyperoxide") an Stelle von Stickstoffmonoxid, fördert dadurch oxidativen Stress im Endothel und schadet so Blutgefäßen mehr als ihnen zu nützen. Aufgrund dieser Doppelfunktion wird eNOS auch als ein [janusköpfiges](https://de.wikipedia.org/wiki/Ianus "Ianus") Enzym bezeichnet.)*
- **Endothelial Cells:**
    
    - Endothelial cell dysfunction plays a crucial role in PAH by affecting vasoconstriction, inflammation, and vascular remodeling. This dysfunction is associated with various signaling pathways, including BMPR2, mTOR, and toll-like receptor signaling[1](https://consensus.app/results/?q=is%20there%20a%20conncetion%20between%20blood%20cells%20and%20pulmanory%20arterial%20hypertension%3F&synthesize=on&copilot=on&lang=en#result-1) [6](https://consensus.app/results/?q=is%20there%20a%20conncetion%20between%20blood%20cells%20and%20pulmanory%20arterial%20hypertension%3F&synthesize=on&copilot=on&lang=en#result-6) [7](https://consensus.app/results/?q=is%20there%20a%20conncetion%20between%20blood%20cells%20and%20pulmanory%20arterial%20hypertension%3F&synthesize=on&copilot=on&lang=en#result-7).
- **Genetic and Molecular Pathways:**
    
    - Genetic predispositions, such as mutations in the BMPR2 signaling pathway, combined with chronic inflammation and metabolic impairments, contribute to the abnormal vascular response in PAH[2](https://consensus.app/results/?q=is%20there%20a%20conncetion%20between%20blood%20cells%20and%20pulmanory%20arterial%20hypertension%3F&synthesize=on&copilot=on&lang=en#result-2) [4](https://consensus.app/results/?q=is%20there%20a%20conncetion%20between%20blood%20cells%20and%20pulmanory%20arterial%20hypertension%3F&synthesize=on&copilot=on&lang=en#result-4) [7](https://consensus.app/results/?q=is%20there%20a%20conncetion%20between%20blood%20cells%20and%20pulmanory%20arterial%20hypertension%3F&synthesize=on&copilot=on&lang=en#result-7).
    - MicroRNAs, such as miR-7110, have been identified in PAH and are involved in regulating gene expression related to vascular pathophysiology[8](https://consensus.app/results/?q=is%20there%20a%20conncetion%20between%20blood%20cells%20and%20pulmanory%20arterial%20hypertension%3F&synthesize=on&copilot=on&lang=en#result-8).

### Conclusion

There is a significant connection between blood cells and pulmonary arterial hypertension. Inflammatory and immune cells, particularly T lymphocytes, play a critical role in the disease's pathogenesis. Additionally, red blood cells exhibit altered nitric oxide production, contributing to vascular dysfunction. Endothelial cell dysfunction and genetic factors further exacerbate the condition. Understanding these cellular and molecular mechanisms can help in developing targeted therapies for PAH.

## ==Keywords==

```
lymphocyte, CD8+, cytotoxic effector-memory cells, nitric oxide synthase (eNOS), nitric oxide (NO), Rho-Kinase (ROCK), vasoconstriction, BMPR2, mTOR, and toll-like receptor signaling, miR-7110
```

# Review 2013
cited: [[@montaniPulmonaryArterialHypertension2013]]
## ==Keywords== Genetic dispositions

```
BMPR2, ACVRL1, Smad8, Smad1, Smad5, caveolin-1
```
## Potentially affected pathway
```
TGF-β signaling pathway
```
This supports the hypothesis that mutations in genes involved in the TGF-β signaling pathway may be a trigger for pulmon ary vascular remodeling. Moreover, this signaling path way controls growth, differentiation and apoptosis of various cell types like pulmonary vascular endothelial cells (ECs) and smooth muscle cells (SMCs). Thereby, mutations in genes involved in the TGF-β signaling pathway may be responsible for abnormal proliferation of pulmonary vascular SMCs and may promote ECs apoptosis, which might lead to the selection of apoptosis resistant cells and formation of plexiform lesions, the hallmark of idiopathic PAH [131,132].
## ==Keywords==
```
BMPR2, ACVRL1, Smad8, Smad1, Smad5, caveolin-1, TGF-β signaling pathway
```
## Cellular Factors
Cellular factors Proliferation of smooth muscular cells in the small per ipheral pulmonary arteries is a common characteristic in all forms of PAH. In hypoxic models, fibroblasts of the adventitia migrate to the media and intima, where prolif eration and production of matrix proteins are observed [138]. Neovascularization, mainly of the adventitia, oc curs concomitantly to the thickening of the vascular walls [139]. In response

The stimuli for endothelial prolif eration is still unknown but several factors have been in criminated such as hypoxia, inflammation, shear stress, drugs, viral infections and genetic susceptibility.
The CXCL12/CXCR4 axis may play an important role in the pulmonary recruitment of these circulating progenitors and can be therapeutically targeted

Inflammatory mechanisms seems to play an important role in certain forms of PAH such as PAH associated with auto immune diseases or HIV infection

Thirty to 40% of patients with PAH have circulating auto-antibodies and elevated plasma concentrations of pro-inflammatory cytokines such as interleukin 1 (IL-1) and interleukin-6 (IL-6), and chemokines such as fractalkine and MCP-1 [145,146]. Inflammatory cells, such as lymphocytes B and T, mac rophages, mastocytes and dendritic cells, can also be found in plexiform lesions of severe PAH [147,148]. Chemokines, like RANTES and fractalkine are also overly expressed in the pulmonary vascular endothelium of PAH patients [145].

In response to certain stimuli, platelets can produce prothrombotic, vasoactive or mitogenic factors, such as thromboxane A2 (TXA2), platelet-derived growth factor (PDGF), serotonin (5 hydroxytryptamine, 5-HT), transforming growth factor beta (TGF-β) and vascular endothelial growth factor (VEGF) that participate in vasoconstriction and vascular remodeling

One can propose that deregulated and unresolved pulmonary inflammation on the background of a genetic predispos ition, could result in persisting vascular remodelling leading to PAH. An initial acute inflammation that is normally expected to resolve with return to homeostasis, could conduct the production of auto-antibodies against vascular wall components, and would shift to chronic persisting and chronic inflammation, endothelial barrier breakdown, infiltration by immune cells, local and chronic autoimmunity, and vascular remodeling culmin ating in PAH

Vasoconstriction has been associated with an abnormal function or Page 10 of 28 Montani et al. Orphanet Journal of Rare Diseases 2013, 8:97 http://www.ojrd.com/content/8/1/97 expression of potassium channels and with endothelial dysfunction [107]. Endothelial dysfunction results in a decreased production of vasodilators such as nitric oxide (NO) and prostacyclin and an increased production of vasoconstrictors such as endothelin-1 [153]. Prostacyclin (prostaglandin I2) is a potent pulmon ary vasodilator that acts via the cyclic adenosine monophosphate (cAMP) pathway. It inhibits the pro liferation of smooth muscle cells and decreases plate let aggregation. Production of prostacyclin is reduced in endothelial cells of patients with PAH [154]

NO is also a pulmonary vasodilator which acts via the cyclic guanosine monophosphate (cGMP) pathway. To increase pulmonary vasodilatation de pendant on NO, a recent therapeutic strategy has targeted type 5 phosphodiesterase which degrades cGMP. Sildenafil or tadalafil, type 5 phosphodiesterase inhibitors, have proven their efficacy in patients with PAH [155]. Vasoactive intestinal peptide (VIP) is a neurotransmitter that has systemic and pulmonary vasodilator properties. It also inhibits smooth cell proliferation and decreases platelet aggregation and acts via the activation of the cAMP and cGMP sys tems [156].

By ligating the ETRA, ET-1 intracellular calcium concentrations increase and activates the protein kinase C pathway [157]. ET-1 is a potent pul monary vasoconstrictor and stimulates mitosis of arterial smooth muscle cells, thus contributing to pulmonary vascu lar remodeling.

Rho protein A and Rho kinases in the vasoconstriction and vascular remod eling of PAH, 5-HTT/RhoA/Rho

Hypoxia inducible factor-1 (HIF-1) is a transcription factor that principally regulates cellular adaptation to hypoxia but also regulates several genes implicated in angiogenesis, erythropoiesis, cellular metabolism and survival [168].

==In conclusion, the pathophysiology of PAH is hetero geneous and multifactorial. The genetic mutations found in familial PAH and in a proportion of sporadic PAH are neither necessary nor sufficient for the development of PAH. Therefore, the current hypothesis is that of a gen etic predisposition for PAH followed by a superimposed environmental factor (infection, inflammation, auto immunity).==

cited: [[@montaniPulmonaryArterialHypertension2013]]

## ==Keywords==
```
CXCL12, CXCR4, auto immune, auto-antibodies, plasma concentrations, interleukin 1 (IL-1), interleukin-6 (IL-6), fractalkine, MCP-1, RANTES, fractalkine, thrombosis, shear stress, A2 (TXA2), PDGF, serotonin, TGF-β, VEGF, potassium, Prostacyclin, vasodilator, vasoconstrictor, cGMP, phosphodiesterase, Vasoactive intestinal peptide, Rho
```

## Known therapeutic targets
Several specific therapeutic agents were developed for the medical man agement of PAH including prostanoids, endothelin receptor antagonists and phosphodiesterase type 5 inhibitors. Furthermore, emerging treatments such as tyrosine kinase inhibitors, soluble guanylate cyclase activators (riociguat) and prostacyclin receptor agonists (selexipag) are cur rently being evaluated in PAH
## ==Keywords==
```
endothelin, phosphodiesterase type 5, tyrosine kinase, soluble guanylate cyclase, prostacyclin (+ receptor)
```

Calcium channel blockers (CCB) are indicated in pa tients with a positive vasodilatation challenge test after inhaled NO. CCB are vasodilators and were initially

Endothelium-derived prostaglandin I2 (PGI2), or prostacyc lin, is an arachidonic acid produced by endothelial cells. Prostacyclin is a powerful systemic and pulmonary vaso dilator and an inhibitor of platelet aggregation through the increase in intracellular cyclic adenosine monophosphate (cAMP)

There are two existing isoforms of ET-1 recep tors: endothelin A (ETRA) and endothelin B (ETRB). Activation of ETRA and ETRB on pulmonary artery smooth muscle cells induce proliferation and vasoconstriction, whereas activation of ETRB on pulmonary endothelial cells leads to release of NO and prostacyclin and partici pate to the clearance of circulating ET-1

Phosphodiesterase type-5 inhibitors NO is a potent pulmonary arteries SMC relaxant that disposed vasodilator activity through up-regulation of its associated down-stream signalling molecule, cyclic GMP (cGMP), metabolism of which is dependent on the activation of a number of PDEs [208]. Phospho diesterase type 3, 4 and 5 are the three main types of this enzyme found in pulmonary artery contractive cells. PDE-5 is the most abundantly expressed isoform in pulmonary circulation which was confirmed by several experimental investigations showing a benefi cial effect of PDE-5 inhibitors on vascular remodel ling and vasodilatation [234,235].

Sildenafil is an oral PDE-5 inhibitor that is available in Europe since 2005

Riociguat is a first in-class drug that augments cGMP biosynthesis through direct stimulation of the enzyme soluble guanylate cyclase (sGC) promoting vasodilatation by direct stimulation of sGC in an NO-independent fashion, and by sensitization of sGC to low endogenous NO levels [250].

Selexipag is an orally active prodrug metabolized to the highly selective prostacyclin IP receptor agonist ACT-333679

One of the most promising targets in PAH is platelet derived growth factor (PDGF). PDGF has been implicated in endothelial cell dysfunction and prolifera tion and migration of smooth muscle cells.

## ==Keywords==
```
Calcium channel, cAMP, 
ET-1 as EDN1 Endothelin 1, ENSG00000078401
ETRA, Endothelin receptor type A, ENSG00000151617
ETRB, Endothelin receptor type B, ENSG00000136160



# All checked and added to the interesting gene results
# Assotiated genes to soluble guanylate cyclase
soluble guanylate cyclase, GUCY1A1 ENSG00000164116, GUCY1A2 ENSG00000152402, GUCY1B1 ENSG00000061918, GUCD1 ENSG00000138867, GUCY1B2 ENSG00000293412, GUCA1B ENSG00000112599,  GUCY1B2 ENSG00000123201, GUCY2C ENSG00000070019, GUCA1C ENSG00000138472, GUCA1A ENSG00000048545, GUCY2F ENSG00000101890, GUCY2D ENSG00000132518, NPPB ENSG00000120937, NPR2 ENSG00000159899, NPPC ENSG00000163273, NPPA ENSG00000175206, NHERF4 ENSG00000172367

# Genes assotiated to prostacyclin/ prostaglandin
# ALL ALREADY CHECKED IF IN TOP CORRELATIONS
# PTGIR ENSG00000160013, # GNG10 ENSG00000242616, # PTGIS ENSG00000124212, # GNG11 ENSG00000127920, # GNG4 ENSG00000168243, # GNG12 ENSG00000172380, # GNG5 ENSG00000174021, # GNGT2 ENSG00000167083, # GNG8 ENSG00000167414, # GNG7 ENSG00000176533, # PTGIR ENSG00000160013, # GNG10 ENSG00000242616, # GNGT1 ENSG00000127928, # IGFBP7 ENSG00000163453, # GNB2 ENSG00000172354, # GNB5 ENSG00000069966, # GNB3 ENSG00000111664, # GNB4 ENSG00000114450, # GNG2 ENSG00000186469, # GNAS ENSG00000087460, # GNB1 ENSG00000078369

# Genes assotiated to Tyrosine kinase
-> noch via ensembl suche hinzufügen
```