### Cross-Tissue-WGCNA-Brain-Blood-Transcriptomic-Signatures-for-Asymptomatic-AD-Diagnosis.-
# Status: Manuscript under review 

---


## Abstract
Background: Alzheimer’s disease (AD) is a neurodegenerative disorder characterised by a prolonged prodromal stage, often resulting in late diagnosis and limited therapeutic efficacy. Early diagnosis at the asymptomatic stage (AsymAD), before irreversible synaptic loss, is therefore critical. Although peripheral blood offers potential for early AD detection, biomarker discovery is limited by systemic signal dilution, dependence on symptomatic cohorts, and conventional differential expression analysis, which fails to detect subtle prodromal alterations.
Methods: To address these methodological gaps, a multi-tiered, cross-tissue Weighted Gene Co-Expression Network Analysis (WGCNA) pipeline was developed. Unbiased co-expression networks anchored in the entorhinal cortex of AsymAD patients were constructed, and cross-tissue preservation in peripheral blood was assessed. Prioritised hub genes from the preserved module were used to train a Random Forest classifier, which was validated on an independent peripheral blood microarray dataset.
Results: Of ten modules identified, one module (green) was preserved in peripheral blood (Z-summary = 15.25). This module, associated with dysregulation of mitochondrial energy and proteostasis pathways, was refined into a five-gene screening panel (ATP6AP2, CKS1B, STAMBPL1, SUB1, and BET1). A Random Forest classifier trained with this panel distinguished AsymAD patients from healthy controls, achieving an AUC of 0.722 and an accuracy of 69.1%.
Conclusion: The findings present anatomically anchored and computationally predictive cross-tissue transcriptomic signatures that warrant further investigation as potential blood-based early biomarkers for AD diagnosis.

## Methodology
To capture subtle biological changes without the bias of traditional DEG thresholds, the pipeline used Median Absolute Deviation to select the most variable genes and carried pou Weighted Gene Co-Expression Network Analysis (WGCNA) witth bulk RNA-seq data from entorinal and frontal corex to identify modules with coexpression networks for aysymptomatic and symptomatic stages of AD. 
Cross-tissue preservation analysis was conducted using the Z-summary statistic to mathematically confirm whether the topological architecture of the brain-derived modules persisted in a peripheral whole-blood cohort.
Prioritized hub genes from the preserved network were used to train a non-linear Random Forest machine-learning classifier.
The predictive efficacy of this model was then validated on an independent peripheral blood microarray dataset.

## Key Findings
Network Preservation: Out of ten structural modules identified in the early-stage entorhinal cortex, the "green" module (containing 973 genes) demonstrated highly significant structural preservation across tissues, achieving a systemic Z-summary score of 15.25 in peripheral blood.
Biological Mechanism: Functional enrichment analysis revealed that the genes within this preserved green module are associated with the systemic dysregulation of mitochondrial energy, protein-folding machinery, and p53-mediated apoptotic signaling pathways.
Biomarker Panel Identification: Through a dual-metric evaluation of Module Membership and Gene Significance, the 973-gene network was refined into a highly targeted, five-gene screening panel: ATP6AP2, CKS1B, STAMBPL1, SUB1, and BET1.
Diagnostic Performance: The Random Forest classifier trained on these five systemic hub genes successfully distinguished AsymAD patients from healthy controls, achieving an overall accuracy of 69.1%.
Predictive Validation: Receiver Operating Characteristic (ROC) analysis of the model yielded an Area Under the Curve (AUC) of 0.722.
Feature Importance: SUB1 and CKS1B emerged as the primary predictive drivers in the blood-based classification, with both dynamic hub genes exhibiting distinct downregulation in the peripheral blood of the AsymAD cohort.


## Tech Stack
Programming language: R
Network construction: WGCNA  
Pathway & Enrichment Analysis: clusterProfiler, GSEA  
Statistical Modelling: prcomp (PCA), PROC package (ROC Curve Analysis)   
Data Visualization: ggplot2  
