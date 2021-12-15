# Agent-based computational modelling of glioblastoma predicts that stromal density is central to oncolytic virus efficacy

This code accompanies the paper submitted to iScience "Agent-based computational modelling of glioblastoma predicts that stromal density is central to oncolytic virus efficacy" a nanoHub simualtion of this code can be found at https://nanohub.org/tools/ohsv1gbmapp.

The code uses PhysiCell: an Open Source Physics-Based Cell Simulator for 3-D Multicellular systems

**Reference:** A Ghaffarizadeh, R Heiland, SH Friedman, SM Mumenthaler, and P Macklin, PhysiCell: an Open Source Physics-Based Cell Simulator for Multicellular Systems, PLoS Comput. Biol. 14(2): e1005991, 2018. DOI: [10.1371/journal.pcbi.1005991](https://dx.doi.org/10.1371/journal.pcbi.1005991)

Visit http://MathCancer.org/blog for the latest tutorials and help. 

**Version:**      1.6.0

**Release date:** 20 August 2019

## Overview: 
The code is an agent-based computational model of glioblastoma tumours predicting the efficacy of oHSV-1 (Herpes Simplex Oncolytic virus). The model describes the individual interactions between glioblastoma, stromal, CD4+ T and CD8+ T cells. Proliferating glioblastoma cells become infected by oHSV-1, changing phenotype to infected cells. CD4+T cells become actived by the presence of infected cells and secrete chemokine which recruits CD8+ T cells thorugh chemotaxis. CD8+ T cells then induce apoptosis in infected cells and/or tumour cells depending on the antigen specificity. Each cell is modelled as an agent and equipped with an independent state and set of rules dictating its behaviour based on load environmental conditions (i.e. chemokine and virus concentrations) and cell-cell interactions wit virus and chemokine modelled as a diffusing field. 
![image](https://user-images.githubusercontent.com/48768705/146128391-5f090985-7f1b-4cd5-aa09-495ead8d8270.png)

### Key makefile rules:

**make**               : compiles the current project. If no 
                     project has been defined, it first 
                     populates the cancer heterogeneity 2D 
                     sample project and compiles it 
   
**make GBM-OV-immune-stroma-patchy**: populates the PhysiCell environment with the GBM-OV-immune-stroma-patchy project. 
                     Use "make" to compile it. 

**make clean**         : removes all .o files and the executable, so that the next "make" recompiles the entire project 

**make data-cleanup**  : clears out all simulation data 

**make reset**         : de-populates the sample project and returns to the original PhysiCell state. Use this when switching to a new PhysiCell sample project. 

* * * 
### Running
Once the model has been compiled, to run simply run .\GBM_OV_immune_stroma_patchy.exe

### Changing dense:sparse configurations
To change the configuration of dense to sparse regions of the GBM, simplfy change the "proportion" and update the corresponding "kappa". 
