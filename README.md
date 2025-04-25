## Github Repository for
# Activation dynamics and assembly of root zone soil bacterial communities in response to stress-associated phytohormones 
## by Sreejata Bandopadhyay, Oishi Bagchi, and Ashley Shade
<i>This work is deposited as a preprint.</i>


### Data
Both the raw read data and metagenome assemblies for this study are available through NCBI under bioproject [PRJNA932434](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA932434/)

### To cite this work or code

Bandopadhyay S, O Bagchi, and A Shade. 2025. Activation dynamics and assembly of root zone soil bacterial communities in response to stress-associated phytohormones. bioRXiv [doi 10.1101/2025.04.23.650272](https://doi.org/10.1101/2025.04.23.650272)   

### Abstract

Plants can “cry for help” to recruit supportive microbiome members during stressful conditions. We evaluated the activation dynamics of root zone soil bacteria in response to phytohormones produced when plants are stressed, hypothesizing that the activated taxa support plant resilience. We conducted a 2-week laboratory experiment using mesocosms of root zone soil collected from two different crops: the annual legume common bean (Phaseolus vulgaris L.) and the perennial grass switchgrass (Panicum virgatum). We inactivated the microbiome by drying and then treated the soils with either abscisic acid, salicylic acid, a carrier control (methanol), or water, and then quantified the reactivation dynamics of bacterial populations over time, at one, 7, and 14 days after phytohormone addition, using amplicon sequencing of 16S rRNA and rDNA. There were several Actinobacterial taxa that switched from an average population-inactive to a population-active state after exposure to abscisic acid and salicylic acid, with Microbispora lineages switching especially noted. Some taxa were activated only in one crop’s soil, and some were activated in both crops’ soils in response to the same phytohormone. This work suggests that different bacteria have different specificities to phytohormones as plant stress signals and provides insights into understanding the mechanisms by which stressed plants may “cry for help” to recruit bacteria from the root zone to the rhizosphere.

### Contents

Code is split up into two directories: [Sequence_processing](https://github.com/ShadeLab/PAPER_Dormancy_resuscitation_phytohormone_mesocosm_Bandopadhyay2025/tree/main/Sequence_processing/dna_cdna_analysis_qiime2_final.txt) and [Analysis](https://github.com/ShadeLab/PAPER_Dormancy_resuscitation_phytohormone_mesocosm_Bandopadhyay2025/tree/main/Analysis/R_analysis_phytohormoneResusc_cleaned_Final.R).

#### Sequence processing
Code used for sequence processing in QIIME2 including denoising, clustering, and taxonomy assignment can be found under [Sequence_processing](https://github.com/ShadeLab/PAPER_Dormancy_resuscitation_phytohormone_mesocosm_Bandopadhyay2025/tree/main/Sequence_processing/dna_cdna_analysis_qiime2_final.txt). Script was run using the MSU HPCC. Output files such as OTU table and taxonomy files were imported in R for subsequent analysis.

#### Analysis
Formal analysis can be found under [Analysis](https://github.com/ShadeLab/PAPER_Dormancy_resuscitation_phytohormone_mesocosm_Bandopadhyay2025/tree/main/Analysis/R_analysis_phytohormoneResusc_cleaned_Final.R). All analysis was run in R version 4.1.2. The analysis directory contains the output figure files in a separate directory along with the R code. A single R file was used for all analyses. 

### Funding
Support for this research was provided by the [United States National Science Foundation under Grant No. MCB #1817377](https://www.nsf.gov/awardsearch/showAward?AWD_ID=1817377&HistoricalAwards=false) to A.S. Additional support was provided by the [Great Lakes Bioenergy Research Center](https://www.glbrc.org/), U.S. Department of Energy, Office of Science, Office of Biological and Environmental Research under Award Number DE-SC0018409 and by the National Science Foundation Long-Term Ecological Research Program (DEB #1832042). Additional support was provided by the Michigan State University Plant Resilience Institute by, the USDA National Institute of Food and Agriculture and Michigan State University AgBioResearch. AS acknowledges project support from the European Union (ERC, [MicroRescue, 101087042](https://cordis.europa.eu/project/id/101087042). Views and opinions expressed are however those of the author(s) only and do not necessarily reflect those of the European Union or the European Research Council. Neither the European Union nor the granting authority can be held responsible for them.

### More info
[ShadeLab](http://ashley17061.wixsite.com/shadelab/home)
[Ecologie Microbienne Lyon](https://www.ecologiemicrobiennelyon.fr/)

