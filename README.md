# DRPPM-PATH-SURVEIOR-Suite

# Introduction

The integration of patient genome expression data, phenotype data, and clinical data can serve as an integral resource for patient prognosis. The DRPPM PATH SURVEIOR Suite: **Path**way level **Surv**ival **E**xam**i**nat**or** serves to do just that, by examining the interaction of gene or gene set pathway expression with cilinical data to discover prominent features that play a role in patient outcome. This utility is comprised of 3 R Shiny apps and a pipeline script which can be employed in a cohesive manor to provide an in-depth analysis towards underlying features affective survival. Through a systematic Cox proportional regression pipeline the expression of individual genes or gene set pathways is linked with a hazard ratio and p.values to determine signifigance and patient risk, which can be followed by interactive visualization in the suite of R Shiny applications.

Users can download and employ all of the tools locally to perform analyses and gain an in-depth perspective of their own data. To start, we provide example skin and kidney cancer data from the PAN ICI (Immune Checkpoint Inhibition) iAtlas database. Additionally, we provied a comprehensive set of pathways which havew been derieved from a variety of databases, such as MSigDB, LINCS L1000, Cell Marker, as well as derived gene sets focusing on ER stress and immune signatures. The example data and provided gene sets will be referenced and utilized throughout the guided overview of the DRPPM-PATH-SURVEIOR Suite of tools. 

## The DRPPM-PATH-SURVEIOR Family

* R Shiny Base Survival App [Interactive Mode]: https://github.com/shawlab-moffitt/DRPPM-PATH-SURVEIOR-Suite/tree/main/2-DRPPM-PATH-SURVEIOR-InteractiveApp
* R Script for Cox Proportional Hazards Ranking [Pipeline Mode]: https://github.com/shawlab-moffitt/DRPPM-PATH-SURVEIOR-Suite/tree/main/3-DRPPM-PATH-SURVEIOR-Pipeline
* R Shiny Jaccard Connectivity App: https://github.com/shawlab-moffitt/DRPPM-PATH-SURVEIOR-Suite/tree/main/4-DRPPM-Jaccard_Connectivity_App
* R Shiny Pre-Ranked GSEA App: https://github.com/shawlab-moffitt/DRPPM-PATH-SURVEIOR-Suite/tree/main/5-DRPPM-Hazard_Ratio_Ranked_GSEA_App

![alt text](https://github.com/shawlab-moffitt/DRPPM-PATH-SURVEIOR-Suite/blob/main/2-DRPPM-PATH-SURVEIOR-InteractiveApp/App_Demo_Pictures/mian_schematic.PNG?raw=true)




# Disclamer

Copyright 2022 Moffitt Cancer Center
Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:
1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.
2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.
3. Neither the name of the copyright holder nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission.
THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
