# SCALLOP
Annotation of results from the SCALLOP consortium INF1 GWAS.

Files are located in `/home/jm2294/rds/rds-jmmh2-projects/olink_proteomics/scallop/jm2294` on CSD3.

`inf1_phenoscanner_corrplot.r` will pull the conditionally independent results from the INF GWAS and query phenoscanner. This then generates a file `efo_list.csv` that is manually annotated (sorry!) to tag immune-related clinical phenotypes. The annotated file `efo_list_annotated.csv` is then read in and used to generate the heatmap matrix. 

The heatmap matrix is then again manually anntoated (sorry again!) to merge together traits which are listed separately in phenoscanner but represent the same thing, e.g. 'Gout' and 'Self-Reported Gout'.

The final part of the code then generates the heatmap.
