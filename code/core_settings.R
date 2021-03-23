###################
## Core Settings ##
###################

library(ggplot2)
library(viridis)
library(patchwork)
library(ggrepel)
library(ggrastr)
 
library(readr)
library(dplyr)
library(tibble)

library(knitr)
library(kableExtra)

# set different default parameters
library(reticulate)
use_python("/usr/bin/python3")

###################
## Plot Settings ##
###################

library(ggplot2)
library(viridis)


theme_publication = theme_bw() +
                      theme(text=element_text(size=5), axis.text=element_text(size=5), axis.title=element_text(size=6), plot.title=element_text(size=6, hjust=0.5)) +
                      theme(strip.text= element_text(size=6)) +
                      theme(legend.title = element_text(size=5), legend.text = element_text(size=5)) +
                      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
                      theme(plot.tag=element_text(size=8, face="bold")) +
                      # remove facet box
                      theme(strip.background = element_blank())

theme_presentation = theme_bw() +
                      theme(text=element_text(size=14), axis.text=element_text(size=12), axis.title=element_text(size=14), plot.title=element_text(size=14, hjust=0.5)) +
                      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
                      theme(plot.tag=element_text(size=16, face="bold")) +
                      # remove facet box
                      theme(strip.background = element_blank())

celltype_colours = c(# neuronal
                    "Embryonic RG" = "#635547",
                    "Juvenile RG & TAPs" = "#DABE99",
                    "Gliogenic precursors (aNSCs)" = "#8DB5CE",
                    "Gliogenic precursors (APCs)" = "#3F84AA",
                    "Quiescent NSCs (dorsal)" = "#005579",
                    "Quiescent NSCs (ventral)" = "#C9EBFB",
                    # "Quiescent NSCs [1]" = "#005579",
                    # "Quiescent NSCs [2]" = "#C9EBFB",
                    "Early EmNBs" = "#9E6762",
                    "GE NBs [1]" = "#C19F70",
                    "GE NBs [2]" = "#139992",
                    "Hippocampal NBs [1]" = "#C594BF",
                    "Hippocampal NBs [2]" = "#0F4A9C",
                    "OPCs" = "#B51D8D",
                    "PreM-ODCs" = "#532C8A",
                    # "ImPreMDs" = "#532C8A",
                    "Oligodendrocytes [1]" = "#8870AD",
                    "Oligodendrocytes [2]" = "#CC7818",
                    "Ependymal cells" = "#FBBE92",
                    "GABAergic INs" = "#FACB12",
                    "ImStNeurons" = "#F397C0",
                    
                    # vascular/immune/blood
                    "VLMC" = "#F6BFCB",
                    "Juvenile microglia" = "#F9DECF",
                    "Adult microglia [1]" = "#C72228",
                    "Adult microglia [2]" = "#F79083",
                    "Macrophages" = "#EF4E22",
                    "T-cells" = "#8EC792",
                    "Erythrocytes" = "#CDE088",
                    "Choroid plexus epithelia" = "#F7F79E",
                    "Endothelial cells" = "#647A4F",
                    "Myeloid-DSCs" = "#354E23"
)



# NSC_RG_lineage_colours = c("Embryonic RG"="#532C8A",
#                            "Juvenile RG"="#3F84AA",
#                            "Active NSCs"="#65A83E",
#                            "aNSCs"="#65A83E",
#                            "Quiescent NSCs [1]"="#DABE99",
#                            "Quiescent NSCs (dorsal)"="#DABE99",
#                            "Quiescent NSCs"="#9F8A70",
#                            "qNSCs"="#9F8A70",
#                            "Quiescent NSCs (ventral)"="#635547",
#                            "Quiescent NSCs [2]"="#635547")

neuronal_glial_colours = c("Embryonic RG" = "#C72228",
                          "Juvenile RG & TAPs" = "#FF891C",
                          "Gliogenic precursors (aNSCs)" = "#FACB12",
                          "Gliogenic precursors (APCs)" = "#F9DECF",
                          "Ependymal cells" = "#139992",
                          "Quiescent NSCs [1]" = "#DABE99",
                          "Quiescent NSCs (dorsal)" = "#DABE99",
                          "Quiescent NSCs [2]" = "#635547",
                          "Quiescent NSCs (ventral)" = "#635547",
                          "Early EmNBs [1]" = "#005579",
                          "Early EmNBs [2]" = "#8DB5CE",
                          "GABAergic INs" = "#CDE088",
                          "GE NBs [1]" = "#3F84AA",
                          "GE NBs [2]" = "#F7F79E",
                          "ImStNeurons" = "#B51D8D",
                          "Hippocampal NBs [1]" = "#C594BF",
                          "Hippocampal NBs [2]" = "#F397C0",
                          "OPCs" = "#C9EBFB",
                          "ImPreMDs" = "#CC7818",
                          "PreM-ODCs" = "#CC7818",
                          "Oligodendrocytes [1]" = "#8870AD",
                          "Oligodendrocytes [2]" = "#532C8A")




NSC_RG_lineage_colours = c("Embryonic RG"="#C72228",
                          "Juvenile RG"="#FF891C",
                          "Active NSCs"="#FACB12",
                          "aNSCs"="#8DB5CE",
                          "Quiescent NSCs [1]"="#DABE99",
                          "Quiescent NSCs (dorsal)"="#DABE99",
                          "Quiescent NSCs"="#9F8A70",
                          "qNSCs"="#9F8A70",
                          "Quiescent NSCs (ventral)"="#635547",
                          "Quiescent NSCs [2]"="#635547")




embryonic_RG_colours = c("Pretectum"="#FACB12",
                         "Epithalamus"="#B51D8D",
                         "Thalamic eminence"="#65A83E",
                         "Gliogenic (ganglionic)"="#139992",
                         "Cortical pallium"="#DABE99",
                         "Cortical hem"="#EF5A9D",
                         "Gliogenic (cortical)"="#FF891C",
                         "Subthalamic nucleus"="#C72228",
                         "Ganglionic eminence"="#354E23")

colour_timepoints = c("E12.5"="#fdcc8a",
                     "E14.5"="#fc8d59",
                     "E16.5"="#d7301f",
                     "P0"="#bdc9e1",
                     "P2"="#74a9cf",
                     "P7"="#0570b0",
                     "P13"="#d9d9d9",
                     "P19"="#bdbdbd",
                     "P39"="#969696",
                     "P111"="#636363",
                     "P365"="#252525")