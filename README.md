# ShinySDAscSeq (ShinSASe) -shin-sa-say (審査せ)

## Introduction

This is a template for a ShinyApp to be used for processing scRNASeq data processed with SDA.

### Input: 

ShinyServerDataLS: a list of data:

"datat" a dataframe with cell barcode + meta data including 1 dimreduc (x,y) like tsne. The other paramters differ per study so custom.
"results" : an SDA object post processing and label correction holds the cell scores and gene loading
"chromosome.lengths" : a dataframe with "chromosome", "length", "length_padded", "genomic_offset", "center" from mapping genes to biomart
"StatFac": A dataframe of "SDA" comp and annotations found in processing, labels to be shown when a comp is selected
"GO_data": a list of GO annotations, 1 per every component in the resuts
"gene_locations": from biomart
"col_vector" : vector of colors to be used
