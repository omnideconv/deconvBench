#!/usr/bin/env Rscript
'
remapCelltypesNF.R
Remap celltypes according to mapping sheet to controlled vocaulary.
Usage:
  remapCelltypesNF.R <mapping_sheet> <celltype_vector>
Options:
  <mapping_sheet>     An an excel document with a "mapping" sheet containing the celltype mapping.
  <celltype_vector>   Vector with celltypes that should be remapped.
' -> doc

#library(docopt)
#arguments = docopt(doc)

##### How to remap celltypes
# get obs$cell_type
# name it input
# replace by real celltypes
### might have to add all this to excel sheet as well --> dropbox: wget + link mit lesezugriff 

#read input -> for singlecell remap celltype annotations
#with rna seq remap facs

###FOR TESTING###
test <- '
remapCelltypesSheet = readxl::read_excel("/nfs/data/omnideconv_benchmarking/cell_type_mapping.xlsx", sheet="mapping")
input <- sample(method_celltypes$method_cell_type, 50, replace = TRUE) ###use real celltypes here###
celltype_annotations <- readRDS("/nfs/data/omnideconv_benchmarking/data/singleCell/maynard/celltype_annotations.rds")
input <- celltype_annotations
method_ds <- "maynard"
'

#remappingSheet = An an excel document with a "mapping" sheet containing the celltype mapping.
#celltype_annotations = Vector with celltypes that should be remapped.
#method_ds = name of dataset (according to mapping sheet)

remapCelltypesWorkflow <- function(remappingPath, celltype_annotations, method_ds){
  remapCelltypesSheet = readxl::read_excel(remappingPath, sheet="mapping")
  input <- celltype_annotations
  if(method_ds %in% unique(remapCelltypesSheet$method_dataset)){
    method_celltypes <- remapCelltypesSheet[remapCelltypesSheet$method_dataset == method_ds,]
    celltypes <- method_celltypes$cell_type[match(input, method_celltypes$method_cell_type)]
  } else {
    print(method_ds)
    print(remapCelltypesSheet$method_dataset)
  }
  result = data.frame(input=input, converted=celltypes)
  celltype_annotations <- celltypes
  return(celltype_annotations)
}

remapCelltypesTree <- function(facs_celltypes, deconv_celltypes){
  #remap eg Tcell CD4 and CD8 etc to the higher level that matches the facs
  df <- data.frame(parent_type=c(), child_type=c())
  for (celltype in facs_celltypes) {
    children <- immunedeconv::get_all_children(cell_type = celltype)
    df <- rbind(df, data.frame(parent_type=rep(celltype, length(children)), child_type=children))
  }
  #what if it's not child but child of child etc?
  #how are we gonna deal with this --> sometimes celltypes just not in facs - need to distinguish if not finished mapping or no parent available
  if(all(deconv_celltypes %in% df$child_type)){
    print(deconv_celltypes)
    print(df$child_type)
  }
  return(df)
}