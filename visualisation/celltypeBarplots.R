celltype_annotations_lambrechts <- readRDS("/nfs/data/omnideconv_benchmarking/data/singleCell/lambrechts/celltype_annotations.rds")
celltype_annotations_lambrechts <- remapCelltypesWorkflow(remappingPath = remapping_sheet, celltype_annotations = celltype_annotations_lambrechts, method_ds = "lambrechts")
ggplot(as.data.frame(table(celltype_annotations_lambrechts)))+
  geom_col(aes(celltype_annotations_lambrechts, Freq, fill=celltype_annotations_lambrechts))+
  theme(axis.text.x = element_text(angle=90), legend.position = "none")+
  labs(x="cell type", y="frequency", title = "Lambrechts")
ggsave("barplotCelltypes_lambrechts.jpeg", width = 7, height = 4)

celltype_annotations_maynard <- readRDS("/nfs/data/omnideconv_benchmarking/data/singleCell/maynard/celltype_annotations.rds")
celltype_annotations_maynard <- remapCelltypesWorkflow(remappingPath = remapping_sheet, celltype_annotations = celltype_annotations_maynard, method_ds = "maynard")
ggplot(as.data.frame(table(celltype_annotations_maynard)))+
  geom_col(aes(celltype_annotations_maynard, Freq, fill=celltype_annotations_maynard))+
  theme(axis.text.x = element_text(angle=90), legend.position = "none")+
  labs(x="cell type", y="frequency", title = "Maynard")
ggsave("barplotCelltypes_maynard.jpeg", width = 7, height = 4)
