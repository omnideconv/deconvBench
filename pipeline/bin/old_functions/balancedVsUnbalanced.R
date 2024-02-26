deconvolution_uniform_lambrechts_bisque <- readRDS("/nfs/proj/omnideconv_benchmarking/benchmark/pipeline/results/deconvolution_uniform_lambrechts_bisque.rds")
deconvolution_bisque_lambrechts_finotello <- readRDS("/nfs/proj/omnideconv_benchmarking/benchmark/pipeline/results/deconvolution_bisque_lambrechts_finotello.rds")

df <- as.data.frame(deconvolution_bisque_lambrechts_finotello) %>% tibble::rownames_to_column(var="sample") %>%
  pivot_longer(!"sample", names_to="celltype", values_to = "predicted_value") %>% mutate(data="original")
df2 <- as.data.frame(deconvolution_uniform_lambrechts_bisque) %>% tibble::rownames_to_column(var="sample") %>%
  pivot_longer(!"sample", names_to="celltype", values_to = "predicted_value") %>% mutate(data="uniform")
full <- rbind(df, df2)
ggplot(full, aes(data, predicted_value, fill=celltype))+geom_boxplot()+facet_wrap(~celltype)

##that does not make sense it's about the single cell dataset not bulk!