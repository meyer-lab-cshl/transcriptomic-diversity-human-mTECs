library(ggplot2)

input = results_df_local
input = mutate(input, ID = sub("\\?", "", ID))
input = mutate(input, class = sub("\\?", "", class))
input = filter(input, significant == TRUE)
input_random_subset = sample_n(input, 50)

##input = filter(input, family == 'ERV1')
#input = filter(input, class == 'Satellite')

input = mutate(input, signed_log10padj = -log10(padj) * sign(log2FoldChange))

## Plot signed p value

ggplot(data = input, aes(x = gene, y = signed_log10padj)) +
  geom_point(alpha = 0.6, size = 0.5) +
  facet_grid(. ~ class)

## Plot fold change

ggplot(data = input, aes(x = gene, y = log2FoldChange)) +
  geom_point() +
  facet_grid(. ~ family)
  
  
  #geom_hline(yintercept = -log10(0.05), linetype = 'dashed') +
  #geom_hline(yintercept = log10(0.05), linetype = 'dashed')

