library(ggplot2)
library(dplyr)

df = data.frame(multi_map_parameter = c(10, 100, 150, 200))

df$multi_map_parameter = as.factor(df$multi_map_parameter)

df$annotated_reads = c(71792494.63372958, 71993044.47464405, 71996200.89465106, 71998061.01108)
df$multi_reads = c(8845166, 9349924, 9355569, 9360374)

df = mutate(df, proportion_multi_reads = (multi_reads / annotated_reads))

plot = ggplot(data = df, aes(x = multi_map_parameter, y = proportion_multi_reads)) +
  geom_bar(stat = 'identity', color = 'black') +
  scale_y_continuous(labels = scales::percent, expand = expansion(mult = c(0, .1))) +
  xlab('outFilterMultimapNmax_array') +
  ylab('Proportion of multi-mapped reads') +
  ggtitle('STAR parameter optimization', 'Patient 226 mTEC-HI sample')

plot + theme_bw() + theme(plot.title = element_text(face = 'bold', size = 20),
                               plot.subtitle = element_text(size = 14),
                               panel.grid.major.x = element_blank(),
                               panel.grid.minor.x = element_blank(),
                               panel.grid.major.y = element_line(color = 'grey'),
                               panel.grid.minor.y = element_blank(),
                               axis.text.x = element_text(size = 13),
                               axis.text.y = element_text(size = 14),
                               axis.title = element_text(size = 14),
                               axis.line = element_line(size = 0.8),
                               panel.border = element_blank(),
                               legend.text = element_text(size = 12),
                               legend.title = element_text(size = 14))