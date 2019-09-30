library(tidyverse)
library(scales)
library(tidylog)
library(egg)
library(ggrastr)


explode_hail_hist = function(data, edge_col, freq_col, above_col) {
  data %>%
    mutate(edge=str_split(str_remove_all(.data[[edge_col]], '[\\[\\]]'), ','),
           data=str_split(str_remove_all(.data[[freq_col]], '[\\[\\]]'), ','),
           data=pmap(list(data, .data[[above_col]]), c)
    ) %>%
    unnest(c('edge', 'data')) %>%
    mutate(edge=as.numeric(edge),
           data=as.numeric(data))
}

explode_hail_column = function(data, col) {
  data %>%
    mutate(data=str_split(str_remove_all(.data[[col]], '[\\[\\]]'), ',')) %>%
    unnest(c('data')) %>%
    mutate(data=as.numeric(data))
}

overall_ab_hist = function(data_type, common_indels_only=T, snps_only=F) {
  hom_data = read_delim(paste0('data/groupings_', data_type, '.txt'), delim='\t')
  
  if (common_indels_only) {
    p = hom_data %>%
      explode_hail_column("overall_ab_hist") %>%
      select(freq, length_group, data) %>%
      filter(freq == 0.2, length_group == -10) %>%
      mutate(x = row_number() / 20) %>% 
      mutate(data = data / sum(data)) %>%
      ggplot + aes(x = x, y = data) +
      geom_bar(stat='identity') + theme_classic() +
      xlab('AB for hets (deletions >= 10 bp, >= 20% frequency)') + ylab('Proportion of het genotypes')
  } else if (snps_only) {
    p = hom_data %>%
      explode_hail_column("overall_ab_hist") %>%
      select(freq, length_group, data) %>%
      filter(freq == 0.01, length_group == 0) %>%
      mutate(x = row_number() / 20) %>% 
      mutate(data = data / sum(data)) %>%
      ggplot + aes(x = x, y = data) +
      geom_bar(stat='identity') + theme_classic() +
      xlab('AB for hets (SNVs, 1-5% frequency)') + ylab('Proportion of het genotypes')
  } else {
    p = hom_data %>%
      explode_hail_column("overall_ab_hist") %>%
      select(freq, length_group, data) %>%
      group_by(freq, length_group) %>%
      mutate(x = row_number() / 20) %>% 
      mutate(data = data / sum(data)) %>%
      ungroup %>%
      ggplot + aes(x = x, y = data) + facet_grid(length_group ~ freq) +
      geom_bar(stat='identity') + theme_classic() +
      xlab('AB for hets') + ylab('Proportion of het genotypes')
  }
}

count_hom_undercalls_by_freq_and_category = function(data_type) {
  data_type = 'exomes'
  hom_data = read_delim(paste0('data/groupings_', data_type, '.txt'), delim='\t')
  
  hom_data %>%
    explode_hail_hist("added_hom_hist.bin_edges", "added_hom_hist.bin_freq", "added_hom_hist.n_larger") %>%
    select(freq, length_group, edge, data) %>%
    group_by(freq, length_group) %>%
    mutate(data = data / sum(data)) %>% ungroup %>%
    ggplot + aes(x = edge, y = data) + facet_grid(length_group ~ freq) +
    geom_bar(stat='identity') + theme_classic() +
    xlab('Number of hets with AB > 0.9') + ylab('Proportion of variants')
  
  hom_data %>%
    explode_hail_hist("prop_added_hom_hist.bin_edges", "prop_added_hom_hist.bin_freq", "prop_added_hom_hist.n_larger") %>%
    select(freq, length_group, edge, data) %>%
    group_by(freq, length_group) %>%
    mutate(data = data / sum(data)) %>% ungroup %>%
    filter(edge < 0.1) %>%
    ggplot + aes(x = edge, y = data) + facet_grid(length_group ~ freq) +
    geom_bar(stat='identity') + theme_classic() +
    xlab('Proportion of hets with AB > 0.9') + ylab('Proportion of variants')
}

figure1 = function() {
  p1 = overall_ab_hist('genomes')
  p2 = overall_ab_hist('genomes', common_indels_only = F, snps_only = T)
  
  f2 = ggarrange(p1, p2, labels = c('a', 'b'), label.args=list(gp=grid::gpar(font=2, cex=1.2)))
  png('figure1.png', height=5, width=4, res=300, units = 'in')
  print(f2)
  dev.off()
}


hom_undercall_summary = function(data_type='genomes', proportion=T, save_plot=F) {
  hom_data = read_delim(paste0('data/groupings_', data_type, '.txt'), delim='\t')
  
  hom_data %>% 
    summarize(total_variants_with_added_hom=sum(total_with_added_hom),
              total_variants_with_1pct_added_hom=sum(total_with_1pct_added_hom), 
              total_variants_with_2pct_added_hom=sum(total_with_2pct_added_hom), 
              total_variants_with_5pct_added_hom=sum(total_with_5pct_added_hom), 
              total_variants=sum(total))
  
  hom_data %>%
    group_by(length_group) %>%
    summarize(total_variants_with_added_hom=sum(total_with_added_hom),
              total_variants_with_1pct_added_hom=sum(total_with_1pct_added_hom), 
              total_variants=sum(total))
  
  if (proportion) {
    p = hom_data %>%
      mutate(freq=fct_inseq(as.factor(freq)),
             prop_with_added_hom=total_with_added_hom / total,
             prop_with_1pct_added_hom=total_with_1pct_added_hom / total) %>%
      ggplot + aes(x = length_group, y = prop_with_1pct_added_hom, color = freq, group = freq) + 
      geom_line() + geom_point() +
      xlab('Variant length category') + ylab('Percent of variants with\nat least 1% hom increase\ndue to AB > 0.9 hets') + 
      theme_classic() + 
      scale_y_continuous(labels=percent_format(accuracy=1)) +
      scale_color_hue(h=c(0, -120) + 15, name='Frequency',
                      labels=function(x) {percent(as.numeric(x), accuracy = 1)}) +
      guides(color = guide_legend(reverse = TRUE))
  } else {
    p = hom_data %>%
      mutate(freq=fct_inseq(as.factor(freq))) %>%
      ggplot + aes(x = length_group, y = added_homs / called_homs, color = freq) +
      geom_line() + geom_point() +
      xlab('Variant length category') + ylab('Mean increase in\nhomozygote count') + 
      theme_classic() + 
      scale_y_continuous(labels=percent_format(accuracy=1)) +
      scale_color_hue(h=c(0, -120) + 15, name='Frequency',
                      labels=function(x) {percent(as.numeric(x), accuracy = 1)})
  }
  
  if (save_plot) {
    pdf(paste0('hom_undercall_', ifelse(proportion, 'proportion_variants_', 'mean_increase_'), data_type, '.pdf'),
        width=6, height=4)
    print(p)
    dev.off()
  }
  return(p)
}

figure2 = function() {
  p1 = hom_undercall_summary(data_type='genomes', proportion=T, save_plot=F)
  p2 = hom_undercall_summary(data_type='genomes', proportion=F, save_plot=F)
  
  f2 = ggarrange(p1, p2, labels = c('a', 'b'), label.args=list(gp=grid::gpar(font=2, cex=1.2)))
  png('figure2.png', height=5, width=6, res=300, units = 'in')
  print(f2)
  dev.off()
}


reclassification_summary = function(data_type='genomes') {
  reclass_data = read_delim(paste0('data/contamination_summary_', data_type, '.txt'), delim='\t')
  
  reclass_data %>%
    filter(freq == 0.2) %>%
    mutate(mean_reclassified_per_variant = mean_number_reclassified / n_variants) %>%
    ggplot + aes(x = length_group, y = mean_reclassified_per_variant, color = segdup) +
    geom_point() + geom_line() + 
    scale_y_continuous(labels=percent_format(accuracy=0.1)) +
    xlab('Variant length category') + 
    ylab('Mean proportion of genotypes affected per variant') + 
    theme_classic()
  
  reclass_data %>%
    filter(!segdup) %>%
    mutate(freq=fct_inseq(as.factor(freq)),
           mean_reclassified_per_variant = mean_number_reclassified / n_variants) %>%
    ggplot + aes(x = length_group, y = mean_reclassified_per_variant, 
                 color = freq) +
    scale_color_hue(h=c(0, -120) + 15, name='Frequency',
                    labels=function(x) {percent(as.numeric(x), accuracy = 1)}) +
    scale_y_continuous(labels=percent_format(accuracy=0.1)) +
    geom_point() + geom_line() + 
    xlab('Variant length category') + 
    ylab('Mean proportion of genotypes affected per variant') + 
    theme_classic() +
    guides(color = guide_legend(reverse = TRUE))
}

contamination_vs_reclassification = function(data_type='genomes', metric='number_reclassified', save_plot=F) {
  contam_data_summary = read_delim(gzfile(paste0('data/num_reclassified_by_sample_all_variants_', data_type, '.txt.bgz')), delim='\t')
  
  p = contam_data_summary %>%
    ggplot + aes(x = freemix) + aes_string(y = metric) +
    geom_point_rast() + theme_classic() +
    xlab('Contamination') + 
    ylab(ifelse(metric == 'number_reclassified', 
                'Number of AB > 0.9 het calls', 'Number of AD[0] > 0 hom alt calls'))
  
  if (save_plot) {
    pdf(paste0('contamination_v_', metric, '_', data_type, '.pdf'),
        width=6, height=4)
    print(p)
    dev.off()
  }
  return(p)
}

figure3 = function() {
  p = contamination_vs_reclassification()
  
  png('figure3.png', height=4, width=6, res=300, units='in')
  print(p)
  dev.off()
}


contamination_summary = function(data_type='genomes') {
  contam_data = read_delim(paste0('data/contamination_summary_', data_type, '.txt'), delim='\t')
  
  contam_data %>%
    filter(!segdup) %>%
    mutate(freq=fct_inseq(as.factor(freq))) %>%
    ggplot + aes(x = length_group, y = linreg_number_reclassified.multiple_r_squared, 
                 color = freq) +
    scale_color_hue(h=c(0, -120) + 15, name='Frequency',
                    labels=function(x) {percent(as.numeric(x), accuracy = 1)}) +
    geom_point() + geom_line() + 
    scale_y_continuous(labels=percent_format(accuracy=1)) +
    xlab('Variant length category') + 
    ylab('R-squared between contamination\nand number of AB > 0.9 het calls') + 
    theme_classic() +
    guides(color = guide_legend(reverse = TRUE))
}

