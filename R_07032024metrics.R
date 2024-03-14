happy_12x_summury = read.table("/home/elodie/Documents/Analysis/Happy_12x/Happy_12x/ref-12x.summary.csv", header=TRUE, sep=',')
happy_24x_summury = read.table("/home/elodie/Documents/Analysis/Happy_24x/Happy_24x/ref-24x.summary.csv", header=TRUE, sep=',')
happy_40x_summury = read.table("/home/elodie/Documents/Analysis/Happy_40x/Happy_40x/ref-40x.summary.csv", header=TRUE, sep=',')

happy_12x_extended = read.table("/home/elodie/Documents/Analysis/Happy_12x/Happy_12x/ref-12x.extended.csv", header=TRUE, sep=',')
happy_24x_extended = read.table("/home/elodie/Documents/Analysis/Happy_24x/Happy_24x/ref-24x.extended.csv", header=TRUE, sep=',')
happy_40x_extended = read.table("/home/elodie/Documents/Analysis/Happy_40x/Happy_40x/ref-40x.extended.csv", header=TRUE, sep=',')



library('ggplot2')
library('tools')

read_single = function(x) {
  nx = strsplit(x, "\\:")[[1]]
  
  if(length(nx) == 1) {
    name = basename(file_path_sans_ext(x))
  } else {
    x = nx[1]
    name = nx[2]
  }
  cat(sprintf("Reading %s as %s\n", x, name))
  
  all_results = list()
  print(all_results)
  all_results$roc_data_snp_all = read.csv(paste(x, "roc.Locations.SNP", "csv", "gz", sep="."))
  all_results$roc_data_indel_all = read.csv(paste(x, "roc.Locations.INDEL", "csv", "gz", sep="."))
  
  all_results$roc_data_snp_pass = read.csv(paste(x, "roc.Locations.SNP.PASS", "csv", "gz", sep="."))
  all_results$roc_data_indel_pass = read.csv(paste(x, "roc.Locations.INDEL.PASS", "csv", "gz", sep="."))
  print(all_results)
  sel_snp_file = paste(x, "roc.Locations.SNP.SEL", "csv", "gz", sep=".")
  
  # if we have a selectively-filtered ROC, don't show "PASS" ROC
  if(file.exists(sel_snp_file)) {
    all_results$roc_data_snp_sel = read.csv(sel_snp_file)
  } else {
    all_results$roc_data_snp_sel = all_results$roc_data_snp_pass
  }
  
  all_results$roc_data_snp_all = head(subset(all_results$roc_data_snp_all,
                                             QQ == min(all_results$roc_data_snp_all["QQ"])),
                                      n=1)
  print(all_results$roc_data_snp_all)
  all_results$roc_data_snp_pass = head(subset(all_results$roc_data_snp_pass,
                                              QQ == min(all_results$roc_data_snp_pass["QQ"])),
                                       n=1)
  sel_min = head(subset(all_results$roc_data_snp_sel,
                        QQ == min(all_results$roc_data_snp_sel["QQ"])),
                 n=1)
  
  all_results$roc_connector_snp =
    rbind(all_results$roc_data_snp_all, sel_min)
  all_results$roc_connector_snp$Filter = "CONN"
  print(all_results$roc_connector_snp)
  all_results$roc_data_snp_sel$Filter = "ROC"
  print(all_results$roc_data_snp_sel)
  sel_indel_file = paste(x, "roc.Locations.INDEL.SEL", "csv", "gz", sep=".")
  if(file.exists(sel_indel_file)) {
    all_results$roc_data_indel_sel = read.csv(sel_indel_file)
  } else {
    # use PASS ROC if no SEL ROC present
    all_results$roc_data_indel_sel = all_results$roc_data_indel_pass
  }
  
  # just keep single ALL and PASS point
  all_results$roc_data_indel_all = head(subset(all_results$roc_data_indel_all,
                                               QQ == min(all_results$roc_data_indel_all["QQ"])),
                                        n=1)
  all_results$roc_data_indel_pass = head(subset(all_results$roc_data_indel_pass,
                                                QQ == min(all_results$roc_data_indel_pass["QQ"])),
                                         n=1)
  
  sel_min = head(subset(all_results$roc_data_indel_sel,
                        QQ == min(all_results$roc_data_indel_sel["QQ"])),
                 n=1)
  
  all_results$roc_connector_indel =
    rbind(all_results$roc_data_indel_all, sel_min)
  all_results$roc_connector_indel$Filter = "CONN"
  all_results$roc_data_indel_sel$Filter = "ROC"
  
  result = do.call(rbind, all_results)
  row.names(result) = NULL
  print(class(result))
  
  result$filename = x
  result$name = name
  result$igroup = paste(result$name,
                        result$Filter,
                        result$Type)
  return(result)
}

vcf12x = read_single("/home/elodie/Documents/Analysis/Happy_12x/Happy_12x/ref-12x")
vcf24x = read_single("/home/elodie/Documents/Analysis/Happy_24x/Happy_24x/ref-24x")
vcf40x = read_single("/home/elodie/Documents/Analysis/Happy_40x/Happy_40x/ref-40x")

#head(vcf12x)

# Plot P/R curves
plot_data = function(pdata, is.PR=FALSE) {
  # precision / recall curve
  if(is.PR) {
    xaxis = "METRIC.Recall"
    yaxis = "METRIC.Precision"
  } else {
    # approximate ROC-style curve (FPR is not correct)
    xaxis = "FPR"
    yaxis = "TPR"
    pdata$FPR = pdata$QUERY.FP / (pdata$QUERY.TOTAL - pdata$QUERY.UNK)
    pdata$TPR = pdata$TRUTH.TP / (pdata$TRUTH.TP + pdata$TRUTH.FN)
    cc = complete.cases(pdata[, c(xaxis, yaxis)])
    print(pdata)
    pdata = pdata[cc, ]
    print(pdata)
  }
  
  plt = ggplot(pdata, aes_string(x=xaxis, y=yaxis, color="name"))
  facet_wrap(~Type)
  
  # ROC lines
  plt = plt +
    geom_line(data = subset(pdata, Filter == "ROC"),
              mapping=aes(group=igroup),
              size=1.5,
              linetype=3)
  
  # Connector between ALL and start of ROC
  plt = plt +
    geom_line(data = subset(pdata, Filter == "CONN"),
              mapping=aes(group=igroup),
              size=1.5,
              linetype=4)
  
  plt = plt +
    geom_point(data = subset(pdata, Filter %in% c("CONN")),
               mapping=aes(group=igroup),
               size=8)
  
  plt = plt +
    geom_point(data = subset(pdata, Filter %in% c("PASS", "ALL")),
               mapping=aes(shape=Filter, group=igroup),
               size=8)
  
  xl_min = max(0,
               min(subset(pdata, Filter %in% c("PASS", "ALL"))[[xaxis]]) - 0.02)
  xl_max = min(1.0,
               max(subset(pdata, Filter %in% c("PASS", "ALL"))[[xaxis]]) + 0.02)
  yl_min = max(0,
               min(subset(pdata, Filter %in% c("PASS", "ALL"))[[yaxis]]) - 0.01)
  yl_max = min(1.0,
               max(subset(pdata, Filter %in% c("PASS", "ALL"))[[yaxis]]) + 0.01)
  
  plt = plt +
    scale_color_brewer("", palette="Set1") +
    #xlim(c(xl_min, xl_max)) +
    xlim(c(0, 1)) +
    #ylim(c(yl_min, yl_max)) +
    ylim(c(0, 1)) +
    theme_bw(base_size=18)
  
  return(plt)
}

vcf12x_indel = filter(vcf12x, Type=="INDEL")
vcf12x_snp = filter(vcf12x, Type=="SNP")
#plot indel+snp_12x
plot_data (vcf12x, is.PR=TRUE)
#plot indel_12x
plot_data (vcf12x_indel, is.PR=TRUE)
#plot snp_12x
plot_data (vcf12x_snp, is.PR=TRUE)

vcf24x_indel = filter(vcf24x, Type=="INDEL")
vcf24x_snp = filter(vcf24x, Type=="SNP")
plot_data (vcf24x, is.PR=TRUE)
plot_data (vcf24x_indel, is.PR=TRUE)
plot_data (vcf24x_snp, is.PR=TRUE)

vcf40x_indel = filter(vcf40x, Type=="INDEL")
vcf40x_snp = filter(vcf40x, Type=="SNP")
plot_data (vcf40x, is.PR=TRUE)
plot_data (vcf40x_indel, is.PR=TRUE)
plot_data (vcf40x_snp, is.PR=TRUE)

vcf_12x_24x_40x_indel = rbind(vcf12x_indel, vcf24x_indel,vcf40x_indel)
vcf_12x_24x_40x_snp = rbind(vcf12x_snp, vcf24x_snp, vcf40x_snp)
plot_data (vcf_12x_24x_40x_indel, is.PR=FALSE)
plot_data (vcf_12x_24x_40x_indel, is.PR=TRUE)
plot_data (vcf_12x_24x_40x_snp, is.PR=FALSE)
plot_data (vcf_12x_24x_40x_snp, is.PR=TRUE)