#Here are some functions which enhance nanotail package

################################################################################
#                   GETMODE
################################################################################

# function to get mode value (based on: https://stackoverflow.com/questions/2547402/how-to-find-the-statistical-mode)
## added as a separate function not within another function - just in case
getmode <- function(x, method ="density", na.rm = FALSE) {
  x <- unlist(x)
  if (na.rm) {
    x <- x[!is.na(x)]
  }
  
  if (method %in% c("value", "density", "") | is.na(method)) {
    # Return actual mode (from the real values in dataset)
    if (method %in% c("density", "")) {
      
      # Return density mode for normalized data - only for numeric!)
      d <- density(x)
      return(d$x[d$y==max(d$y)]) 
      #return(modeest::mlv(x,na.rm=TRUE,method="parzen", abc=T)) #in some cases this method produces weird output
      
    } else if (method %in% c("value")) {
      
      uniqx <- unique(x)  
      n <- length(uniqx)
      freqx <- tabulate(match(x, uniqx))
      modex <- freqx == max(freqx)
      return(uniqx[which(modex)])
    }
  }
}  

################################################################################
#                   SUMMARIZE_POLYA2
################################################################################
#slightly changed nanotails summarize function - now it also returns most frequent values (mode)
summarize_polya2 <- function(polya_data,summary_factors = c("group"),transcript_id_column = c("transcript")) {
  polya_data_summarized2 <-
    polya_data %>% dplyr::ungroup() %>% dplyr::group_by(.dots = c(transcript_id_column,summary_factors)) %>% dplyr::summarise(
      counts = dplyr::n(),
      polya_mean = mean(polya_length),
      polya_sd = sd(polya_length),
      polya_median = median(polya_length),
      polya_gm_mean = gm_mean(polya_length),
      polya_sem = polya_sd/sqrt(counts),
      polya_mode = getmode(polya_length, "value"),
      polya_dens_mode = getmode(polya_length, "density"))
  
  return(polya_data_summarized2)
}


################################################################################
#                   ANNOTATE_WITH_BIOMART_2
################################################################################
#Slightly changed nanotails' function - drop biotype attribute, as it causes problems (duplicates) in further analysis - this code is for editing diffexp data only!

annotate_with_biomart2 <- function(polya_data,attributes_to_get=c('external_gene_name','description', 'ensembl_gene_id'),filters='ensembl_gene_id',mart_to_use=NA) {
  
  ensembl_ids <- unique(polya_data$ensembl_gene_id)
  ensembl_ids <- ensembl_ids[!is.na(ensembl_ids)]
  
  #polya_data <- polya_data %>% dplyr::rename(ensembl_transcript_id = ensembl_transcript_id_short)
  
  {
    annotation_data<-biomaRt::getBM(attributes=attributes_to_get, filters = filters, values = ensembl_ids, mart = mart_to_use)
  }
  polya_data_annotated <-  polya_data %>% dplyr::left_join(annotation_data)
  
  return(polya_data_annotated)
}



################################################################################
#                   PLOT_DISTRIBUTION_WITH_VALUES
################################################################################
#This modified version of plot_distribution function has options to show or hide center values; it is intended to use in every density plots across the document
#usage: plot_distribution_with_values(input_data, value_to_show="mode", factor="group", number = T, title= "Mitochondrial transcripts", origin="mouse")
#This modified version of plot_distribution function has options to show or hide center values
#Based on PK's snippets

plot_distribution_with_values <- function(input_data,value_to_show="mode", grouping_factor = "group", title="", number = T, origin="", limit="") {
  
  
  #determine the summary level factor
  if (grouping_factor=="group") {
    factor_label <- grouping_factor
  }
  else if (grouping_factor=="sample_name") {
    factor_label <- paste0("sample")
  }
  else {
    stop("Wrong grouping factor specified")
  }
  #determine the center values to be plotted as x intercepting line(s)
  center_values = input_data %>% dplyr::group_by(!!rlang::sym(grouping_factor)) %>% dplyr::summarize(median_value = median(polya_length,na.rm = TRUE), mode_value=getmode(polya_length,na.rm=TRUE,method="density"),mean_value=mean(polya_length,na.rm=T))
  
  #call nanotail
  plot_distribution<-nanotail::plot_polya_distribution(input_data,show_center_values = "none",groupingFactor = grouping_factor,color_palette="nejm") + coord_cartesian(ylim = c(0, 1.05), xlim = c(0, as.numeric(limit)))
  
  if (value_to_show=="median") {
    center_value="median_value"
  }
  else if (value_to_show=="mode") {
    center_value="mode_value"
  }
  else if (value_to_show=="mean") {
    center_value="mean_value"
  }
  else if (value_to_show=="") {
    center_value=""
  }
  else {
    stop("Wrong center value to show specified (should be median, mode or mean")
  }
  
  
  if (origin=="") {
    p= paste0("")
  }
  else if (origin!="") {
    p= paste0("Source: ", origin)
  }
  else {
    stop("An unexpected error regarding origin argument ocurred. Strange.")
  }
  
  #Plot settings (aesthetics, geoms, axes behavior etc.):
  g.line <- ggplot2::geom_vline(data=center_values,aes(xintercept=!!rlang::sym(center_value),color=!!rlang::sym(grouping_factor)),linetype="longdash")
  g.labs <- ggplot2::labs(title= paste0("Poly(A) lengths distribution per ", factor_label),
                          subtitle= paste0(title),
                          x="poly(A) length",
                          y= "normalized density ",
                          caption = paste0("Dashed lines - ", value_to_show, " values\n", p),
                          color=factor_label)
  g.values <- ggrepel::geom_text_repel(data=center_values,aes(x=round(!!rlang::sym(center_value)),y=1.07,color=!!rlang::sym(grouping_factor),label=formatC(round(!!rlang::sym(center_value)),digits=1,format = "d")),size=4, force=1, nudge_x = 2, direction = "x", hjust= 0, show.legend =F)
  
  #Overall plotting configuration:
  plot <- plot_distribution + mytheme + g.line + g.labs 
  
  if (number==TRUE) {
    plot_values<-plot + g.values 
  }
  else if (number==FALSE) {
    plot_values<-plot
  }
  else {
    stop("Specify whether to show values or not (number= F or T)")
  }
  
  return(plot_values)
}


################################################################################
#                   PLOT_BOX_WITH_VALUES
################################################################################
#THIS ONE works for plotting distribution boxplots for data grouped per group and sample, respectively. It allows for plotting with or without numeric values.
#This is a bit time consuming, equip with patience
plot_box_with_values <- function(input_data, grouping_factor = "group", value_to_show ="mode", title="", origin="", limit="", number=T) {
  
  #determine the summary level factor
  if (grouping_factor=="group") {
    factor_label <- grouping_factor
    facet <- facet_grid(. ~ group,drop = T,scales = "free_x")
  }
  else if (grouping_factor=="sample_name") {
    factor_label <- paste0("sample")
    facet <- facet_grid(. ~ sample_name,drop = T,scales = "free_x")
  }
  else {
    stop("Wrong grouping factor specified")
  }
  
  if (value_to_show=="median") {
    center_value="median_value"
    
  }
  else if (value_to_show=="mode") {
    center_value="mode_value"
    
  }
  else if (value_to_show=="mean") {
    center_value="mean_value"
    
  }
  else if (value_to_show=="") {
    center_value=""
    
  }
  else {
    stop("Wrong center value to show specified (should be median, mode or mean")
  }
  
  
  if (origin=="") {
    p= paste0("")
  }
  else if (origin!="") {
    p= paste0("Source: ", origin)
  }
  else {
    stop("An unexpected error regarding origin argument ocurred. Strange.")
  }
  
  #mutate the initial dataset:
  data <- input_data %>% dplyr::group_by(!!rlang::sym(grouping_factor)) %>% dplyr::add_count() %>% dplyr::ungroup() %>% dplyr::mutate(label = paste0("(n=",n,")"))
  
  #determine the center values to be plotted as intercepting line(s):
  center_values <- input_data %>% dplyr::group_by(!!rlang::sym(grouping_factor)) %>% dplyr::summarize(median_value = median(polya_length,na.rm = TRUE), mode_value=getmode(polya_length,na.rm=TRUE,method="density"),mean_value=mean(polya_length,na.rm=T))
  # merge mutated data with computed center values: 
  final_data <- dplyr::left_join(data, center_values, by=grouping_factor)
  
  
  #Plot settings (aesthetics, geoms, axes behavior etc.):
    g.plot <- ggplot2::ggplot(final_data, aes(x=label,y=polya_length, fill=group)) + coord_cartesian(ylim = c(0, as.numeric(limit))) + ggsci::scale_fill_simpsons(alpha=0.8) + scale_x_discrete() 
    g.violin <- ggplot2::geom_violin(aes())
    g.boxplot <- ggplot2::geom_boxplot(width=0.2,fill="white",notch = TRUE,outlier.size = 0)
    g.labs <- ggplot2::labs(title= paste0("Poly(A) lengths distribution per ", factor_label),
                            subtitle= paste0(title),
                            x="\npoly(A) length",
                            y= "normalized density ",
                            caption = paste0("Dashed lines - ", value_to_show, " values\n", p),
                            color=factor_label)
    g.line <- ggplot2::geom_hline(data=final_data,aes(yintercept=!!rlang::sym(center_value)),linetype="longdash", size = 0.6, show.legend =F, color ="#db0d14")
    g.values <- ggplot2::geom_text(data = final_data, aes(1, !!rlang::sym(center_value), label = formatC(round(!!rlang::sym(center_value)),digits=1,format = "d")), vjust = 2, hjust = 1, color="#db0d14", nudge_x = 0.6, check_overlap = TRUE)
    
  
  if (number==TRUE) {
    plot <- g.plot + g.violin + g.boxplot + mytheme + g.labs + g.line + g.values + facet 
  }
  else if (number==FALSE) {
    plot <- g.plot + g.violin + g.boxplot + mytheme + g.labs + g.line + facet
  }
  else {
    stop("Specify whether to show values or not (number= F or T)")
  }
  
  return(plot)
  
}



################################################################################
#                   PLOT_LENGTH_LOG2FOLDCHANGE
################################################################################
# Plotting log2fold change vs poly(A) length diff:

##this function plots log2fold change plots for poly(A) length gm means and modes. There are 2 colouring codes available: potential TENT substrates and significantly differently expressed genes.

plot_length_log2foldchange <- function(input_data, x="length_diff", origin="", color="") {
  
  #determine whether the x would be gm or mode
  if (x=="length_diff") {
    x_label <- "(gm)"
  }
  else if (x=="mode_length_diff") {
    x_label <- "(mode)"
  }
  else {
    stop("Wrong x variable provided")
  }  
  
  if (origin=="") {
    p= paste0("")
  }
  else if (origin!="") {
    p= paste0("Source: ", origin)
  }
  else {
    stop("An unexpected error regarding origin argument ocurred. Strange.")
  }
  
  #determine coloring variable
  if (color=="TENT") {
    g.point <- ggplot2::geom_point(aes(color = !!rlang::sym(x) <=-5, size = effect_size), alpha =0.5)
    g.color <- ggplot2::scale_colour_manual(name = "Transcript type:",values = c('darkgrey', '#d98818'), labels = c("other", "TENT5C substrate"), na.translate=FALSE)
    
  }
  else if (color=="diffexp") {
    g.point <- ggplot2::geom_point(aes(color = diffexp_sig, size = effect_size), alpha =0.5)
    g.color <- ggplot2::scale_colour_manual(name = "Differential expression:", values = c('#08a3a8', 'darkgrey'), na.translate=FALSE)
    
  }
  else if (color=="") {
    color_value=""
    
  }
  else {
    stop("Wrong value to be coloured specified (should be TENT, diffexp or nothing")
  }
  #PLOT CONFIGURATION:
  g.plot <- ggplot2::ggplot(input_data, aes(x = !!rlang::sym(x), y = log2FoldChange))
  g.vline <- ggplot2::geom_vline(xintercept = 0, linetype="dashed", color = "#db0d14") 
  g.hline <- ggplot2::geom_hline(yintercept = 0, linetype="dashed", color = "#db0d14")
  g.labs <- ggplot2::labs(title= paste0("Expression change vs. difference of poly(A) length ", x_label),
                          x=paste0("\u0394 poly(A) length ", x_label, " between CKO and WT"),
                          y= "Effect size: " ~log[2]~"(fold-change) in expression (nanopore)",
                          caption = p,
                          size="Effect size (tail length):")
  
  plot <- g.plot + g.point + g.vline + g.hline + theme_bw() + g.color + g.labs
  
  return(plot)  
}

################################################################################
#                   PLOT_LENGTH_LOG2FOLDCHANGE
################################################################################
#splitviolin - code of 'jan-glx' was obtained from Stackoverflow https://stackoverflow.com/questions/35717353/split-violin-plot-with-ggplot2

GeomSplitViolin <- ggproto(
  "GeomSplitViolin",
  GeomViolin,
  draw_group = function(self, data, ..., draw_quantiles = NULL) {
    data <- transform(data,xminv = x - violinwidth * (x - xmin),xmaxv = x + violinwidth * (xmax - x)
    )
    grp <- data[1, "group"]
    newdata <- plyr::arrange(transform(data, x = if (grp %% 2 == 1) xminv else xmaxv),
                             if (grp %% 2 == 1) y else -y
    )
    newdata <- rbind(newdata[1, ], newdata, newdata[nrow(newdata), ], newdata[1, ])
    newdata[c(1, nrow(newdata) - 1, nrow(newdata)), "x"] <- round(newdata[1, "x"])
    if (length(draw_quantiles) > 0 & !scales::zero_range(range(data$y))) {
      stopifnot(all(draw_quantiles >= 0), all(draw_quantiles <= 1))
      quantiles <- ggplot2:::create_quantile_segment_frame(data, draw_quantiles)
      aesthetics <- data[rep(1, nrow(quantiles)), setdiff(names(data), c("x", "y")), drop = FALSE]
      aesthetics$alpha <- rep(1, nrow(quantiles))
      both <- cbind(quantiles, aesthetics)
      quantile_grob <- GeomPath$draw_panel(both, ...)
      ggplot2:::ggname(
        "geom_split_violin",
        grid::grobTree(GeomPolygon$draw_panel(newdata, ...), quantile_grob)
      )
    } else {
      ggplot2:::ggname("geom_split_violin", GeomPolygon$draw_panel(newdata, ...))
    }
  }
)

geom_split_violin <- function(mapping = NULL,data = NULL,stat = "ydensity",position = "identity", ..., draw_quantiles = NULL,trim = TRUE,scale = "area",na.rm = FALSE,show.legend = NA,inherit.aes = TRUE) {
  layer(data = data,
        mapping = mapping,
        stat = stat,
        geom = GeomSplitViolin,
        position = position,
        show.legend = show.legend,
        inherit.aes = inherit.aes,
        params = list(
          trim = trim,
          scale = scale,
          draw_quantiles = draw_quantiles,
          na.rm = na.rm, ...)
  )
}


################################################################################
#             PLOT_ABUNDANT_VIOPLOT
################################################################################

#function to plot split violin plot with boxes - most abundant
plot_abundant_vioplot <- function(abundant_to_box) {
  input_data <- polya_data_filtered_annotated %>% group_by(group) %>% add_count(transcript) %>% filter(dense_rank(-n) <=15)
  g.plot <- ggplot(input_data, aes(x = reorder(transcript, -n), y = polya_length, stat=transcript, fill=group)) + coord_cartesian(ylim = c(0, 100)) + ggsci::scale_fill_d3(alpha = 0.5)
  g.violin <- geom_split_violin(aes(trim=FALSE)) 
  g.labels <- labs(x="\ntranscripts", y="poly(A) length\n", title = "Most abundant transcripts", caption = "Red lines - mode values \nSource: murine pancreas")
  g.theme <- theme_bw() + theme(axis.text.x = element_text(angle = 90, size = 10, vjust=0.5, hjust=1 ), 
                                axis.text.y = element_text(size = 10), 
                                axis.title=element_text(size=10, face="bold"), 
                                legend.title = element_text(size = 10), 
                                legend.text = element_text(size = 10))
  #g.boxplot <- geom_boxplot(width = 0.25, notch = FALSE, outlier.shape = NA, coef=0, show.legend =FALSE)
  g.mode <- stat_summary(fun = getmode, fun.min = getmode, fun.max = getmode, show.legend =FALSE, geom = "crossbar", color = "red", width = 0.25, position = position_dodge(width = .25),)
  final.plot <- g.plot + g.violin + g.labels + g.theme + g.mode 
  plot(final.plot)
}
plot_abundant_vioplot()



################################################################################
#             PLOT_LONG_VIOPLOT
################################################################################

plot_long_vioplot <- function(abundant_to_box) {
  input_data <- polya_data_filtered_annotated  %>% dplyr::inner_join(head(summarized_polya_table_per_transcript %>%  dplyr::filter(counts >= 100) %>% dplyr::arrange(desc(polya_mean)), 21))
  g.plot <- ggplot(input_data, aes(x = reorder(transcript, -polya_mean), y = polya_length, stat=transcript, fill=group)) + coord_cartesian(ylim = c(0, 150)) + ggsci::scale_fill_d3(alpha = 0.5)
  g.violin <- geom_split_violin() 
  g.labels <- labs(x="\n transcripts", y="poly(A) length\n", title = "Transcripts with longest poly(A) tails", caption = "Red lines - mode values \nSource: murine pancreas")
  g.theme <- theme_bw() + theme(axis.text.x = element_text(angle = 90, size = 10 ), 
                                axis.text.y = element_text(size = 10), 
                                axis.title=element_text(size=10, face="bold"), 
                                legend.title = element_text(size = 10), 
                                legend.text = element_text(size = 10))
  #g.boxplot <- geom_boxplot(width = 0.2, notch = FALSE, outlier.shape = NA, coef=0, show.legend =FALSE)
  g.mode <- stat_summary(fun = getmode, fun.min = getmode, fun.max = getmode, show.legend =FALSE, geom = "crossbar", color = "red", width = 0.25, position = position_dodge(width = .25),)
  final.plot <- g.plot + g.violin + g.labels + g.theme + g.mode #+ g.boxplot 
  plot(final.plot)
}
plot_long_vioplot()



################################################################################
#             RAW_TO_SQUIGGLE
################################################################################
#This function plots raw signal (squiggle) with segmentation based on nanopolish output
# It plots a single signal of interest based on its position within the polya_data output
# to run this, one must provide the workspace with basecalled .fast5 files & sequencing summary within config.yml (most convenient way) or otherwise 
# this is only a lookup function; not the most efficient, but it works

raw_to_squiggle <- function(row = NA){
  
  #polya_table <- read_polya_multiple(samples_table = input_samples_table) - in case if
  #sequencing_summary<-vroom::vroom(input_samples_table[1,]$seq_summ,col_select=c("filename","read_id"))
  
  selected_row <- row #numerical - position within the polya output table
  selected_read <- polya_table[selected_row,]$read_id #this is a read name without an extension (.fast5)
  filename <-sequencing_summary[sequencing_summary$read_id==selected_read,]$filename #this is a read name with an extension (.fast5)
  selected_workspace <- polya_table[selected_row,]$workspace # it is actually a dir to fast5s (in workspace) turned into a factor
  filename2 <- paste0(selected_workspace,'/',filename) # it is actually a fast5 directory
  read_id2<-paste0("read_",selected_read) # this has to be done, bcause of the convention within fast5
  
  hdf5_file <- hdf5r::h5file(filename2)
  raw_read <- hdf5_file[[read_id2]][["Raw"]]
  raw_signal <- raw_read[["Signal"]]
  no_events <- raw_signal$dims
  signal_df <- data.frame(pos=seq(1,no_events,1),signal=raw_signal[1:no_events]) # this is signal converted to dframe actually -> this is PK mode of action
  hdf5_file$close_all()
  
  polya_start <- polya_table[selected_row,]$polya_start
  transcript_start <- polya_table[selected_row,]$transcript_start
  adapter_start <- polya_table[selected_row,]$adapter_start
  leader_start <- polya_table[selected_row,]$leader_start
  read_id_title <- polya_table[selected_row,]$read_id
  transcript_title <- polya_table[selected_row,]$transcript
  
  polya_start_pos <- polya_table[selected_row,18] %>% pull(polya_start)
  polya_end_pos <- polya_table[selected_row,19] %>% pull(transcript_start)
  
  #adaptor sequence
  signal_df$segment[signal_df$pos > polya_start_pos] <- "transcript"
  #transcript sequence
  signal_df$segment[signal_df$pos < polya_end_pos +1] <- "adapter"
  #polya sequence
  signal_df$segment[signal_df$pos>polya_start_pos -1 & signal_df$pos<polya_end_pos +1] <- "poly(A)"
  signal_df$segment <- as.factor(signal_df$segment)
  
  
  
  plot <- ggplot2::ggplot(data=signal_df, aes(x=pos, y=signal, color=segment)) + geom_line() + geom_vline(xintercept = polya_start, color = "red") + geom_vline(xintercept = transcript_start, color = "blue") + theme_bw() + labs(title= paste0("Raw signal of ", selected_read, " read")) + scale_colour_manual(values = c('#089bcc', "#f56042",'#3a414d'))
  
  plotly <- plotly::ggplotly(plot)
  
  return(plot)
  
}









