#NINTH

## Prior to running this R markdown, you must run DCC on your samples to test for circular RNA. See https://github.com/dieterich-lab/DCC for instructions on how to perform this using Unix.

## I have also uploaded my own DCC log on the GitHub repository for this project.

```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE)
```

```{r}
#install.packages("devtools")
library("devtools")
require(devtools)
#install.packages("remotes")
#install_github('dieterich-lab/CircTest')
library("CircTest")
library("utils")
library("dplyr")
library("ggplot2")
library("pheatmap")
library("RColorBrewer")
library("tidyverse")
library("rstatix")
library("stats")
library("data.table")
library("DEXSeq")
library("ggpubr")
#install_github("jdstorey/qvalue")
library("qvalue")
library("aod")
library("rtracklayer")
library("GenomicRanges")
library("gdata")
```

# Make a list of your cell types and replicate numbers:
```{r}
singletypes <- c("ADL","AFD","AIM","AIN","AIY","ASEL","ASER","ASG","ASI","ASK","AVA","AVE","AVG","AVH","AVK","AVL","AVM","AWA","AWB","AWC","BAG","CAN","DA","DD","DVC","I5","IL1","IL2","LUA","NSM","OLL","OLQ","PHA","PVC","PVD","PVM","RIA","RIC","RIM","RIS","RMD","SMB","SMD","VB","VC","VD")

singletypes_withreps <- c("ADL1", "ADL2", "ADL3", "ADL4", "AFD1", "AFD2", "AFD3", "AFD4", "AFD5", "AIM1", "AIM2", "AIM3", "AIM4", "AIN1", "AIN2", "AIN3", "AIN4", "AIN5", "AIN6", "AIY1", "AIY2", "AIY3", "ASEL1", "ASEL2", "ASEL3", "ASER1", "ASER2", "ASER3", "ASER4", "ASG1", "ASG2", "ASG3", "ASG4", "ASI1", "ASI2", "ASI3", "ASI4", "ASK1", "ASK2", "ASK3", "ASK4", "AVA1", "AVA2", "AVA3", "AVA4", "AVA5", "AVA6", "AVE1", "AVE2", "AVE3", "AVG1", "AVG2", "AVG3", "AVH1", "AVH2", "AVH3", "AVH4", "AVK1", "AVK2", "AVK3", "AVK4", "AVL1", "AVL2", "AVL3", "AVM1", "AVM2", "AVM3", "AWA1", "AWA2", "AWA3", "AWA4", "AWB1", "AWB2", "AWB3", "AWB4", "AWB5", "AWC1", "AWC2", "AWC3", "AWC4", "BAG1", "BAG2", "BAG3", "BAG4", "CAN1", "CAN2", "CAN3", "DA1", "DA2", "DA3", "DA4", "DD1", "DD2", "DD3", "DD4", "DVC1", "DVC2", "DVC3", "DVC4", "I51", "I52", "I53", "I54", "IL11", "IL12", "IL13", "IL21", "IL22", "IL23", "IL24", "LUA1", "LUA2", "LUA3", "LUA4", "NSM1", "NSM2", "NSM3", "OLL1", "OLL2", "OLQ1", "OLQ2", "OLQ3", "PHA1", "PHA2", "PHA3", "PHA4", "PVC1", "PVC2", "PVC3", "PVC4", "PVC5", "PVD1", "PVD2", "PVM1", "PVM2", "RIA1", "RIA2", "RIA3", "RIA4", "RIA5", "RIA6", "RIC1", "RIC2", "RIC3", "RIC4", "RIM1", "RIM2", "RIM3", "RIM4", "RIS1", "RIS2", "RIS3", "RMD1", "RMD2", "RMD3", "RMD4", "RMD5", "SMB1", "SMB2", "SMB3", "SMB4", "SMB5", "SMD1", "SMD2", "SMD3", "SMD4", "VB1", "VB2", "VB3", "VB4", "VC1", "VC2", "VC3", "VC4", "VC5", "VC6", "VD1", "VD2", "VD3", "VD4")

##NOTE: Since OLL2 and PVC1 had no circRNA detected (No circRNA passed the expression threshold filtering), no files could be generated and renamed using the DCC scripts
```

# Import of DCC output files into R:
```{r, include=FALSE}
setwd("D:/Zach Wolfe's CircRNA analysis") # set your working directory here

singletypes_withreps <- singletypes_withreps[!(singletypes_withreps %in% c("OLL1", "OLL2", "PVC1"))] # OLL2 and PVC1 had no circRNA detected, which means they will have no output files for DCC, which means we have to exclude them from our analysis
## additionally, since there are only 2 reps of OLL, this means we will have to eliminate OLL entirely from our circRNA analysis at this point. Otherwise we will have a problem with Circ.test in our next chunk
singletypes <- singletypes[!(singletypes %in% "OLL")]

circRNA_table <- matrix(nrow = 177, ncol = 2, dimnames = list(x = singletypes_withreps, y = c(paste("circRNA count"), paste("Most frequently circularized gene"))))
circRNA_table[is.na(circRNA_table)] <- 0

CircRNACount_filtered_combined <- data.frame(matrix(ncol = 3, nrow = 0))
colnames <- c("Chr", "Start", "End")
colnames(CircRNACount_filtered_combined) <- colnames
CircRNACount_filtered_combined$Chr <- as.character()
CircRNACount_filtered_combined$Start <- as.numeric()
CircRNACount_filtered_combined$End <- as.numeric()

LinearCount_filtered_combined <- data.frame(matrix(ncol = 3, nrow = 0))
colnames <- c("Chr", "Start", "End")
colnames(LinearCount_filtered_combined) <- colnames
LinearCount_filtered_combined$Chr <- as.character()
LinearCount_filtered_combined$Start <- as.numeric()
LinearCount_filtered_combined$End <- as.numeric()

CircCoordinates_filtered_combined <- data.frame(matrix(ncol = 8, nrow = 0))
colnames <- c("Chr", "Start", "End", "Gene", "JunctionType", "Strand", "Start.End.Region", "OverallRegion")
colnames(CircCoordinates_filtered_combined) <- colnames
CircCoordinates_filtered_combined$Chr <- as.character()
CircCoordinates_filtered_combined$Start <- as.numeric()
CircCoordinates_filtered_combined$End <- as.numeric()
CircCoordinates_filtered_combined$Gene <- as.character()
CircCoordinates_filtered_combined$JunctionType <- as.numeric()
CircCoordinates_filtered_combined$Strand <- as.character()
CircCoordinates_filtered_combined$Start.End.Region <- as.character()
CircCoordinates_filtered_combined$OverallRegion <- as.character()

for (i in 1:length(singletypes_withreps)){
CircRNACount <- read.delim(paste0(singletypes_withreps[i], 'CircRNACount'))
LinearCount <- read.delim(paste0(singletypes_withreps[i], 'LinearCount'))
CircCoordinates <- read.delim(paste0(singletypes_withreps[i], 'CircCoordinates'))

CircRNACount_filtered <- Circ.filter(circ = CircRNACount,
                                     linear = LinearCount,
                                     Nreplicates = 1,
                                     filter.sample = 1,
                                     filter.count = 5, # default is 5
                                     percentage = 0, # default is 0.1
                                    )


CircCoordinates_filtered <- CircCoordinates[rownames(CircRNACount_filtered),]
LinearCount_filtered <- LinearCount[rownames(CircRNACount_filtered),]

# We need to combine chimeric.out.junctions columns, then put them on one big table

names(CircRNACount_filtered)[names(CircRNACount_filtered) == "Chimeric.out.junction"] <- paste0("Chimeric.out.junction", singletypes_withreps[i])
CircRNACount_filtered_combined <- full_join(CircRNACount_filtered_combined, CircRNACount_filtered)

names(LinearCount_filtered)[names(LinearCount_filtered) == "Chimeric.out.junction"] <- paste0("Chimeric.out.junction", singletypes_withreps[i])
LinearCount_filtered_combined <- full_join(LinearCount_filtered_combined, LinearCount_filtered)

CircCoordinates_filtered_combined <- full_join(CircCoordinates_filtered_combined, CircCoordinates_filtered)


circRNA_table[i,1] <- nrow(CircRNACount_filtered)

if (circRNA_table[i,1] != 0){
gene_count <- as.data.frame(table(CircCoordinates_filtered$Gene))
gene_count <- gene_count[rev(order(gene_count$Freq)),]
circRNA_table[i,2] <- as.character(first(gene_count$Var1))
 }
if (circRNA_table[i,1] == 0){
circRNA_table[i,2] <- as.character(paste("NA"))
 }
}

# For now, we need NA values to be zero in order to perform Circ.test later.
CircRNACount_filtered_combined[is.na(CircRNACount_filtered_combined)] <- 0
LinearCount_filtered_combined[is.na(LinearCount_filtered_combined)] <- 0
CircCoordinates_filtered_combined[is.na(CircCoordinates_filtered_combined)] <- 0
```


# Circular RNA test (Circ.test):
```{r, message=FALSE}
group <- c(rep(1:length(singletypes), times = sapply(singletypes, function(x) sum(grepl(paste0("^", x), singletypes_withreps)))))
group <- as.numeric(group, levels = unique(group))

test <- Circ.test(CircRNACount_filtered_combined,
                 LinearCount_filtered_combined,
                 CircCoordinates_filtered_combined,
                 group=group, # I have 45 groups with variable numbers of replicates in each group
                 alpha = 1 # p-value cutoff (default is 0.05)
)
```

# Significant results are shown in a summary table:
```{r}
colnames(test$summary_table) <- c("Chr", "Start", "End", "Gene", "JunctionType", "Strand", "Start.End.Region", "OverallRegion", "sig_p", paste0(singletypes, "_ratio_mean"))

#test$summary_table[is.na(test$summary_table)] <- 0

# p-value cutoff:
test_p_summary_table <- test$summary_table[test$summary_table$sig_p < 0.05,]

# optional: re-order by row names:
test_summary_table_ordered <- test$summary_table[order(as.numeric(rownames(test$summary_table))),]

View(test$summary_table)

# unless EVERY SINGLE replicate of a given cell type has a value present for a given RNA strand, the Circ.ratioplot won't be able to plot ANY ratios for that cell type. Here is my solution:
## Change all 0s in LinearCount file to 0.00001

LinearCount_filtered_combined_0_as_fraction <- LinearCount_filtered_combined
LinearCount_filtered_combined_0_as_fraction[LinearCount_filtered_combined == 0] <- 0.00001
```

# Make a summarySE function (this is from DCC's github):
```{r}
summarySE <- function(data=NULL, measurevar, groupvars=NULL, na.rm=FALSE,
                      conf.interval=.95, .drop=TRUE) {
  require(plyr)
  
  # New version of length which can handle NA's: if na.rm==T, don't count them
  length2 <- function (x, na.rm=FALSE) {
    if (na.rm) sum(!is.na(x))
    else       length(x)
  }
  
  # This does the summary. For each group's data frame, return a vector with
  # N, mean, and sd
  datac <- ddply(data, groupvars, .drop=.drop,
                 .fun = function(xx, col) {
                   c(N    = length2(xx[[col]], na.rm=na.rm),
                     mean = mean   (xx[[col]], na.rm=na.rm),
                     sd   = sd     (xx[[col]], na.rm=na.rm)
                   )
                 },
                 measurevar
  )
  
  # Rename the "mean" column    
  datac <- rename(datac, c("mean" = measurevar))
  
  datac$se <- datac$sd / sqrt(datac$N)  # Calculate standard error of the mean
  
  # Confidence interval multiplier for standard error
  # Calculate t-statistic for confidence interval: 
  # e.g., if conf.interval is .95, use .975 (above/below), and use df=N-1
  ciMult <- qt(conf.interval/2 + .5, datac$N-1)
  datac$ci <- datac$se * ciMult
  
  return(datac)
}
```

# Edit DCC's Circ.ratioplot function (this will allow us to manually adjust parameters and aesthetics):
```{r}
my.Circ.ratioplot <- function(Circ,Linear,CircCoordinates = None,plotrow='1',size=12,ncol=2,groupindicator1=NULL,groupindicator2=NULL,x='Conditions',y='circRNA/(circRNA+Linear)',lab_legend='groupindicator1', circle_description = c(1:3), gene_column = None, y_axis_range = 1, colour_mode = "colour"){

  if( !is.null(groupindicator1) & length(groupindicator1) != ncol(Circ)-length(circle_description) ){
    stop("If provided, the length of groupindicator1 should be equal to the number of samples.")
  }
  if( !is.null(groupindicator2) & length(groupindicator2) != ncol(Circ)-length(circle_description) ){
    stop("If provided, the length of groupindicator2 should be equal to the number of samples.")
  }
  if(is.null(groupindicator1)){
    stop("At least one grouping should be provided through groupindicator1.")
  }
  if(!is.null(groupindicator2)){
    twolevel <- TRUE
  }else{
    twolevel <- FALSE
  }
  
  rownames.circ <- rownames(Circ)
  Circ <- data.frame(lapply(Circ, as.character), stringsAsFactors=FALSE)
  rownames(Circ) <- rownames.circ
  
  rownames.linear <- rownames(Linear)
  Linear <- data.frame(lapply(Linear, as.character), stringsAsFactors=FALSE)
  rownames(Linear) <- rownames.linear
  
  if(!missing(CircCoordinates)){
    rownames.circCoordinates <- rownames(CircCoordinates)
    CircCoordinates <- data.frame(lapply(CircCoordinates, as.character), stringsAsFactors=FALSE)
    rownames(CircCoordinates) <- rownames.circCoordinates
  }else{
    CircCoordinates <- data.frame(Circ[,circle_description])
    rownames(CircCoordinates) <- rownames.circ
    rownames.circCoordinates <- rownames(CircCoordinates)
    CircCoordinates <- data.frame(lapply(CircCoordinates, as.character), stringsAsFactors=FALSE)
    rownames(CircCoordinates) <- rownames.circCoordinates  
  }
  
  groupindicator1 <- factor(groupindicator1,levels=unique(groupindicator1))
  groupindicator2 <- factor(groupindicator2,levels=unique(groupindicator2))
  
  # Get gene name, if no annotation, output NULL
  if (is.character(plotrow)){
    if ( ! plotrow %in% rownames(CircCoordinates) ){
      stop("Specified 'plotrow' not found.")
    }
  }else{
    if ( is.numeric(plotrow) ){
      if ( ! plotrow %in% 1:nrow(CircCoordinates) ){
        stop("Specified 'plotrow' not found.")
      }
    }else{
      stop("Specified plotrow should be ONE rowname or ONE rownumber.")
    }
  }
  # Choose your own column containing the gene name using gene_column. The genename will be displayed in the plot title if available
  if (missing(gene_column)){
    genename = NULL
  }else{
    genename <- as.character(CircCoordinates[plotrow,gene_column])
    if (genename == '.'){
      genename = NULL
    }
  }
  if(twolevel){
    plotdat <- summarySE( data.frame(Ratio=as.numeric(Circ[plotrow,-circle_description])/(as.numeric(Linear[plotrow,-circle_description])+as.numeric(Circ[plotrow,-circle_description])),
                                    groupindicator1,
                                    groupindicator2),
                         measurevar='Ratio',groupvars=c('groupindicator1','groupindicator2') )
  }else{
    plotdat <- summarySE( data.frame(Ratio=as.numeric(Circ[plotrow,-circle_description])/(as.numeric(Linear[plotrow,-circle_description])+as.numeric(Circ[plotrow,-circle_description])),
                                     groupindicator1),
                                     measurevar='Ratio',groupvars=c('groupindicator1') )
  }
# construct plot
  Q <- ggplot(plotdat, aes(x=groupindicator1, y=Ratio)) +
       geom_boxplot() + theme_classic() +
       theme(axis.text.x = element_blank())+
       theme(axis.text.y = element_text(size=8))+
       theme(axis.ticks = element_line(colour = 'black', linewidth = 0.7)) +
       theme(axis.ticks.x = element_blank())+
       theme(legend.title=element_blank()) + 
       theme(text=element_text(size=9))+
       theme(legend.text=element_text(size=8)) +
       theme(plot.title = element_text(size=10)) + 
       theme(axis.text.y = element_text(margin=margin(5,5,10,5,"pt")))+
       ggtitle(paste("Annotation: ", genename, "\nChr ", toString(Circ[plotrow,circle_description]),sep="")) +
       ylab("Percentage of Circular RNA vs All RNA") + 
       xlab("Cell Type") +
       geom_bar(stat="identity",aes(fill=groupindicator1), color = "black", size=1) +         
       geom_errorbar(aes(ymin=Ratio-se, ymax=Ratio+se), width=1 , size=0.5)

  if (colour_mode == "bw"){
      Q <- Q + scale_fill_grey(start = 0.0, end = 1)
  } else {
      Q <- Q + scale_fill_discrete(name=lab_legend)
  }

       Q <- Q +
       theme(legend.position="bottom") +
       theme(axis.ticks.length = unit(0.4, "cm")) +
       theme(panel.background = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        panel.border = element_rect(colour = "black", fill=NA, size=2)) + 
        guides(fill=guide_legend(
                 keywidth=0.3,
                 keyheight=0.2,
                 default.unit="inch")
      ) + scale_y_continuous(expand=c(0,0), limits= c(0, y_axis_range), labels = scales::percent_format(scale = 100))

  if(twolevel){
    Q <- Q + facet_wrap( ~ groupindicator2,ncol=ncol )
  }

  print(Q)
}
```


# Read in an updated and correct gtf file (from WormBase):
```{r}
gtf <- fread("caenorhabditis_elegans.PRJNA13758.WBPS18.canonical_geneset.gtf")
colnames(gtf) <- c("chr", "source", "feature", "start", "end", "score", "strand", "frame", "attribute")
# replace the combined attribute column with multiple columns that are more appropriate for naming:
gtf$gene_id <- sub('.*gene_id "(.*?)".*', '\\1', gtf$attribute)
gtf$gene_version <- sub('.*gene_version "(.*?)".*', '\\1', gtf$attribute)
gtf$gene_biotype <- sub('.*gene_biotype "(.*?)".*', '\\1', gtf$attribute)
gtf$gene_name <- sub('.*gene_name "(.*?)".*', '\\1', gtf$attribute)

# Remove the original 'attribute' column:
gtf <- gtf[, c("chr", "source", "feature", "start", "end", "score", "strand", "frame", "gene_id", "gene_version", "gene_biotype", "gene_name")]
```

# Extract the coordinates of each Circ event. This will be used to produce the correct and updated gene names for CircRNA coordinates with more than one gene name:
```{r}
for (i in seq_len(nrow(CircCoordinates_filtered_combined))) {
  
  # Extract chromosome, strand, start, and end coordinates from the Coordinates data frame:
  chr_var <- CircCoordinates_filtered_combined$Chr[i]
  strand_var <- CircCoordinates_filtered_combined$Strand[i]
  start_var <- as.numeric(CircCoordinates_filtered_combined$Start[i])
  end_var <- as.numeric(CircCoordinates_filtered_combined$End[i])

  # Find the corresponding gene in the GTF dataframe:
  matching_gene <- subset(gtf, chr == chr_var & feature == "gene" & start <= start_var & end >= end_var)
  
  if (nrow(matching_gene) == 0){
    matching_gene_start <- subset(gtf, chr == chr_var & feature == "gene" & start <= start_var)
    matching_gene_end <- subset(gtf, chr == chr_var & feature == "gene" & end >= end_var)
    matching_gene_start <- matching_gene_start[nrow(matching_gene_start),]
    matching_gene_end <- matching_gene_end[1,]
    if (matching_gene_start$end >= start_var & matching_gene_end$start <= end_var){
      gene_name <- paste0(matching_gene_start$gene_name, "; ", matching_gene_end$gene_name)
    }
    else if (matching_gene_start$end >= start_var & matching_gene_end$start >= end_var){
      gene_name <- matching_gene_start$gene_name
    }
    else if (matching_gene_start$end <= start_var & matching_gene_end$start <= end_var){
      gene_name <- matching_gene_end$gene_name
    }
    else {
      gene_name <- "not_annotated"
    }

    # Replace the original gene name with the updated one:
    CircCoordinates_filtered_combined$Gene[i] <- gene_name
  }
  
  if (nrow(matching_gene) > 0){
    gene_name <- matching_gene$gene_name

    # Replace the original gene name with the updated one:
    CircCoordinates_filtered_combined$Gene[i] <- gene_name
  }
}

for (i in seq_len(nrow(CircRNACount_filtered_combined))) {
  
  # Extract chromosome, start, and end coordinates from the data frame:
  chr_var <- CircRNACount_filtered_combined$Chr[i]
  start_var <- as.numeric(CircRNACount_filtered_combined$Start[i])
  end_var <- as.numeric(CircRNACount_filtered_combined$End[i])

  # Find the corresponding gene in the GTF dataframe:
  matching_gene <- subset(gtf, chr == chr_var & feature == "gene" & start <= start_var & end >= end_var)
  
  if (nrow(matching_gene) == 0){
    matching_gene_start <- subset(gtf, chr == chr_var & feature == "gene" & start <= start_var)
    matching_gene_end <- subset(gtf, chr == chr_var & feature == "gene" & end >= end_var)
    matching_gene_start <- matching_gene_start[nrow(matching_gene_start),]
    matching_gene_end <- matching_gene_end[1,]
    if (matching_gene_start$end >= start_var & matching_gene_end$start <= end_var){
      gene_name <- paste0(matching_gene_start$gene_name, "; ", matching_gene_end$gene_name)
    }
    else if (matching_gene_start$end >= start_var & matching_gene_end$start >= end_var){
      gene_name <- matching_gene_start$gene_name
    }
    else if (matching_gene_start$end <= start_var & matching_gene_end$start <= end_var){
      gene_name <- matching_gene_end$gene_name
    }
    else {
      gene_name <- "not_annotated"
    }

    # Replace the original gene name with the updated one:
    CircRNACount_filtered_combined$Gene[i] <- gene_name
  }
  
  if (nrow(matching_gene) > 0){
    gene_name <- matching_gene$gene_name

    # Replace the original gene name with the updated one:
    CircRNACount_filtered_combined$Gene[i] <- gene_name
  }
}

for (i in seq_len(nrow(LinearCount_filtered_combined_0_as_fraction))) {
  
  # Extract chromosome, start, and end coordinates from the data frame:
  chr_var <- LinearCount_filtered_combined_0_as_fraction$Chr[i]
  start_var <- as.numeric(LinearCount_filtered_combined_0_as_fraction$Start[i])
  end_var <- as.numeric(LinearCount_filtered_combined_0_as_fraction$End[i])

  # Find the corresponding gene in the GTF dataframe:
  matching_gene <- subset(gtf, chr == chr_var & feature == "gene" & start <= start_var & end >= end_var)
  
  if (nrow(matching_gene) == 0){
    matching_gene_start <- subset(gtf, chr == chr_var & feature == "gene" & start <= start_var)
    matching_gene_end <- subset(gtf, chr == chr_var & feature == "gene" & end >= end_var)
    matching_gene_start <- matching_gene_start[nrow(matching_gene_start),]
    matching_gene_end <- matching_gene_end[1,]
    if (matching_gene_start$end >= start_var & matching_gene_end$start <= end_var){
      gene_name <- paste0(matching_gene_start$gene_name, "; ", matching_gene_end$gene_name)
    }
    else if (matching_gene_start$end >= start_var & matching_gene_end$start >= end_var){
      gene_name <- matching_gene_start$gene_name
    }
    else if (matching_gene_start$end <= start_var & matching_gene_end$start <= end_var){
      gene_name <- matching_gene_end$gene_name
    }
    else {
      gene_name <- "not_annotated"
    }

    # Replace the original gene name with the updated one:
    LinearCount_filtered_combined_0_as_fraction$Gene[i] <- gene_name
  }
  
  if (nrow(matching_gene) > 0){
    gene_name <- matching_gene$gene_name

    # Replace the original gene name with the updated one:
    LinearCount_filtered_combined_0_as_fraction$Gene[i] <- gene_name
  }
}

for (i in seq_len(nrow(test_summary_table_ordered))) {
  
  # Extract chromosome, strand, start, and end coordinates from the Coordinates data frame:
  chr_var <- test_summary_table_ordered$Chr[i]
  strand_var <- test_summary_table_ordered$Strand[i]
  start_var <- as.numeric(test_summary_table_ordered$Start[i])
  end_var <- as.numeric(test_summary_table_ordered$End[i])

  # Find the corresponding gene in the GTF dataframe:
  matching_gene <- subset(gtf, chr == chr_var & feature == "gene" & start <= start_var & end >= end_var)
  
  if (nrow(matching_gene) == 0){
    matching_gene_start <- subset(gtf, chr == chr_var & feature == "gene" & start <= start_var)
    matching_gene_end <- subset(gtf, chr == chr_var & feature == "gene" & end >= end_var)
    matching_gene_start <- matching_gene_start[nrow(matching_gene_start),]
    matching_gene_end <- matching_gene_end[1,]
    if (matching_gene_start$end >= start_var & matching_gene_end$start <= end_var){
      gene_name <- paste0(matching_gene_start$gene_name, "; ", matching_gene_end$gene_name)
    }
    else if (matching_gene_start$end >= start_var & matching_gene_end$start >= end_var){
      gene_name <- matching_gene_start$gene_name
    }
    else if (matching_gene_start$end <= start_var & matching_gene_end$start <= end_var){
      gene_name <- matching_gene_end$gene_name
    }
    else {
      gene_name <- "not_annotated"
    }

    # Replace the original gene name with the updated one:
    test_summary_table_ordered$Gene[i] <- gene_name
  }
  
  if (nrow(matching_gene) > 0){
    gene_name <- matching_gene$gene_name

    # Replace the original gene name with the updated one:
    test_summary_table_ordered$Gene[i] <- gene_name
  }
}
```

# Align and order all the dataframes with each other so we plot them correctly:
```{r}
test_p_summary_table <- test_summary_table_ordered[test_summary_table_ordered$sig_p < 0.05,]
#test_p_summary_table <- test$summary_table[test$summary_table$sig_p < 0.05,]

matching_rows <- rownames(test_p_summary_table)

matching_features <- test_p_summary_table[matching_rows, c("Chr", "Start", "End", "Gene")]

matching_rows <- match(paste(matching_features$Chr, matching_features$Start, matching_features$End, matching_features$Gene), paste(CircRNACount_filtered_combined$Chr, CircRNACount_filtered_combined$Start, CircRNACount_filtered_combined$End, CircRNACount_filtered_combined$Gene))

# Subset and reorder
subset_CircRNACount <- CircRNACount_filtered_combined[matching_rows, ]
subset_LinearCount <- LinearCount_filtered_combined_0_as_fraction[matching_rows, ]
subset_CircRNACount <- subset_CircRNACount[order(subset_CircRNACount[,1],subset_CircRNACount[,2],subset_CircRNACount[,3]),]
subset_LinearCount <- subset_LinearCount[order(subset_LinearCount[,1],subset_LinearCount[,2],subset_LinearCount[,3]),]
subset_test <- test_p_summary_table[order(test_p_summary_table[,1],test_p_summary_table[,2],test_p_summary_table[,3]),]
rownames(subset_CircRNACount) <- NULL
rownames(subset_LinearCount) <- NULL
rownames(subset_test) <- NULL
```

# Circ.ratioplot:
```{r}
## remember that these bar graphs will represent Circ and Linear counts, not necessarily the ratio mean from test$summary_table:
## the ratio mean from test$summary_table and subset_test EXCLUDES NA values
## we changed all linear NA values to 0.00001 in order to generate these graphs
for (j in rownames(subset_test)){
my.Circ.ratioplot(subset_CircRNACount[,-ncol(subset_CircRNACount)],   ## the last column of CircRNACount_filtered_combined is the "gene" column we just created 2 chunks ago
                 subset_LinearCount[,-ncol(subset_LinearCount)],
                 subset_CircRNACount[,c(1:3, ncol(subset_CircRNACount))],
                 plotrow = j,
                 size = 14,
                 groupindicator1 = c(rep(singletypes, times = sapply(singletypes, function(x) sum(grepl(paste0("^", x), singletypes_withreps))))),
                 lab_legend = "groupindicator1",
                 gene_column = 4
                 )
 }
```

# Circ.ratioplot plots the true Circ and Linear values as opposed to the ratio means values for each cell type. Here is a ggplot function if you actually want to plot the ratio means:
```{r}
for (i in 1:nrow(subset_test)) {
  row_data <- subset_test[i,]
  
  ratio_columns <- grepl("ratio_mean", names(row_data))
  row_df <- data.frame(Ratio_Type = names(row_data)[ratio_columns],
                       Ratio_Mean = unlist(row_data[ratio_columns]))
  
  ratio_names <- c(singletypes)
  
# Create a separate plot for each row
 p <- ggplot(row_df, aes(x = Ratio_Type, y = Ratio_Mean)) +
  geom_bar(stat = "identity", aes(fill = Ratio_Type), color = "black", size = 1) +
  scale_x_discrete(labels = setNames(ratio_names, row_df$Ratio_Type)) +  # Set x-axis labels
  scale_fill_discrete(labels = ratio_names[!is.nan(row_df$Ratio_Mean)]) +  # Set legend labels for non-NaN values
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  theme(axis.text.y = element_text(size = 8)) +
  theme(legend.title = element_blank()) +
  theme(text = element_text(size = 9)) +
  theme(legend.text = element_text(size = 8)) +
  theme(plot.title = element_text(size = 10)) +
  scale_y_continuous(labels = percent_format(scale = 100), expand = c(0, 0)) +  # Format y-axis as percentage
  ggtitle(paste("Annotation: ", row_data$Gene, "\nChr ", row_data$Chr, ":", row_data$Start, "-", row_data$End, sep = "")) +
  ylab("Percentage of Circular RNA vs All RNA") +
  theme(axis.title.x = element_blank())
  
  print(p)
}
```

# Save your work:
```{r}
write.csv(test, file = "CircRNA test results.csv")
write.csv(test$summary_table, file = "CircRNA test results summary.csv")
write.csv(circRNA_table, file = "circRNA table.csv")
write.csv(CircRNACount_filtered_combined, file = "CircRNACount_filtered_combined.csv")
write.csv(LinearCount_filtered_combined, file = "LinearCount_filtered_combined.csv")
write.csv(CircCoordinates_filtered_combined, file = "CircCoordinates_filtered_combined.csv")
```

# Let's continue filtering independently to find the "most interesting" circular RNAs:
```{r, message=FALSE}
# consolidate test table:
circ_test_consolidated <- test$summary_table[,c(1:4,10:54)]

# add a p-value cutoff:
# circ_test_consolidated <- circ_test_consolidated[circ_test_consolidated$sig_p < 0.05,]

# consolidate the table by gene name:
unique_circ_genes <- unique(circ_test_consolidated$Gene)
  
circ_gene_means <- data.frame()
for(i in unique_circ_genes){
  for (j in 5:ncol(circ_test_consolidated)){
   circ_gene_means[i,(j-4)] <- mean(circ_test_consolidated[circ_test_consolidated$Gene == i, j])
 }
}

# add column names:
names(circ_gene_means) <- c(singletypes)

circ_gene_count <- as.data.frame(table(circ_test_consolidated$Gene))
circ_gene_count <- circ_gene_count[rev(order(circ_gene_count$Freq)),]
names(circ_gene_count) <- c("Gene name", "Freq")

circ_gene_means[is.na(circ_gene_means)] <- 0
```

# Instead of consolidating by gene, consolidate by start coordinate:
```{r}
unique_circ_coordinates <- unique(as.character(circ_test_consolidated$Start))

circ_coordinates_means <- data.frame()
for(i in unique_circ_coordinates){
  for (j in 5:ncol(circ_test_consolidated)){
   circ_coordinates_means[i,(j-4)] <- mean(circ_test_consolidated[circ_test_consolidated$Start == i, j])
 }
}

# add column names:
names(circ_coordinates_means) <- c(singletypes)

circ_coordinates_count <- as.data.frame(table(circ_test_consolidated$Start))
circ_coordinates_count <- circ_coordinates_count[rev(order(circ_coordinates_count$Freq)),]
names(circ_coordinates_count) <- c("Start coordinate", "Freq")

circ_coordinates_means[is.na(circ_coordinates_means)] <- 0
```

# Heatmap plots:
```{r}
colors <- colorRampPalette(brewer.pal(6, "Oranges"))(250)
pheatmap(circ_gene_means,
         border_color = "black",
         cluster_rows = F,
         cluster_cols = F,
         display_numbers = F,
         fontsize = 9,
         fontsize_number = 16,
         angle_col = 90,
         col = colors,
         main = "Mean Proportion of Circular to Total RNAs for Each Cell Type")

circ_gene_means_RBPs <- circ_gene_means[rownames(circ_gene_means) %in% RBP_common_names$Public.Name,]
circ_gene_means_RBPs <- circ_gene_means_RBPs[order(circ_gene_means_RBPs$AWC, decreasing = TRUE),] # change cell type of interest here

colors <- colorRampPalette(c("white", "darkblue"))(250)
pheatmap(circ_gene_means_RBPs,
         border_color = "black",
         cluster_rows = F,
         cluster_cols = F,
         display_numbers = F,
         fontsize = 8,
         fontsize_number = 16,
         angle_col = 90,
         col = colors,
         main = "Mean Proportion of Circular to Total RNAs of RNA Binding Proteins for Each Cell Type")

colors <- colorRampPalette(brewer.pal(6, "PuRd"))(250)
pheatmap(circ_coordinates_means,
         border_color = "black",
         cluster_rows = F,
         cluster_cols = F,
         display_numbers = F,
         fontsize = 9,
         fontsize_number = 16,
         angle_col = 90,
         col = colors,
         main = "Mean Proportion of Circular to Total RNAs for Each Cell Type")
```

Copyright 2024 The Regents of the University of California

All Rights Reserved

Created by Zachery Wolfe

Department of Biochemistry

This file is part of Differential Expression in C. elegans. \
Differential Expression in C. elegans is a free pipeline: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version. \
Differential Expression in C. elegans is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details. \
You should have received a copy of the GNU General Public License along with Differential Expression in C. elegans. If not, see <https://www.gnu.org/licenses/>.
