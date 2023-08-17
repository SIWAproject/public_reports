# https://github.com/umerijaz/microbiomeSeq/blob/ffa009067f6f3fa229c9c453b82b027e07540474/R/perform_anova.R
library(dplyr)
library(ggplot2)
library(phyloseq)
library(tidyMicro)


diff_geneexp <-
  function(df_gene, GroupRef, variable="KitID", print=FALSE) {
    # IN: df with DeltaCq - out: DD y RQ for cec and ileum separated given the KitRef
    # required columns: KitID, SampleLocation
    targets = c("IL10", "IL1B", "MUC2")
    cec <- df_gene[df_gene$SampleLocation == "C",]
    il <- df_gene[df_gene$SampleLocation == "I",]
    
    if (GroupRef %in% df_gene[[variable]]){
      refs <-
        df_gene[df_gene[[variable]] == GroupRef, ] %>% dplyr::group_by(SampleLocation) %>% dplyr::summarise(
          DeltaCq_IL10 = mean(DeltaCq_IL10, na.rm=TRUE),
          DeltaCq_IL1B = mean(DeltaCq_IL1B, na.rm=TRUE),
          DeltaCq_MUC2 = mean(DeltaCq_MUC2, na.rm=TRUE)
        )
      if (print){print(refs)}
      
      for (t in targets) {
        cec[paste0("DD_", t)] = cec[, paste0("DeltaCq_", t)] - as.numeric(refs[refs$SampleLocation == "C", paste0("DeltaCq_", t)])
        il[paste0("DD_", t)] = il[, paste0("DeltaCq_", t)] - as.numeric(refs[refs$SampleLocation == "I", paste0("DeltaCq_", t)])
        cec[paste0("RQ_", t)] = 2 ^ -(cec[, paste0("DD_", t)])
        il[paste0("RQ_", t)] = 2 ^ -(il[, paste0("DD_", t)])
        cec[[paste0("neg_DD_", t)]] <- -cec[[paste0("DD_", t)]]
        il[[paste0("neg_DD_", t)]] <- -il[[paste0("DD_", t)]]
      }
      list("cec" = cec, "il" = il)
    } else {
      print("Reference Group not present in Variable column")
    }
    
  }


aggregate_taxa_siwa <-
  function(phylo_relative, taxa_level) {
    # se agregan todos los OTUs del Genus en la OTU table - 
    # coge Genus unkown como 1 solo genus
    OTUS <- as.data.frame(phyloseq::otu_table(phylo_relative))
    OTUS$otu <- rownames(OTUS)
    taxa <- as.data.frame(phyloseq::tax_table(phylo_relative)) 
    taxa$otu <- rownames(taxa)
    taxa_rank <- taxa[,c("otu", taxa_level)]
    OTUS_TAXA <- dplyr::left_join(OTUS, taxa_rank, by="otu")
    OTUS_TAXA <- OTUS_TAXA %>% dplyr::select(-otu)
    grouped <- aggregate(as.formula(paste0(". ~ ", taxa_level)), data = OTUS_TAXA, FUN = sum)
    #grouped <- aggregate(. ~ Family, data = OTUS_TAXA, FUN = sum)
    rownames(grouped) <- grouped[,taxa_level]
    grouped[,taxa_level] <- NULL
    #melted <- psmelt(otu_table(grouped,taxa_are_rows=TRUE))
    grouped
  }
  
aggregate_otus_by_taxa <-
  function(metadata, features, taxonomy, taxa_level) {
    rownames(features) <- features$OTU
    rownames(taxonomy) <- taxonomy$OTU
    rownames(metadata) <- metadata$SampleID
    features <- dplyr::select(features,-c(OTU))
    #taxonomy <- dplyr::select(taxonomy,-c(OTU, SciName))
    taxonomy <- dplyr::select(taxonomy,-c(OTU))
    taxa_matrix <- as.matrix(taxonomy)
    featuresM <- as(features, "matrix")
    class(featuresM) <- "numeric"
    ASVobj = otu_table(featuresM, taxa_are_rows = TRUE)
    TAXobj = tax_table(taxa_matrix)
    SAMobj = sample_data(metadata)
    phyloseq_ob = phyloseq(ASVobj, TAXobj, SAMobj)
    phyloseq_ob_agg <-
      microbiome::aggregate_taxa(phyloseq_ob, taxa_level)
    keep <- list(unique(taxonomy[taxa_level]))
    phyloseq_ob_agg <-
      phyloseq::prune_taxa(unique(taxonomy %>% pull(taxa_level)), phyloseq_ob_agg)
    df <- as.data.frame(otu_table(phyloseq_ob_agg))
    return(df)
  }


alpha_div_edit <- function(physeq,method){
  #==check for validity of selected methods
  ## MAKE SURE THE GROUPING VAR IS FACTOR
  method<- match.arg(method,c("richness", "fisher", "simpson", "shannon", "evenness","pd"), several.ok = TRUE)
  
  abund_table <- as.data.frame(as.matrix(otu_table(physeq)))
  df <- NULL
  if("richness"%in%method){
    tab <- otu_table(physeq)
    class(tab) <- "matrix" 
    #tab <- t(tab) - ya se hizo arriba
    R <- vegan::rarefy(tab,min(rowSums(tab)))
    #R<- vegan::rarefy(abund_table,min(rowSums(abund_table)))
    df_R<-data.frame(sample=names(R),value=R,measure=rep("Richness",length(R)))
    if(is.null(df)){
      df<-df_R}
    else {
      df<-rbind(df,df_R)}
  }
  if("fisher"%in%method){
    alpha <- vegan::fisher.alpha(abund_table, index = "shannon", base=2)
    df_alpha<-data.frame(sample=names(alpha),value=alpha,measure=rep("Fisher alpha",length(alpha)))
    if(is.null(df)){
      df<-df_alpha}
    else {
      df<-rbind(df,df_alpha)}
  }
  if("simpson"%in%method){
    simp <- vegan::diversity(abund_table, "simpson")
    df_simp<-data.frame(sample=names(simp),value=simp,measure=rep("Simpson",length(simp)))
    if(is.null(df)){
      df<-df_simp}
    else {
      df<-rbind(df,df_simp)}
  }
  if("shannon"%in%method){
    H<- vegan::diversity(abund_table, base=2)
    df_H<-data.frame(sample=names(H),value=H,measure=rep("Shannon",length(H)))
    if(is.null(df)){
      df<-df_H}
    else {
      df<-rbind(df,df_H)}
  }
  if("evenness"%in%method){
    H<-vegan::diversity(abund_table)
    S <- specnumber(abund_table)
    J <- H/log(S)
    df_J<-data.frame(sample=names(J),value=J,measure=rep("Pielou's evenness",length(J)))
    if(is.null(df)){
      df<-df_J}
    else {
      df<-rbind(df,df_J)}
  }
  if("pd"%in%method){
    otu_tree <- phyloseq::phy_tree(physeq)
    PD <- pd(abund_table, otu_tree ,include.root = TRUE)
    df_PD<-data.frame(sample=names(PD),value=PD,measure=rep("PD",length(PD)))
    if(is.null(df)){
      df<-df_PD}
    else {
      df<-rbind(df,df_PD)}
  }
  return(df)
}

#perform_anova_edit(df, meta_table,grouping_column,pValueCutoff)

perform_anova_edit <- function(df, meta_table,grouping_column,pValueCutoff){
  dt<-data.table::data.table(data.frame(df,.group.=meta_table[,grouping_column]))
  #specifying a p-value cutoff for the ggplot2 strips
  pval<-dt[, list(pvalue = sprintf("%.2g",
                                   tryCatch(summary(aov(value ~ .group.))[[1]][["Pr(>F)"]][1],error=function(e) NULL))),
           by=list(measure)]
  #Filter out pvals that we are not significant
  pval<-pval[!pval$pvalue=="",]
  pval<-pval[as.numeric(pval$pvalue)<=pValueCutoff,]
  
  #using sapply to generate significances for pval$pvalue using the cut function.
  pval$pvalue<-sapply(as.numeric(pval$pvalue),function(x){as.character(cut(x,breaks=c(-Inf, 0.001, 0.01, pValueCutoff, Inf),label=c("***", "**", "*", "")))})
  #Update df$measure to change the measure names if the grouping_column has more than three classes
  if(length(unique(as.character(meta_table[,grouping_column])))>2){
    df$measure<-as.character(df$measure)
    if(dim(pval)[1]>0){
      for(i in seq(1:dim(pval)[1])){
        df[df$measure==as.character(pval[i,measure]),"measure"]=paste(as.character(pval[i,measure]),as.character(pval[i,pvalue]))
      }
    }
    df$measure<-as.factor(df$measure)
  }
  #Get all possible pairwise combination of values in the grouping_column
  s<-combn(unique(as.character(df[,grouping_column])),2)
  
  #df_pw will store the pair-wise p-values
  df_pw<-NULL
  pValueCutoff=0.05
  for(k in unique(as.character(df$measure))){
    #We need to calculate the coordinate to draw pair-wise significance lines
    #for this we calculate bas as the maximum value
    bas<-max(df[(df$measure==k),"value"])
    #Calculate increments as 10% of the maximum values
    inc<-0.1*bas
    #Give an initial increment
    bas<-bas+inc
    for(l in 1:dim(s)[2]){
      #Do a pair-wise anova
      tmp<-c(k,s[1,l],s[2,l],bas,paste(sprintf("%.2g",
                                               tryCatch(summary(aov(as.formula(paste("value ~",grouping_column)),data=df[(df$measure==k) & (df[,grouping_column]==s[1,l] | df[,grouping_column]==s[2,l]),] ))[[1]][["Pr(>F)"]][1],error=function(e) NULL)),"",sep=""))
      #Ignore if anova fails
      if(!is.na(as.numeric(tmp[length(tmp)]))){
        #Only retain those pairs where the p-values are significant
        if(as.numeric(tmp[length(tmp)])<pValueCutoff){
          if(is.null(df_pw)){df_pw<-tmp}else{df_pw<-rbind(df_pw,tmp)}
          #Generate the next position
          bas<-bas+inc
        }
      }
    }
  }
  if(!is.null(df_pw)){
    if (!is.null(dim(df_pw))){
      df_pw<-data.frame(row.names=NULL,df_pw)
      names(df_pw)<-c("measure","from","to","y","p")
    } else {
      df_pw<-data.frame(row.names=NULL, t(df_pw))
      names(df_pw)<-c("measure","from","to","y","p")
    }
  }
  out <- list("df_pw"=df_pw, "df"=df)
  return(out)
}


plot_anova_diversity_edit <- function(physeq, method, grouping_column,
                                      pValueCutoff=0.05, outfile="anova_diversity.csv")
{
  #enforce orientation
  if(taxa_are_rows(physeq)){
    physeq <- t(physeq)
  }
  abund_table <- otu_table(physeq)
  meta_table <- sample_data(physeq)
  
  #get diversity measure using selected methods
  div.df <- alpha_div_edit(physeq,method)
  
  #=add grouping information to alpha diversity measures
  df<-data.frame(div.df,(meta_table[,grouping_column])[as.character(div.df$sample),])
  #perform anova of diversity measure between groups
  anova_res <- perform_anova_edit(df, meta_table,grouping_column,pValueCutoff)
  df_pw <- anova_res$df_pw #get pairwise p-values
  #write.csv(df_pw, file = outfile)
  
  #Draw the boxplots
  p<-ggplot(aes_string(x=grouping_column,y="value",color=grouping_column),data=df)
  p<-p+geom_boxplot()+geom_jitter(position = position_jitter(height = 0, width=0))
  p<-p+theme_bw()
  p<-p+theme(axis.text.x = element_text(angle = 90, hjust = 1))
  p<-p+facet_wrap(~measure,scales="free_y",nrow=1)+ylab("Observed Values")+xlab("Samples")
  p<-p+theme(strip.background = element_rect(fill = "white"))+xlab("Groups")
  
  #This loop will generate the lines and signficances
  if(!is.null(df_pw)){ #this only happens when we have significant pairwise anova results
    for(i in 1:dim(df_pw)[1]){
      p<-p+geom_path(inherit.aes=F,aes(x,y),data = data.frame(x = c(which(levels(df[,grouping_column])==as.character(df_pw[i,"from"])),which(levels(df[,grouping_column])==as.character(df_pw[i,"to"]))), y = c(as.numeric(as.character(df_pw[i,"y"])),as.numeric(as.character(df_pw[i,"y"]))), measure=c(as.character(df_pw[i,"measure"]),as.character(df_pw[i,"measure"]))), color="black",lineend = "butt",arrow = arrow(angle = 90, ends = "both", length = unit(0.1, "inches")))
      p<-p+geom_text(inherit.aes=F,aes(x=x,y=y,label=label),data=data.frame(x=(which(levels(df[,grouping_column])==as.character(df_pw[i,"from"]))+which(levels(df[,grouping_column])==as.character(df_pw[i,"to"])))/2,y=as.numeric(as.character(df_pw[i,"y"])),measure=as.character(df_pw[i,"measure"]),label=as.character(cut(as.numeric(as.character(df_pw[i,"p"])),breaks=c(-Inf, 0.001, 0.01, 0.05, Inf),label=c("***", "**", "*", "")))))
    }
  }
  return(p)
}


plot_anova_diversity_edit_treatment_order <- function(physeq, method, grouping_column,pValueCutoff=0.05, outfile="anova_diversity.csv", order_levels)
{
  #enforce orientation
  if(taxa_are_rows(physeq)){
    physeq <- t(physeq)
  }
  abund_table <- otu_table(physeq)
  meta_table <- sample_data(physeq)
  
  #get diversity measure using selected methods
  div.df <- alpha_div_edit(physeq,method)
  
  #=add grouping information to alpha diversity measures
  df<-data.frame(div.df,(meta_table[,grouping_column])[as.character(div.df$sample),])
  
  #perform anova of diversity measure between groups
  anova_res <- perform_anova_edit(df, meta_table,grouping_column,pValueCutoff)
  df_pw <- anova_res$df_pw #get pairwise p-values
  #write.csv(df_pw, file = outfile)
  
  #order of treatments to use in graph
  df <- df %>%
    mutate(Treatment=factor(Treatment,levels=order_levels))
  
  #Draw the boxplots
  p<-ggplot(aes_string(x=grouping_column,y="value",color=grouping_column),data=df)
  p<-p+geom_boxplot()+geom_jitter(position = position_jitter(height = 0, width=0))
  p<-p+theme_bw()
  p<-p+theme(axis.text.x = element_text(angle = 90, hjust = 1))
  p<-p+facet_wrap(~measure,scales="free_y",nrow=1)+ylab("Observed Values")+xlab("Samples")
  p<-p+theme(strip.background = element_rect(fill = "white"))+xlab("Groups")
  
  #This loop will generate the lines and signficances
  if(!is.null(df_pw)){ #this only happens when we have significant pairwise anova results
    for(i in 1:dim(df_pw)[1]){
      p<-p+geom_path(inherit.aes=F,aes(x,y),data = data.frame(x = c(which(levels(df[,grouping_column])==as.character(df_pw[i,"from"])),which(levels(df[,grouping_column])==as.character(df_pw[i,"to"]))), y = c(as.numeric(as.character(df_pw[i,"y"])),as.numeric(as.character(df_pw[i,"y"]))), measure=c(as.character(df_pw[i,"measure"]),as.character(df_pw[i,"measure"]))), color="black",lineend = "butt",arrow = arrow(angle = 90, ends = "both", length = unit(0.1, "inches")))
      p<-p+geom_text(inherit.aes=F,aes(x=x,y=y,label=label),data=data.frame(x=(which(levels(df[,grouping_column])==as.character(df_pw[i,"from"]))+which(levels(df[,grouping_column])==as.character(df_pw[i,"to"])))/2,y=as.numeric(as.character(df_pw[i,"y"])),measure=as.character(df_pw[i,"measure"]),label=as.character(cut(as.numeric(as.character(df_pw[i,"p"])),breaks=c(-Inf, 0.001, 0.01, 0.05, Inf),label=c("***", "**", "*", "")))))
    }
  }
  return(p)
}



filter_deseq <- function(deseq_list, phyloseqobj, alpha=0.1){
  ## Le entran los resultados de deseq en la lista, y filtra los que pasen el threshold
  results <- deseq_list[["results"]]
  print(paste0("IN = ", length(results)))
  groupsNames <- deseq_list[["groupnames"]]
  groups <- deseq_list[["groupscontrol"]] 
  dfs_alfaST <- lapply(results, function(x) {x[which(x$padj < alpha), ]})
  filter_outST <- sapply(dfs_alfaST, function(x) dim(x)[1])>0
  print(paste0("DIDN'T PASS THE FILTER = ",  matrixStats::count(filter_outST == FALSE)))
  groups<-groups[filter_outST]
  groupsNames<-groupsNames[filter_outST]
  dfs_alfaST <- dfs_alfaST[filter_outST]
  dfs_alfaST <- lapply(dfs_alfaST, function(x) {cbind(as(x, "data.frame"), as(tax_table(phyloseqobj)[rownames(x), ], "matrix"))}) 
  l <- list()
  l[["groups"]] <- groups
  l[["df"]] <- dfs_alfaST
  l
  
}


#simple function to plot only wanted comparison for the report
plot_deseq_report <-
  function(dfs_filtered, comparisonINDEX, Title, color_location, size=8) {
    theme_set(theme_bw())
    df <- dfs_filtered[[comparisonINDEX]]
    df$taxa <- df$SciName
    df$taxa <- as.character(df$taxa)
    x = sort(tapply(df$log2FoldChange, df$taxa, function(x)
      max(x)), TRUE)
    df$taxa = factor(as.character(df$taxa), levels = names(x))
    df$abs <- abs(df$log2FoldChange)
    df <- df[order(df$abs, decreasing = TRUE), ]
    if (dim(df)[[1]] > 10) {
      df <- df[1:10, ]
    }
    plot = ggplot(df, aes(x = log2FoldChange, y = taxa)) + geom_point(size = size, color = color_location) + 
      labs(fill = "taxa", color = "taxa")  +
      theme(legend.position = "none") + ggtitle(Title)
    return(plot)
  }


plot_deseq_report_for_genera <-
  function(dfs_filtered, comparisonINDEX, Title, color_location, size=8) {
    theme_set(theme_bw())
    df <- dfs_filtered[[comparisonINDEX]]
    df$taxa <- df$Genus
    df$taxa <- as.character(df$taxa)
    x = sort(tapply(df$log2FoldChange, df$taxa, function(x)
      max(x)), TRUE)
    df$taxa = factor(as.character(df$taxa), levels = names(x))
    df$abs <- abs(df$log2FoldChange)
    df <- df[order(df$abs, decreasing = TRUE), ]
    if (dim(df)[[1]] > 10) {
      df <- df[1:10, ]
    }
    plot = ggplot(df, aes(x = log2FoldChange, y = taxa)) + geom_point(size = size, color = color_location) + 
      labs(fill = "taxa", color = "taxa")  +
      theme(legend.position = "none") + ggtitle(Title)
    return(plot)
  }


plot_deseq_report_with_color <-
  # function to plot only wanted comparison for the report colored by probiotic score
  function(dfs_filtered, comparisonINDEX, Title, species_taxonomy_info, genera_taxonomy_info,  size=8) {
    theme_set(theme_bw())
    df <- dfs_filtered[[comparisonINDEX]]
    df$taxa <- df$SciName
    df$taxa <- as.character(df$taxa)
    x = sort(tapply(df$log2FoldChange, df$taxa, function(x)
      max(x)), TRUE)
    df$taxa = factor(as.character(df$taxa), levels = names(x))
    df$abs <- abs(df$log2FoldChange)
    df <- df[order(df$abs, decreasing = TRUE), ]
    if (dim(df)[[1]] > 10) {
      df <- df[1:10, ]
    }
    
    #Add to df column of probiotic score
    df_sp_score=species_taxonomy_info[c("Species", "Probiotic_potential")]
    df_sp_score$Species <- gsub(" ", "_", df_sp_score$Species)
    df_gen_score=genera_taxonomy_info[c("Genus", "Probiotic_potential")]
    df_gen_score$Genus <- gsub(" ", "_", df_gen_score$Genus)
    
    scorelist <- list()
    for(i in 1:nrow(df)) {
      row <- df[i,]
      sp <- row$Species
      gen <- row$Genus
      if (sp %in% df_sp_score$Species) {
        score = df_sp_score$Probiotic_potential[df_sp_score$Species==sp]
        
      } else if (gen %in% df_gen_score$Genus) {
        score =df_gen_score$Probiotic_potential[df_gen_score$Genus==gen]
        
      } else {score = "none"}
      scorelist <- append(scorelist, score)
    }
    
    df$ProbioticScore <- scorelist
    df$ ProbioticScore <- as.character(df$ProbioticScore)
    
    plot = ggplot(df, aes(x = log2FoldChange, y = taxa, color=ProbioticScore)) + geom_point(size = size) + 
      scale_color_manual(values = c("-6" = "#CB756E","-5" = "#CB766E", "-4" = "#CB766E","-3" = "#C87A6F","-2" = "#C67F71","-1" = "#c48473","0" = "#be9379","1" = "#b4ab81","2" = "#ADB987","3" = "#A9C38B","4" = "#A5CD8F","5" = "#A3D291"))+
      labs(fill = "taxa", color = "taxa")  +
      theme(legend.position = "none") + ggtitle(Title)+
      geom_vline(xintercept=c(0,0), linetype="dotted")+
      ylab("")
    return(plot)
  }


plot_deseq_report_with_color_for_genera <-
  # function to plot only wanted comparison for the report colored by probiotic score
  function(dfs_filtered, comparisonINDEX, Title, species_taxonomy_info, genera_taxonomy_info,  size=8) {
    theme_set(theme_bw())
    df <- dfs_filtered[[comparisonINDEX]]
    df$taxa <- df$Genus
    df$taxa <- as.character(df$taxa)
    x = sort(tapply(df$log2FoldChange, df$taxa, function(x)
      max(x)), TRUE)
    df$taxa = factor(as.character(df$taxa), levels = names(x))
    df$abs <- abs(df$log2FoldChange)
    df <- df[order(df$abs, decreasing = TRUE), ]
    if (dim(df)[[1]] > 10) {
      df <- df[1:10, ]
    }
    
    #Add to df column of probiotic score
    df_sp_score=species_taxonomy_info[c("Species", "Probiotic_potential")]
    df_sp_score$Species <- gsub(" ", "_", df_sp_score$Species)
    df_gen_score=genera_taxonomy_info[c("Genus", "Probiotic_potential")]
    df_gen_score$Genus <- gsub(" ", "_", df_gen_score$Genus)
    
    scorelist <- list()
    for(i in 1:nrow(df)) {
      row <- df[i,]
      sp <- row$Species
      gen <- row$Genus
      if (sp %in% df_sp_score$Species) {
        score = df_sp_score$Probiotic_potential[df_sp_score$Species==sp]
        
      } else if (gen %in% df_gen_score$Genus) {
        score =df_gen_score$Probiotic_potential[df_gen_score$Genus==gen]
        
      } else {score = "none"}
      scorelist <- append(scorelist, score)
    }
    
    df$ProbioticScore <- scorelist
    df$ ProbioticScore <- as.character(df$ProbioticScore)
    
    plot = ggplot(df, aes(x = log2FoldChange, y = taxa, color=ProbioticScore)) + geom_point(size = size) + 
      scale_color_manual(values = c("-6" = "#CB756E","-5" = "#CB766E", "-4" = "#CB766E","-3" = "#C87A6F","-2" = "#C67F71","-1" = "#c48473","0" = "#be9379","1" = "#b4ab81","2" = "#ADB987","3" = "#A9C38B","4" = "#A5CD8F","5" = "#A3D291"))+
      labs(fill = "taxa", color = "taxa")  +
      theme(legend.position = "none") + ggtitle(Title)+
      geom_vline(xintercept=c(0,0), linetype="dotted")+
      ylab("")
    return(plot)
  }


plot_deseq_report_with_color_easyread_for_genera <-
  # INPUT= LIST OF DF AND GROUPS 
  # function to plot only wanted comparison for the report colored by probiotic score
  function(list_dfs_filtered, comparisonINDEX, Title, species_taxonomy_info, genera_taxonomy_info,  size=8, top=10) {
    theme_set(theme_bw())
    df <- list_dfs_filtered$df[[comparisonINDEX]]
    df$taxa <- df$Genus
    df$taxa <- as.character(df$taxa)
    x = sort(tapply(df$log2FoldChange, df$taxa, function(x)
      max(x)), TRUE)
    df$taxa = factor(as.character(df$taxa), levels = names(x))
    df$abs <- abs(df$log2FoldChange)
    df <- df[order(df$abs, decreasing = TRUE), ]
    if (dim(df)[[1]] > top) {
      df <- df[1:top, ]
    }
    
    #Add to df column of probiotic score
    df_sp_score=species_taxonomy_info[c("Species", "Probiotic_potential")]
    df_sp_score$Species <- gsub(" ", "_", df_sp_score$Species)
    df_gen_score=genera_taxonomy_info[c("Genus", "Probiotic_potential")]
    df_gen_score$Genus <- gsub(" ", "_", df_gen_score$Genus)
    
    scorelist <- list()
    for(i in 1:nrow(df)) {
      row <- df[i,]
      sp <- row$Species
      gen <- row$Genus
      if (isTRUE(sp) && sp %in% df_sp_score$Species) {
        score = df_sp_score$Probiotic_potential[df_sp_score$Species==sp]
        
      } else if (gen %in% df_gen_score$Genus) {
        score =df_gen_score$Probiotic_potential[df_gen_score$Genus==gen]
        
      } else {score = "none"}
      scorelist <- append(scorelist, score)
    }
    
    df$ProbioticScore <- scorelist
    df$ProbioticScore <- as.character(df$ProbioticScore)
    
    plot = ggplot(df, aes(x = log2FoldChange, y = taxa, color=ProbioticScore)) + geom_point(size = size) + 
      scale_color_manual(values = c("-6" = "#CB756E","-5" = "#CB766E", "-4" = "#CB766E","-3" = "#C87A6F","-2" = "#C67F71","-1" = "#c48473","0" = "#be9379","1" = "#b4ab81","2" = "#ADB987","3" = "#A9C38B","4" = "#A5CD8F","5" = "#A3D291"), na.value="dark grey") +
      labs(fill = "taxa", color = "taxa", subtitle = c(paste("Higher in",list_dfs_filtered$groups[[comparisonINDEX]][2],sep=" "), paste("Higher in",list_dfs_filtered$groups[[comparisonINDEX]][1], sep=" ")))  +
      theme(legend.position = "none", plot.subtitle = element_text(size = 20, hjust = c(0,1), vjust = 1)) + 
      geom_vline(xintercept=0,linetype="dashed",color="dark grey") +
      ggtitle(Title) +
      ylab("")
    return(plot)
  }



plot_deseq <- function(dfs_filtered, groupsNames, taxonomy, save=FALSE, specific=0){
  theme_set(theme_bw())
  if (specific != 0){
    print("Just one plot")
    df <- dfs_filtered[[specific]]
    df$taxa <-lapply(rownames(df), function(x) {taxonomy[taxonomy$OTU == x, "SciName"]})
    df$taxa <- as.character(df$taxa)
    x = sort(tapply(df$log2FoldChange, df$taxa, function(x) max(x)), TRUE)
    df$taxa = factor(as.character(df$taxa), levels=names(x))
    df$abs <- abs(df$log2FoldChange)
    df <- df[ order( df$abs, decreasing = TRUE ),]
    if (dim(df)[[1]] > 10 ){ df <- df[1:10,] }
    plot = ggplot(df, aes(x=log2FoldChange, y=taxa, color=Species)) + geom_point(size=1)+ labs(fill="Species",color="Species") + 
      theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5)) + ggtitle(paste(unlist(groupsNames[[specific]]),collapse=' vs ' ))
    
    if (save){
        print(paste0(groupsNames[[specific]],".png"))
        png(paste0(groupsNames[[specific]], ".png"),res=130, width = 1000, height = 750)
        dev.off()
    }
    return(plot)
  }
  else{
    plotsST <- list()
    for (i in 1:length(dfs_filtered)){
      df <- dfs_filtered[[i]]
      df$taxa <-lapply(rownames(df), function(x) {taxonomy[taxonomy$OTU == x, "SciName"]})
      df$taxa <- as.character(df$taxa)
      x = sort(tapply(df$log2FoldChange, df$taxa, function(x) max(x)), TRUE)
      df$taxa = factor(as.character(df$taxa), levels=names(x))
      df$abs <- abs(df$log2FoldChange)
      df <- df[ order( df$abs, decreasing = TRUE ),]
      if (dim(df)[[1]] > 10 ){ df <- df[1:10,] }
      plot <- ggplot(df, aes(x=log2FoldChange, y=taxa, color=Species)) + geom_point(size=1)+ labs(fill="Species",color="Species") + 
        theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5)) + ggtitle(paste(unlist(groupsNames[[i]]),collapse=' vs ' ))
      plotsST[[i]] <- plot
    }
    if (save){
      for (i in 1:length(plotsST)){
        print(paste0(groupsNames[[i]],".png"))
        png(paste0(groupsNames[[i]],".png"),res=130, width = 1000, height = 750)
        dev.off()
      }
    }
  }
  return(plotsST)
}




plot_deseq_ancombc <-
  function(deseqF, ancomF, taxonomyEdit, numotu = 0) {
    "Si se pone numotu=0, se muestran todos"
    deseqF$log2FoldChangeScaled <- scale(deseqF$log2FoldChange)[, 1]
    df_fig =  deseqF %>%  transmute(otu, log2FoldChangeScaled) %>%
      filter(log2FoldChangeScaled != 0) %>% arrange(desc(log2FoldChangeScaled)) %>%
      mutate(group = ifelse(log2FoldChangeScaled > 0, "pos", "neg"))
    df_fig$label <-
      taxonomyEdit[match(df_fig$otu, taxonomyEdit$OTU), "SciName"]
    df_fig$label <-
      lapply(seq_along(1:nrow(df_fig)), function(i)
        if (sum(df_fig$label[[i]] == df_fig$label) > 1)
        {
          paste0(df_fig$label[[i]], "-", substring(df_fig$otu[[i]], 1, 10))
        } else {
          df_fig$label[[i]]
        })
    df_fig$label = factor(df_fig$label, levels = df_fig$label)
    df_fig = cbind(df_fig, type = "deseq")
    df_fig <- df_fig %>% dplyr::rename(value = log2FoldChangeScaled)
    ### parse ancombc
    ancomF$betaScaled <- scale(ancomF$beta)[, 1]
    df_fig2 =  ancomF %>%  transmute(otus, betaScaled) %>%
      filter(betaScaled != 0) %>% arrange(desc(betaScaled)) %>%
      mutate(group = ifelse(betaScaled > 0, "g1", "g2"))
    df_fig2$label <-
      taxonomyEdit[match(df_fig2$otus, taxonomyEdit$OTU), "SciName"]
    df_fig2$label <-
      lapply(seq_along(1:nrow(df_fig2)), function(i)
        if (sum(df_fig2$label[[i]] == df_fig2$label) > 1)
        {
          paste0(df_fig2$label[[i]], "-", substring(df_fig2$otu[[i]], 1, 10))
        } else {
          df_fig2$label[[i]]
        })
    df_fig2$label = factor(df_fig2$label, levels = df_fig2$label)
    df_fig2 = cbind(df_fig2, type = "ancombc")
    df_fig2 <-
      df_fig2 %>% dplyr::rename(value = betaScaled, otu = otus)
    ## merge
    df_full <- rbind(df_fig, df_fig2)
    df_full_both <- df_full %>%  group_by(otu) %>%  filter(n() > 1)
    if (dim(df_full_both)[[1]] == 0) {
      print("No intersection")
      return(NULL)
    }
    df_full_both =  df_full_both %>%  arrange(desc(otu))
    print(df_full_both)
    df_full_both$label <-
      taxonomyEdit[match(df_full_both$otu, taxonomyEdit$OTU), "SciName"]
    df_full_both$lastLevel <-
      taxonomyEdit[match(df_full_both$otu, taxonomyEdit$OTU), "granular_taxa"]
    df_full_both$finallabel <-
      unlist(lapply(seq_along(1:nrow(df_full_both)), function(i)
        if (df_full_both$label[[i]] == "UNKOWN-UNKOWN")
        {
          paste0(df_full_both$lastLevel[[i]], "-U-U")
        } else {
          df_full_both$label[[i]]
        }))
    
    df_full_both$finallabel <-
      lapply(seq_along(1:nrow(df_full_both)), function(i)
        if (sum(df_full_both$finallabel[[i]] == df_full_both$finallabel) > 2)
        {
          paste0(df_full_both$finallabel[[i]],
                 "-",
                 substring(df_full_both$otu[[i]], 1, 10))
        } else {
          df_full_both$finallabel[[i]]
        })
    df_full_both$finallabel = factor(df_full_both$finallabel, levels = unique(df_full_both$finallabel))
    if (numotu != 0) {
      if (numotu * 2 < dim(df_full_both)[[1]]) {
        print(paste(
          "Se muestran únicamente ",
          numotu,
          " OTUS de ",
          dim(df_full_both)[[1]] / 2
        )) ## otu=0  muestra todos los disponibles
        numfilas = 2 * numotu
        df_full_both <- df_full_both[1:numfilas,]
      } else {
        print(paste(
          "Numero de otus seleccionado es mayor que el total presentes -->",
          dim(df_full_both)[[1]]
        ))
      }
      
    }
    plot = ggplot(df_full_both) + geom_bar(
      aes(x = value, y = finallabel, fill = type),
      width = 0.7,
      colour = "black",
      stat = "identity",
      position = position_dodge(width = 0.4)
    )  + scale_fill_manual(values = c(deseq = "#343aeb", ancombc = "#BA4FC8")) +
      labs(y = NULL, x = "Normalized difference") + theme_bw() + guides(fill = guide_legend(title =
                                                                                              "Analysis")) + geom_vline(
                                                                                                xintercept = 0,
                                                                                                linetype = "solid",
                                                                                                color = " dark gray",
                                                                                                size = 1
                                                                                              ) + theme(legend.position = "bottom")
    return(plot)
    
  }


cvi_palettes = function(name, n, all_palettes = cvi_colours, type = c("discrete", "continuous")) {
  palette = all_palettes[[name]]
  if (missing(n)) {
    n = length(palette)
  }
  type = match.arg(type)
  out = switch(type,continuous = grDevices::colorRampPalette(palette)(n),discrete = palette[1:n]
  )
  structure(out, name = name, class = "palette")
}

