##################################################
## R script for MicrobiomeAnalyst
## Description: report generation using Sweave
## Note: most analyses were already performed, only need to embedding
## the results to the right place without rerun the whole analysis
## through Sweave. Only some auxilliary info (i.e. time, version etc need to
## run R throught Sweave
##
## Author: Jeff Xia, jeff.xia@mcgill.ca
## McGill University, Canada
##
## License: GNU GPL (>= 2)
###################################################

#'Function to create PDF report
#'@description This function creates a PDF report.
#'@param mbSetObj Input the name of the mbSetObj.
#'@param usrName Input the preferred user name for the Analysis Report.
#'@author Jeff Xia \email{jeff.xia@mcgill.ca}
#'McGill University, Canada
#'License: GNU GPL (>= 2)
#'@export
PreparePDFReport <- function(mbSetObj, usrName){
  
  mbSetObj <- .get.mbSetObj(mbSetObj);

  # create the Rnw file
  file.create("Analysis_Report.Rnw");
  # open for write
  rnwFile <<- file("Analysis_Report.Rnw", "w")

  # create a global counter to label figures
  fig.count <<- 0;
  table.count <<- 0;

  if(mbSetObj$module.type == "mdp" ){
    CreateMDPRnwReport(mbSetObj, usrName);
  }else if(mbSetObj$module.type == "sdp"){
    CreateSDPRnwReport(mbSetObj, usrName);
  }else if(mbSetObj$module.type == "ppd"){
    CreatePPDRnwReport(mbSetObj, usrName);
  }else{
    CreateTaxaEnrichRnwReport(mbSetObj, usrName);
  }
    
  # close opened files
  close(rnwFile);

  # all call from bash external to get around Sweave Fatigue 
  if(!.on.public.web){
    Sweave("Analysis_Report.Rnw", encoding="utf8");
    res <- try(tools::texi2dvi("Analysis_Report.tex", pdf = TRUE, quiet=TRUE));
  }

  return(.set.mbSetObj(mbSetObj));
  
}

# this is for PDF report generation from bash
SaveCurrentSession <- function(){
  file.copy("../../libs/Sweave.sty", ".")
  save.image("SweaveImage.RData");
}

################################
## MicrobiomeAnalyst statistics
###############################

#'Function to create PDF report for MDP module
#'@description This function creates a PDF report
#'for the MDP module - writes .Rnw file template.
#'@param mbSetObj Input the name of the mbSetObj.
#'@param usrName Input the preferred user name for the Analysis Report.
#'@author Jeff Xia \email{jeff.xia@mcgill.ca}
#'McGill University, Canada
#'License: GNU GPL (>= 2)
#'@export
CreateMDPRnwReport<-function(mbSetObj, usrName){
  
  mbSetObj <- .get.mbSetObj(mbSetObj);
  
  CreateHeader(usrName);
  CreateIntr();
  CreateIOdoc();
  CreateNORMdoc();
  Init16SAnalMode();
    
  if(!is.null(mbSetObj[["analSet"]]) & (length(mbSetObj$analSet)>0)){
    CreateVISEXPLRdoc(mbSetObj);
    CreateRAREFCTIONCURVEdoc(mbSetObj);
    CreatePHYLOGENETICTREEdoc(mbSetObj);
    CreateHEATTREEdoc(mbSetObj);
    CreateALPHDIVdoc(mbSetObj);
    CreateBETADIVdoc(mbSetObj);
    CreateHCdoc(mbSetObj);
    CreateCOREdoc(mbSetObj);
    CreateCORRdoc(mbSetObj);
    CreateFUNCPREDdoc(mbSetObj);
    CreateUNIVARdoc(mbSetObj);
    CreateRNASEQdoc(mbSetObj);
    CreateMETAGENOSEQdoc();
    CreateLEFSEdoc(mbSetObj);
    CreateRFdoc(mbSetObj);
  } else {
    CreateAnalNullMsg();
    }
  CreateFooter(mbSetObj);
}

#'Function to create PDF report for SDP module
#'@description This function creates a PDF report
#'for the SDP module - writes .Rnw file template.
#'@param mbSetObj Input the name of the mbSetObj.
#'@param usrName Input the preferred user name for the Analysis Report.
#'@author Jeff Xia \email{jeff.xia@mcgill.ca}
#'McGill University, Canada
#'License: GNU GPL (>= 2)
#'@export
CreateSDPRnwReport<-function(mbSetObj, usrName){
  
  mbSetObj <- .get.mbSetObj(mbSetObj);
    
  CreateHeader(usrName);
  CreateIntr();
  CreateIOdoc(mbSetObj);
  CreateNORMdoc(mbSetObj);
  InitShotAnalMode();
    
  if(!is.null(mbSetObj[["analSet"]]) & (length(mbSetObj$analSet)>0)){
    CreateFUNCPROFdoc(mbSetObj);
    CreateHCdoc(mbSetObj);
    CreateCORRdoc(mbSetObj);
    CreatePCAdoc(mbSetObj);
    CreateUNIVARdoc(mbSetObj);
    CreateRNASEQdoc(mbSetObj);
    CreateMETAGENOSEQdoc(mbSetObj);
    CreateLEFSEdoc(mbSetObj);
    CreateRFdoc(mbSetObj);
  } else {
    CreateAnalNullMsg();
  }
  CreateFooter(mbSetObj);
}

# create header
CreateHeader <- function(usrName){
    
  header <- c("\\documentclass[a4paper]{article}",
            "\\usepackage[margin=1.0in]{geometry}",
            "\\usepackage{longtable}",
            "\\SweaveOpts{eps=FALSE,pdf=TRUE}",
            "\\title{Microbiome Data Analysis with MicrobiomeAnalyst}",
            paste("\\author{ User ID: ", usrName, " }", sep=""),
            "\\begin{document}",
            "\\parskip=.3cm",
            "\\maketitle");
  cat(header, file=rnwFile, sep="\n", append=TRUE);
}

CreateIntr <- function(){
  descr <- c("\\section{Data Processing and Normalization}\n");
  cat(descr, file=rnwFile, append=TRUE);
}

# read and process the raw data
CreateIOdoc <- function(mbSetObj){
  
  mbSetObj <- .get.mbSetObj(mbSetObj);
  
  if(mbSetObj$module.type == "sdp"){
    descr <- c("\\subsection{Reading and Processing the Raw Data}\n",
              "MicrobiomeAnalyst accepts gene abundance profile in two formats generated from shotgun metagenomics or",
              "metatranscriptomics studies including plain text table and biom format.",
              "Users need to upload data in either of these two formats and specify the gene ID",
              "when uploading their data in order for MicrobiomeAnalyst to  process the information correctly.",
              "Currently, genes annotated as KEGG Orthology IDs (KO), Enzyme Commission Numbers (EC),", 
              "and Cluster of Orthologous Groups (COG) are supported in MicrobiomeAnalyst.\n",
              "Also, uploading metadata file is required as plain text (.txt or.csv) with both these formats.\n"
              );
  } else {
    descr <- c("\\subsection{Reading and Processing the Raw Data}\n",
              "MicrobiomeAnalyst accepts count data in a variety of formats generated in microbiome studies",
              "including plain text table, biom format as well as output from mothur pipeline.",
              "User need to upload their data in one of three available formats and specify the taxonomic labels",
              "when uploading their data in order for MicrobiomeAnalyst to  process the taxonomic information correctly.",
              "The hierarchical information for taxa can either be present within abundance data or uploaded as a", 
              "separate taxonomy table file (.txt or .csv format).\n",
              "Also, uploading metadata file is required as plain text (.txt or.csv) with all three formats.\n"
              );
  }
        
  cat(descr, file=rnwFile, append=TRUE);

  if(mbSetObj$dataSet$data.type=="text"){
    descr<-c("\\subsubsection{Reading abundance count data table}\n",
            "The abundance count data should be uploaded in tab-delimited text (.txt) or comma separated values (.csv)",
            "format. Samples are represented in columns, while rows contains the information about the features.",
            "Metadata file contains additional information about samples such as experimental factors or sample", 
            "grouping.\n");

  }else if(mbSetObj$dataSet$data.type=="biom"){
    descr<-c("\\subsubsection{Reading BIOM data}\n",
            "BIOM format is the standard for representing the taxa abundance profiles.",
            "The BIOM file is supported in both sparse or dense (without zeros) format.",
            "The metadata information can be present within a biom file or uploaded separately in plain text format.",
            "Metadata file contains additional information about samples such as experimental factors or sample",
            "grouping.\n");

  }else if(mbSetObj$dataSet$data.type=="mothur"){
    descr<-c("\\subsubsection{Reading mothur data}\n",
            "The mothur pipeline generate abundance and annotation information in its unique format.",
            "It contains these data in two separate plain text files (.shared and .taxonomy).",
            "Metadata file contains additional information about samples such as experimental factors or sample",
            "grouping.\n");
  }
        
  cat(descr, file=rnwFile, append=TRUE);
  cat("\n\n", file=rnwFile, append=TRUE);
        
  #user data info
  cat(mbSetObj$dataSet$read.msg, file=rnwFile, append=TRUE, sep="\n");
  cat(mbSetObj$dataSet$smpl.msg, file=rnwFile, append=TRUE, sep=" ");
        
  if(mbSetObj$module.type == "sdp"){
    cat("The genes are annotated as ", mbSetObj$dataSet$gene.id,"label.", file=rnwFile, append=TRUE, sep=" ");
  } else {
    cat("The OTUs are annotated as ", mbSetObj$dataSet$taxa.type,"label.", file=rnwFile, append=TRUE, sep=" ");
  }

  # the last step is sanity check
  descr<-c("\\subsubsection{Data Integrity Check}\n",
          "Before data analysis, a data integrity check is performed to make sure that all the necessary",
          "information has been collected. The sample variable should contain atleast two groups to perform most ",
          "of the comparative analysis. \\textit{By default, sample variables which are found to be constant and ",
          "continuous in nature will be removed from further analysis. Additionally, features just present in",
          "one sample will also be discarded from the data.\n}",
          paste("Figure", fig.count<<-fig.count+1,"shows the library size for inspection of each sample."),
          "\n ");
  cat(descr, file=rnwFile, append=TRUE);
        
  cmdhist<-c("\\begin{figure}[htp]",
            "\\begin{center}",
            paste("\\includegraphics[width=1.0\\textwidth]{", mbSetObj$imgSet$lib.size,"}", sep=""),
            "\\caption{Library size Overview}",      
            "\\end{center}",
            paste("\\label{",mbSetObj$imgSet$lib.size,"}", sep=""),
            "\\end{figure}");
        
  cat(cmdhist, file=rnwFile, append=TRUE);
  cat("\\clearpage", file=rnwFile, append=TRUE);
        
  if(!mbSetObj$module.type == "ppd"){
    # the data filtering
    # need to check if this process is executed
      if(is.null(mbSetObj$dataSet$filt.data)){
          errorMsg<- c("Error occured during filtering of your data ....",
                      "Fail to proceed. Please check if the data format you uploaded is correct.",
                      "Please visit our FAQs, Data Formats, and TroubleShooting pages for more information!");
          cat(errorMsg, file=rnwFile, append=TRUE);
          return();
      }
            
    descr<-c("\\subsubsection{Data Filtering}\n",
            "The purpose of the data filtering is to identify and remove features that are unlikely to be of",
            "use when modeling the data. No phenotype information are used in the filtering process, so the result",
            "can be used with any downstream analysis. This step can usually improves the results. \n",
            "Features having low count and variance can be removed during the filtration step.",               
            "Features having very few counts are filtered based on their abundance levels (minimum counts) across ",
            "samples (prevalence). Other than sample prevalence, such features can also be detected using minimum",
            " count cutoff based on their mean and median values.\n",
            "Features or taxa with constant or less variable abundances are invaluable for comparative analysis.",
            "Such features are filtered based on their inter-quantile ranges, standard deviations or coefficient of variations.",
            "\\textit{By default, features having zero counts across all the samples, or only appears in one sample",
            "will be removed from further analysis.}");
            
    cat(descr, file=rnwFile, append=TRUE);
    cat("\n\n", file=rnwFile, append=TRUE);

    filt.msg <- mbSetObj$dataSet$filt.msg;
            
    if(is.null(filt.msg)){
      filt.msg <- "No data filtering was performed.";
    }
    cat(filt.msg, file=rnwFile, append=TRUE);
  }
  cat("\n\n", file=rnwFile, append=TRUE);
}

# create normalization doc 
CreateNORMdoc <- function(mbSetObj){
  
  mbSetObj <- .get.mbSetObj(mbSetObj);

  # need to check if this process is executed
  if(is.null(mbSetObj$dataSet$norm.phyobj)){
    errorMsg<- c("Error occured during normalization of your data ....",
                 "Fail to proceed. Please check if the data format you uploaded is correct.",
                "Please visit our FAQs and Data Format pages for more information!");
            
    cat(errorMsg, file=rnwFile, append=TRUE);
    return();
  }
        
  descr1 <- c("\\subsection{Data Normalization}\n",
              "The data is stored as a table with one sample per column and one variable (taxa or ",
              "OTU) per row. The normalization procedures implemented below are grouped into three categories.",
              "Data rarefaction and scaling based methods deal with uneven sequencing depths by bringing",
              "samples to the same scale for comparison. While transformation based",
              "methods account for sparsity, compositionality, and large variations within the data.",
              "You can use one or combine all three to achieve better results. For more information",
              "about these methods, please refer to the paper by Weiss et al.\\footnote{Weiss et al.",
              "\\textit{Normalization and microbial differential abundance strategies depend upon",
              "data characteristics}, Microbiome 2017}",
              "\n",
              "The normalization consists of the following options:");

  cat(descr1, file=rnwFile, append=TRUE);
        
  descr2 <- c("\\begin{enumerate}",
              "\\item{Data rarefying (with or without replacement) }",
              "\\item{Data scaling: }",
              "\\begin{itemize}",
              "\\item{Total sum scaling (TSS) }",
              "\\item{Cumulative sum scaling (CSS) }",
              "\\item{Upper-quantile normalization (UQ)}",
              "\\end{itemize}",
              "\\item{Data transformation : }",
              "\\begin{itemize}",
              "\\item{Relative log expression (RLE) }",
              "\\item{Trimmed mean of M-values (TMM)}",
              "\\item{Centered log ratio (CLR)}",
              "\\end{itemize}",
              "\\end{enumerate}",
              "\n\n"
              );
  cat(descr2, file=rnwFile, append=TRUE, sep="\n");
  cat("\n\n", file=rnwFile, append=TRUE);

  norm.msg <- mbSetObj$dataSet$norm.msg;
        
  if(is.null(norm.msg)){
    norm.msg <- "No data normalization was performed.";
  }
  cat(norm.msg, file=rnwFile, append=TRUE);
  cat("\n", file=rnwFile, append=TRUE);
}

Init16SAnalMode<-function(){
        
  descr <- c("\\section{Marker Gene Analysis}",
             "MicrobiomeAnalyst offers a variety of methods commonly used in microbiome data analysis.",
             "They include:\n");
        
  cat(descr, file=rnwFile, append=TRUE, sep="\n");

  descr2 <- c("\\begin{enumerate}",
              "\\item{Visual exploration: }",
              "\\begin{itemize}",
              "\\item{Stacked bar/area plot }",
              "\\item{Rarefaction curve }",
              "\\item{Phylogenetic tree }",
              "\\item{Heat tree }",
              "\\end{itemize}",
              "\\item{Community profiling: }",
              "\\begin{itemize}",
              "\\item{Alpha diversity analysis }",
              "\\item{Beta Diversity analysis }",
              "\\item{Core microbiome analysis }",
              "\\end{itemize}",
              "\\item{Clustering analysis: }",
              "\\begin{itemize}",
              "\\item{Heatmap}",
              "\\item{Dendrogram}",
              "\\item{Correlation analysis}",
              "\\item{Pattern Search}",
              "\\end{itemize}",
              "\\item{Differential abundance analysis: }",
              "\\begin{itemize}",
              "\\item{Univariate analysis}",
              "\\item{metagenomeSeq}",
              "\\item{RNAseq methods}",
              "\\end{itemize}",
              "\\item{Biomarker analysis: }",
              "\\begin{itemize}",
              "\\item{LEfSe}",
              "\\item{Random Forests}",
              "\\end{itemize}",
              "\\item{Predictive functional profiling: }",
              "\\begin{itemize}",
              "\\item{PICRUSt}",
              "\\item{Tax4Fun}",
              "\\end{itemize}",
              "\\end{enumerate}");
        
  cat(descr2, file=rnwFile, append=TRUE, sep="\n");
  cat("\n\n", file=rnwFile, append=TRUE, sep="\n");
}

InitShotAnalMode<-function(){
        
  descr <- c("\\section{Shotgun data Profiling}",
             "MicrobiomeAnalyst offers a variety of methods commonly used in shotgun data analysis.",
             "They include:\n");
        
  cat(descr, file=rnwFile, append=TRUE, sep="\n");
  
  descr2 <- c("\\begin{enumerate}",
              "\\item{Overall functional profiling: }",
              "\\begin{itemize}",
              "\\item{Diversity overview }",
              "\\item{Association  analysis}",
              "\\end{itemize}",
              "\\item{Clustering analysis: }",
              "\\begin{itemize}",
              "\\item{Heatmap clustering}",
              "\\item{Dendrogram clustering}",
              "\\item{Correlation analysis}",
              "\\item{Pattern Search}",
              "\\end{itemize}",
              "\\item{Differential abundance analysis: }",
              "\\begin{itemize}",
              "\\item{Univariate analysis}",
              "\\item{metagenomeSeq}",
              "\\item{RNAseq methods}",
              "\\end{itemize}",
              "\\item{Biomarker analysis: }",
              "\\begin{itemize}",
              "\\item{LEfSe}",
              "\\item{Random Forests}",
              "\\end{itemize}",
             "\\end{enumerate}");
        
  cat(descr2, file=rnwFile, append=TRUE, sep="\n");
  cat("\n\n", file=rnwFile, append=TRUE, sep="\n");
}

CreateAnalNullMsg<-function(){
  descr <- c("No analysis was performed on your data.\n");
  cat(descr, file=rnwFile, append=TRUE, sep="\n");
}

# create Stacked Area/Bar plot/Piechart doc
CreateVISEXPLRdoc <- function(mbSetObj){
  
  mbSetObj <- .get.mbSetObj(mbSetObj);

  # need to check if this process is executed
  if(is.null(mbSetObj$analSet$stack) & is.null(mbSetObj$analSet$pie)){
    return();
  }
        
  descr <- c("\\subsection{Visual Exploration}\n",
             "These methods are used to visualize the taxonomic composition of community through direct quantitative comparison of abundances.",
             "MicrobiomeAnalyst provides an option to view this composition at various taxonomic levels (phylum, class, order) using either",
             "stacked bar/stacked area plot or piechart. Viewing composition at higher-levels (phylum) provides a better picture than lower-levels (species)",
             "when the number of species in a community is large and diversified.\n",
             "Additionally, such taxonomic abundance or composition can be viewed at community-level (all samples), sample-group",
             "level (based on experimental factor) or at individual sample-level.\n",
             "Taxa with very low read counts can also be collapsed into \\texttt{Others} category using a count cutoff",
             "based on either sum or median of their counts across all samples or all groups. Merging such minor",
             "taxa will help in better visualization of significant taxonomic patterns in data.",
             "\n");

  cat(descr, file=rnwFile, append=TRUE);

  # STACKPLOT
  if(!is.null(mbSetObj$analSet$stack)){
            
    descr<- paste("Figure", fig.count<<-fig.count+1,"shows the taxonomic composition using Stacked bar/area plot.");
            
    cat(descr, file=rnwFile, append=TRUE);
            
    cmdhist<-c("\\begin{figure}[htp]",
               "\\begin{center}",
               paste("\\includegraphics[width=1.0\\textwidth]{", mbSetObj$imgSet$stack,"}", sep=""),
               paste("\\caption{Taxonomic composition of", 
               " community at ", "\\texttt{", mbSetObj$analSet$stack.taxalvl, "} level using ", "\\texttt{", mbSetObj$analSet$plot, "} plot}", sep=""),      
               "\\end{center}",
               paste("\\label{",mbSetObj$imgSet$stack,"}", sep=""),
               "\\end{figure}");
            
    cat(cmdhist, file=rnwFile, append=TRUE);
  }

  # PIECHART  
  if(!is.null(mbSetObj$analSet$pie)){ 
    descr <- paste("Figure", fig.count<<-fig.count+1,"shows the taxonomic composition using piechart.");
    cat(descr, file=rnwFile, append=TRUE);
            
    cmdhist<-c("\\begin{figure}[htp]",
               "\\begin{center}",
               paste("\\includegraphics[width=1.0\\textwidth]{", mbSetObj$imgSet$pie,"}", sep=""),
               paste("\\caption{Taxonomic composition of", 
               " community at ", "\\texttt{", mbSetObj$analSet$pie.taxalvl, "} level using piechart}", sep=""),      
               "\\end{center}",
               paste("\\label{",mbSetObj$imgSet$pie,"}", sep=""),
               "\\end{figure}");
            cat(cmdhist, file=rnwFile, append=TRUE);
  }
  cat("\\clearpage", file=rnwFile, append=TRUE);
}

###################################
# Rarefaction curve analysis#######
####################################


CreateRAREFCTIONCURVEdoc <- function(mbSetObj){
  
  mbSetObj <- .get.mbSetObj(mbSetObj);
  
  # need to check if this process is executed
  if(is.null(mbSetObj$analSet$rarefaction_curve)){
    return();
  }
  
  descr <- c("\\subsection{Rarefaction curve analysis }\n",
             "This method is used to present relationship between number of OTUs and number of sequences",
             "It can infer if the reads of a sample is enough to reach plateau, which means that with increasing of sequences, the gain of newly discovered OTUs is limited.",
             "If sequence depth of some samples are not enough, you may consider to resequence these samples or removed from downstream analysis.",
             "User can choose different metadata variables as group, line colors and line types.",
             "Rarefaction curve analysis is performed using the modified function ggrare originated from \\texttt{ranacapa} package\\footnote{Gaurav S. Kandlikar",
             "\\textit{ranacapa: An R package and Shiny web app to explore environmental DNA data with exploratory statistics and interactive visualizations.}, 2018.}\n",
             paste("Figure", fig.count<<-fig.count+1,"shows the rarefaction curve."),
             "\n "
  );
  cat(descr, file=rnwFile, append=TRUE,sep="\n");

   cmdhist<-c("\\begin{figure}[htp]",
               "\\begin{center}",
               paste("\\includegraphics[width=1.0\\textwidth]{", mbSetObj$analSet$rarefaction_curve,"}", sep=""),
               paste("\\caption{Rarefaction curve", 
                     " using ", "\\texttt{", mbSetObj$analSet$rarefaction_curve_data.src, "} dataset }", sep=""),      
               "\\end{center}",
               paste("\\label{", mbSetObj$analSet$rarefaction_curve,"}", sep=""),
               "\\end{figure}");
    cat(cmdhist, file=rnwFile, append=TRUE);

  cat("\\clearpage", file=rnwFile, append=TRUE);
}


###################################
# Phylogenetic tree analysis#######
####################################


CreatePHYLOGENETICTREEdoc <- function(mbSetObj){
  
  mbSetObj <- .get.mbSetObj(mbSetObj);
  
  # need to check if this process is executed
  if(is.null(mbSetObj$analSet$phylogenetic_tree_curve)){
    return();
  }
  
  descr <- c("\\subsection{Phylogenetic tree }\n",
             "This method is used to determine the evolutionary relationship among taxonomic groups",
             "A phylogenetic tree is a diagram which represents evolutionary relationships among species.",
             "It reflects how they evolved from common ancestors. If two organisms share a more recent common ancestor, they are more related.",
             "Therefore, phylogenetic tree can be used to represent distances between them, which will be used for UniFrac distance based analysis, such as phylogenetic beta diversity.",
             "Phylogenetic tree analysis is performed using R package \\texttt{phyloseq} package\\footnote{McMurdie",
             "\\textit{ phyloseq: An R Package for Reproducible Interactive Analysis and Graphics of Microbiome Census Data.}, 2013.}\n",
             paste("Figure", fig.count<<-fig.count+1,"shows the phylogenetic tree."),
             "\n ");
  cat(descr, file=rnwFile, append=TRUE,sep="\n");

   cmdhist<-c("\\begin{figure}[htp]",
               "\\begin{center}",
               paste("\\includegraphics[width=1.0\\textwidth]{", mbSetObj$analSet$phylogenetic_tree_curve,"}", sep=""),
               paste("\\caption{Phylogenetic tree", 
                     " at ", "\\texttt{", mbSetObj$analSet$phylogenetic_tree_curve_tax_level, "} level }", sep=""),      
               "\\end{center}",
               paste("\\label{", mbSetObj$analSet$phylogenetic_tree_curve,"}", sep=""),
               "\\end{figure}");
    cat(cmdhist, file=rnwFile, append=TRUE);

  cat("\\clearpage", file=rnwFile, append=TRUE);
}

###################################
# Heat tree analysis#######
####################################


CreateHEATTREEdoc <- function(mbSetObj){
  
  mbSetObj <- .get.mbSetObj(mbSetObj);
  # need to check if this process is executed
  if(is.null(mbSetObj$analSet$heat_tree_plot)){
    return();
  }
  #print(mbSetObj$analSet$heat_tree_tax);print(mbSetObj$analSet$heat_tree_meta);print(mbSetObj$analSet$heat_tree_comparison);
descr <- c("\\subsection{Heat tree }\n",
           "This method is used to compare abundance of different taxonomic levels for each pair of factors in a metadata variable.",
           "Heat Tree uses hierarchical structure of taxonomic classifications to quantitatively (median abundance) and statistically (non-parameter Wilcoxon Rank Sum test ) depict taxon differences among communities.",
           "It generates a differential heat tree to show which taxa are more abundance in each treatment.",
           "User can choose different taxonomic levels from phylum to species.",
           "User can also select the metadata variable and subsequently determined which pair of factors in that metadata variable they want to compare.",
           "Heat tree analysis is performed using R package \\texttt{metacoder} package\\footnote{Zachary S. L. Foster",
           "\\textit{ metacoder: An R package for visualization and manipulation of community taxonomic diversity data.}, 2017, R package version 0.3.2.}\n",
           paste("Figure", fig.count,"shows the heat tree for pair-wise comparison."),
           "\n ");

  cat(descr, file=rnwFile, append=TRUE,sep="\n");

cmdhist<-c("\\begin{figure}[htp]",
           "\\begin{center}",
           paste("\\includegraphics[width=1.0\\textwidth]{", mbSetObj$analSet$heat_tree_plot,"}", sep=""),
           paste("\\caption{Phylogenetic tree", 
                 " at ", "\\texttt{", mbSetObj$analSet$heat_tree_tax, "} level ",
                 "in group ", "\\texttt{", mbSetObj$analSet$heat_tree_meta, "}", " }", sep=""),      
           "\\end{center}",
           paste("\\label{", mbSetObj$analSet$heat_tree_plot,"}", sep=""),
           "\\end{figure}");

  cat(cmdhist, file=rnwFile, append=TRUE);

  cat("\\clearpage", file=rnwFile, append=TRUE);
}

# Alpha-diversity
CreateALPHDIVdoc <- function(mbSetObj){
  
  mbSetObj <- .get.mbSetObj(mbSetObj);
  
  # need to check if this process is executed
  if(is.null(mbSetObj$analSet$alpha)){
    return();
  }
  
  descr <- c("\\subsection{Alpha diversity analysis }\n",
             "This method is used to measure the diversity present within a sample or community.",
             "Alpha diversity can be characterized via the total number of species (richness),",
             "the abundances of the species (evenness) or measures that considered",
             "both richness and evenness. How these measures estimates the diversity",
             "is need to be considered when performing alpha-diversity analysis.",
             "User can choose from richness based measure such as Observed index which calculates the actual number of ",
             "unique taxa observed in each sample. While the Chao1 and ACE measures",
             "estimate the richness by inferring out the number of rare organisms that may",
             "have lost due to undersampling. Also,there are indices such as Shannon, Simpson",
             "and Fisher in which along with the number (richness), the abundance of organisms (evenness)",
             "is also measured to describe the actual diversity of a community.\n",
             "Alpha diversity analysis is performed using the \\texttt{phyloseq} package\\footnote{Paul J. McMurdie",
             "\\textit{phyloseq: An R package for reproducible interactive analysis and graphics of microbiome census data.}, 2013, R package version 1.19}.",
             "The results are plotted across samples and reviewed as box plots for each group or experimental factor.",
             "Further, the statistical significance of grouping based on experimental factor is also estimated using",
             "either parametric or non-parametric test.",
             paste("Figure", fig.count<<-fig.count+1,"shows the alpha diversity measure across all the samples for given diversity index."),
             paste("Figure", fig.count<<-fig.count+1,"shows the diversity distribution using box plot for a given group or experimental factor.\n"));
  cat(descr, file=rnwFile, append=TRUE,sep="\n");
  cat("\n\n", file=rnwFile, append=TRUE);
  
  cmdhist<-c(
    "\\begin{figure}[htp]",
    "\\begin{center}",
    paste("\\includegraphics[width=1.0\\textwidth]{", mbSetObj$imgSet$alpha,"}", sep=""),
    "\\caption{Alpha-diversity measure using","\\texttt{",mbSetObj$analSet$alpha.meth, "} at","\\texttt{",mbSetObj$analSet$alpha.taxalvl, "} level across all the samples.",
                    "The samples are represented on X-axis and their estimated diversity on Y-axis. Each sample
                     is colored based on","\\texttt{",mbSetObj$analSet$alpha.metadata, "} class}",
    "\\end{center}",
    paste("\\label{",mbSetObj$imgSet$alpha,"}", sep=""),
    "\\end{figure}"
    );
    cat(cmdhist, file=rnwFile, append=TRUE);
  cmdhist2<-c( 
    "\\begin{figure}[htp]",
    "\\begin{center}",
    paste("\\includegraphics[width=1.0\\textwidth]{", mbSetObj$imgSet$alpha.box,"}", sep=""),
    "\\caption{Alpha-diversity measure using","\\texttt{",mbSetObj$analSet$alpha.meth, "} at","\\texttt{", mbSetObj$analSet$alpha.taxalvl, "} level represented as boxplot.",
             "Each boxplot represents the diversity distribution of a group present within", "\\texttt{",mbSetObj$analSet$alpha.metadata, "} class
              [Statistical significance:" ,"\\texttt{", mbSetObj$analSet$alpha.stat.info,"}]}",
    "\\end{center}",
    paste("\\label{",mbSetObj$imgSet$alpha.box,"}", sep=""),
    "\\end{figure}"
  );
  cat(cmdhist2, file=rnwFile, append=TRUE, sep="\n");
  cat("\\clearpage", file=rnwFile, append=TRUE, sep="\n");
}

# Beta-diversity
CreateBETADIVdoc<-function(mbSetObj){
  
  mbSetObj <- .get.mbSetObj(mbSetObj);
  
  # need to check if this process is executed
  if(is.null(mbSetObj$analSet$beta)){
    return();
  }
  
  descr <- c("\\subsection{Beta diversity Analysis}\n",
             "This method provides a way to compare the diversity or composition between two samples or microbial communities.",
             "These methods compare the changes in the presence/absence or abundance of thousands of taxa present in a dataset",
             "and summarize these into how 'similar' or 'dissimilar' two samples. Each sample gets compared to every other sample", 
             "generating a distance or dissimilarity matrix. Two parameters need to be considered when performing",
             "beta diversity analysis. The first one is how similarity or distance between sample is measured which includes",
             "non-phylogenetic (Bray-Curtis distance, Shannon index, Jaccard index) and phylogenetic-based (weighted and unweighted UniFrac)",
             "distances. The other parameter is how to visualize such dissimilarity matrix in lower dimensions. Ordination-based methods such as Principle Coordinate",
             "Analysis (PCoA) and non-metric multidimensional scaling (NMDS) are used to visualize these matrix in 2 or 3-D plot where",
             "each point represents the entire microbiome of a single sample. Each axis reflects the percent of the variation between", 
             "the samples with the X-axis representing the highest dimension of variation and the Y-axis representing the second",
             "highest dimension of variation. Further,each point or sample displayed on PCoA or NMDS plots is colored based on either",
             "sample group, features alpha diversity measures, or the abundance levels of a specific feature.\n",
             "Also, the statistical significance of the clustering pattern in ordination plots can be evaluated using anyone among",
             "Permutational ANOVA (PERMANOVA), Analysis of group Similarities (ANOSIM) and Homogeneity of Group Dispersions (PERMDISP).",
             "\n\n",
             "Beta diversity analysis is performed using the \\texttt{phyloseq} package\\footnote{Paul J. McMurdie",
             "\\textit{phyloseq: An R package for reproducible interactive analysis and graphics of microbiome census data.}, 2013, R package version 1.19}.",
             paste("Figure", fig.count<<-fig.count+1,"shows the ordination plot represented in 2-D;"),
             paste("Statistical significance is found out using","\\texttt{", mbSetObj$analSet$stat.info, "}."));
  
  cat(descr, file=rnwFile, append=TRUE,sep="\n");
  
  cmdhist<-c(
    "\\begin{figure}[htp]",
    "\\begin{center}",
    paste("\\includegraphics[width=1.0\\textwidth]{", mbSetObj$imgSet$beta2d,"}", sep=""),
    "\\caption{2-D","\\texttt{", mbSetObj$analSet$beta.meth, "} plot using \\texttt{", mbSetObj$analSet$beta.dist, "} distance. The explained",
                          "variances are shown in brackets.}",
    "\\end{center}",
    paste("\\label{",mbSetObj$imgSet$beta2d,"}", sep=""),
    "\\end{figure}"
  );
  cat(cmdhist, file=rnwFile, append=TRUE, sep="\n");
  cat("\n\n", file=rnwFile, append=TRUE, sep="\n");
}



# Hierarchical clustering
CreateHCdoc <- function(mbSetObj){
  
  mbSetObj <- .get.mbSetObj(mbSetObj);
  
  # need to check if this process is executed
  if(is.null(mbSetObj$analSet$tree) & is.null(mbSetObj$analSet$heatmap)){
    return();
  }
  
  descr <- c("\\subsection{Hierarchical Clustering}\n",
             "In hierarchical cluster analysis, each sample begins",
             "as a separate cluster and the algorithm proceeds to combine them until all",
             "samples belong to one cluster. Two parameters need to be considered when performing",
             "hierarchical clustering. The first one is how similarity or distance between sample is measured which includes Bray-Curtis distance,",
             "Shannon index, Jaccard index, weighted and unweighted UniFrac. The other parameter is clustering",
             "algorithms, including average linkage (clustering uses the centroids of the observations),",
             "complete linkage (clustering uses the farthest pair of observations between the two groups),",
             "single linkage (clustering uses the closest pair of observations) and Ward's linkage",
             "(clustering to minimize the sum of squares of any two clusters). In MicrobiomeAnalyst, the result of clustering analysis are",
             "supported using Heatmap and dendrogram.\n",
             "Hierarchical clustering is performed with the \\texttt{hclust} function in package \\texttt{stat}."
  );
    
  cat(descr, file=rnwFile, append=TRUE,sep="\n");
  
  if(!is.null(mbSetObj$analSet$tree)){
    descr<-paste("Figure", fig.count<<-fig.count+1,"shows the clustering result in the form of a dendrogram.");
    cat(descr, file=rnwFile, append=TRUE);
    #taxonomic class will be replaced with gene id in SDP
    
    if(mbSetObj$module.type == "sdp"){
      mbSetObj$analSet$tree.taxalvl <- mbSetObj$dataSet$gene.id;
    }
    
    cmdhist<-c(
        "\\begin{figure}[htp]",
        "\\begin{center}",
        paste("\\includegraphics[width=1.0\\textwidth]{", mbSetObj$imgSet$tree,"}", sep=""),
        paste("\\caption{Clustering result shown as dendrogram (", 
            "distance measure using ","\\texttt{",mbSetObj$analSet$tree.dist, "} and clustering algorithm using ","\\texttt{", mbSetObj$analSet$tree.clust, "} at ","\\texttt{", mbSetObj$analSet$tree.taxalvl,"} level)}", sep=""),
        "\\end{center}",
        paste("\\label{",mbSetObj$imgSet$tree,"}", sep=""),
        "\\end{figure}"
    );
    cat(cmdhist, file=rnwFile, append=TRUE, sep="\n");
  }
  
  if(!is.null(mbSetObj$analSet$heatmap)){
    descr <- paste("Figure", fig.count<<-fig.count+1,"shows the clustering result in the form of a heatmap.");
    cat(descr, file=rnwFile, append=TRUE);
    #taxonomic class will be replaced with gene id in SDP
    if(mbSetObj$module.type == "sdp"){
      mbSetObj$analSet$heat.taxalvl <- mbSetObj$dataSet$gene.id;
    }
    cmdhist<-c(
        "\\begin{figure}[htp]",
        "\\begin{center}",
        paste("\\includegraphics[width=1.0\\textwidth]{", mbSetObj$imgSet$heatmap,"}", sep=""),
        paste("\\caption{Clustering result shown as heatmap (", 
            "distance measure using ","\\texttt{", mbSetObj$analSet$heatmap.dist, "} and clustering algorithm using ","\\texttt{", mbSetObj$analSet$heatmap.clust, "} at ","\\texttt{", mbSetObj$analSet$heat.taxalvl,"} level)}", sep=""),
        "\\end{center}",
        paste("\\label{",mbSetObj$imgSet$heatmap,"}", sep=""),
        "\\end{figure}"
    );                    
    cat(cmdhist, file=rnwFile, append=TRUE, sep="\n");
  }
   cat("\\clearpage", file=rnwFile, append=TRUE);
}

# Core microbiome analysis
CreateCOREdoc <- function(mbSetObj){
  
  mbSetObj <- .get.mbSetObj(mbSetObj);
  
  # need to check if this process is executed
  if(is.null(mbSetObj$analSet$core)){
    return();
  }
  
  if(isEmptyMatrix(mbSetObj$analSet$core)){
      core.tab<-NULL;
  }else{
      core.tab<-paste("Table", table.count<<-table.count+1,"shows the details of these features.");
  }
  
  descr <- c("\\subsection{Core Microbiome analysis}\n",
             "This method helps in identifying core taxa or features that remain unchanged in their",
             "composition across the whole microbial community. Two parameters need to be considered when performing",
             "core microbiome analysis. The first one is sample prevalence which is defined as minimum fractions (percentage) of",
             "samples that a taxa or feature must be observed. The other parameter is ",
             "the relative abundance (fractions) of a taxa or features in order to consider them as a part of core member.",
             "\n",
             "Core microbiome analysis is adopted from the \\texttt{core} function in R package \\texttt{microbiome}.",
             "The result of this analysis is represented in the form of heatmap of core taxa or features where",
             "Y-axis represent the prevalence level of core features across the detection threshold (Relative abundance) range on X-axis.",
             paste("Figure", fig.count<<-fig.count+1,"shows the core microbiome result represented in the form of a heatmap."),
             core.tab,
             "\n"
  );
  
  cat(descr, file=rnwFile, append=TRUE,sep="\n");
  
  if(!is.null(mbSetObj$analSet$core)){
    cmdhist<-c("\\begin{figure}[htp]",
               "\\begin{center}",
               paste("\\includegraphics[width=1.0\\textwidth]{", mbSetObj$imgSet$core,"}", sep=""),
               paste("\\caption{Heatmap representing the core microbiome ", 
               "at ","\\texttt{", mbSetObj$analSet$core.taxalvl, "} level}", sep=""),
               "\\end{center}",
               paste("\\label{",mbSetObj$imgSet$core,"}", sep=""),
               "\\end{figure}"); 
        
    cat(cmdhist, file=rnwFile, append=TRUE, sep="\n");
        
  }
  cat("\\clearpage", file=rnwFile, append=TRUE, sep="\n");
}


# create Correlation doc
CreateCORRdoc <- function(mbSetObj){
  
  mbSetObj <- .get.mbSetObj(mbSetObj);

  # need to check if this process is executed
  if(is.null(mbSetObj$analSet$cor.mat) & is.null(mbSetObj$analSet$cor.heatmat)){
    return();
  }

  # need to check if this process is executed
  if(!is.null(mbSetObj$analSet$cor.heatmat)){
    descr <- c("\\subsection{Correlation Analysis}\n",
               "Correlation analysis can be used to visualize the overall correlations between different features",
               "It can also be used to identify which features are correlated with a feature of interest.",
               "Correlation analysis can also be used to identify if certain features show particular patterns",
               "under different conditions. Users first need to define a pattern in the form of a series of hyphenated numbers.",
               "For example, in a time-series study with four time points, a pattern of of",
               "\\texttt{1-2-3-4} is used to search taxa with increasing the abundance as",
               "time changes; while a pattern of \\texttt{3-2-1-3} can be used to search feature",
               "that decrease at first, then bounce back to the original abundance. \n",
               paste("Figure", fig.count<<-fig.count+1, "shows the overall correlation heatmap."),
               "\n");

    cat(descr, file=rnwFile, append=TRUE);
    #taxonomic class will be replaced with gene id in SDP
    if(mbSetObj$module.type == "sdp"){
      mbSetObj$analSet$corheat.taxalvl<-mbSetObj$dataSet$gene.id;
    }
            
    cmdhist<-c("\\begin{figure}[htp]",
               "\\begin{center}",
               paste("\\includegraphics[width=1.0\\textwidth]{", mbSetObj$imgSet$cor.heat,"}", sep=""),
               "\\caption{Correlation Heatmap at","\\texttt{", mbSetObj$analSet$corheat.taxalvl, "}level}",
               "\\end{center}",
               paste("\\label{",mbSetObj$imgSet$cor.heat,"}", sep=""),
               "\\end{figure}");

    cat(cmdhist, file=rnwFile, append=TRUE, sep="\n");
    cat("\n\n", file=rnwFile, append=TRUE, sep="\n");
  }

  if(!is.null(mbSetObj$analSet$cor.mat)){
    if(isEmptyMatrix(mbSetObj$analSet$cor.mat)){
      cor.tab<-NULL;
    }else{
      cor.tab<-paste("Table", table.count<<-table.count+1,"shows the details of these features.");
    }

    descr <- c("\\subsection{Pattern Search}\n",
               "Paatern Search can be used to identify which features are correlated with a feature of interest.",
               "Basically,correlation analysis is used to identify if certain features show particular patterns",
               "under different conditions. Users first need to define a pattern in the form of a series of hyphenated numbers.",
               "For example, in a time-series study with four time points, a pattern of of",
               "\\texttt{1-2-3-4} is used to search taxa with increasing the abundance as",
               "time changes; while a pattern of \\texttt{3-2-1-3} can be used to search feature",
               "that decrease at first, then bounce back to the original abundance level.",
               "\n\n",
               paste("Figure", fig.count<<-fig.count+1, "shows the important features identified by correlation analysis."),
               cor.tab,
               "\n");

    cat(descr, file=rnwFile, append=TRUE);
      
    #taxonomic class will be replaced with gene id in SDP
    if(mbSetObj$module.type == "sdp"){
      mbSetObj$analSet$corph.taxalvl<-mbSetObj$dataSet$gene.id;
    }
            
    cmdhist<-c("\\begin{figure}[htp]",
               "\\begin{center}",
               paste("\\includegraphics[width=1.0\\textwidth]{", mbSetObj$imgSet$cor.ph,"}", sep=""),
               "\\caption{Important features (at","\\texttt{", mbSetObj$analSet$corph.taxalvl, "} level) selected by correlation analysis with light",
               "purple indicates positive correlation and blue indicate negative correlations.}",
               "\\end{center}",
               paste("\\label{",mbSetObj$imgSet$cor.ph,"}", sep=""),
               "\\end{figure}");
            cat(cmdhist, file=rnwFile, append=TRUE,sep="\n");
            cmdhist2<-c("<<echo=false, results=tex>>=",
                   "GetSigTable.Corr(mbSet)",
                   "@");
            cat(cmdhist2, file=rnwFile, append=TRUE,sep="\n");
    }
  cat("\\clearpage", file=rnwFile, append=TRUE);
}


# create Univariate doc
CreateUNIVARdoc<-function(mbSetObj){
  
  mbSetObj <- .get.mbSetObj(mbSetObj);

  # need to check if this process is executed
  if(is.null(mbSetObj$analSet$Univar$resTable)){
    return();
  }

  if(isEmptyMatrix(mbSetObj$analSet$Univar$resTable)){
    univar.tab<-NULL;
  }else{
    #taxonomic class will be replaced with gene id in SDP
    if(mbSetObj$module.type == "sdp"){
      mbSetObj$analSet$univar.taxalvl<-mbSetObj$dataSet$gene.id;
    }
    univar.tab<-paste("Table", table.count<<-table.count+1, "shows the important features identified by Univariate analysis at","\\texttt{", mbSetObj$analSet$univar.taxalvl, "}level");
  }
        
  descr <- c("\\subsection{Univariate analysis}\n",
             "Univariate analysis methods are the most common methods used for exploratory data analysis.",
             "This method is used of identify differentially abundant features in microbiome data analysis.",
             "MicrobiomeAnalyst provides an option to perform such statistical comparisons using either common",
             "parametric (T-test/ANOVA) or non-parametric (Mann-Whitney/Kruskal-Wallis) tests based on single",
             "grouping variable or experimental factor.\n",
             "Features are considered to  be significant based on their",
             "adjusted p-value. The default is \\texttt{adj.p-value cutoff} = 0.05.",
             "When the number of samples is high (>50 samples), rarefying or proportion normalized data paired",
             "with a non-parametric test, yield as high sensitivity as other methods for identifying differential features.",
             "\n\n",
             univar.tab,
             "\n");

  cat(descr, file=rnwFile, append=TRUE);
  cat("\n\n", file=rnwFile, append=TRUE);
         
  cmdhist<-c("<<echo=false, results=tex>>=",
             "GetSigTable.UNIVAR(mbSet)",
             "@");
  cat(cmdhist, file=rnwFile, append=TRUE, sep="\n");
  cat("\\clearpage", file=rnwFile, append=TRUE);
}

# create MetagenomeSeq doc
CreateMETAGENOSEQdoc <- function(mbSetObj){
  
  mbSetObj <- .get.mbSetObj(mbSetObj);

  # need to check if this process is executed
  if(is.null(mbSetObj$analSet$metagenoseq$resTable)){
    return();
  }

  if(isEmptyMatrix(mbSetObj$analSet$metagenoseq$resTable)){
    metagenoseq.tab<-NULL;
  }else{
    #taxonomic class will be replaced with gene id in SDP
    if(mbSetObj$module.type == "sdp"){
      mbSetObj$analSet$metageno.taxalvl<-mbSetObj$dataSet$gene.id;
    }
      metagenoseq.tab<-paste("Table", table.count<<-table.count+1,"shows the important features identified by metagenomeSeq at","\\texttt{", mbSetObj$analSet$metageno.taxalvl, "}level");
  }

  descr <- c("\\subsection{metagenomeSeq}\n",
             "This method is specifically designed to evaluate differential abundance ",
             "in sparse marker-gene survey data. This method combines Cumulative Sum",
             "Scaling (CSS) normalization with zero-inflated Log-Normal mixture model",
             "(fitFeature) or zero-inflated Gaussian (fitZIG) distribution mixture to account",
             "for undersampling and sparsity in OTU count data. In the fitZIG model, the count",
             "distribution is modeled as a mixture of two distributions. The zeros present in count",
             "data are modeled using point mass at zero, while remaining log transformed counts follows",
             "a normal distribution. On the other hand, fitFeature model shapes the count distribution using",
             "zero-inflated lognormal model. fitFeature model can be only used and recommended for two groups,",
             "whereas fitZIG model works on multiple groups for differential abundance testing. This method",
             "outperforms other methods in the detection of differentially abundant rare features. Features are",
             "considered to  be significant based on their adjusted p-value. The default is",
             "\\texttt{adj.p-value cutoff} = 0.05.\n\n",
             "This method is available as a \\texttt{metagenomeSeq} R package \\footnote{JN Paulson \\textit{metagenomeSeq: Statistical",
             "analysis for sparse high-throughput sequencing},2013, R Bioconductor package version 1.20.1}.",
             metagenoseq.tab,
             "\n");

  cat(descr, file=rnwFile, append=TRUE);
  cat("\n\n", file=rnwFile, append=TRUE);
  cmdhist<-c("<<echo=false, results=tex>>=",
             "GetSigTable.METAGENOSEQ(mbSet)",
             "@");
  cat(cmdhist, file=rnwFile, append=TRUE, sep="\n");
  cat("\\clearpage", file=rnwFile, append=TRUE);
}


# RNASeq method
CreateRNASEQdoc <- function(mbSetObj){
  
  mbSetObj <- .get.mbSetObj(mbSetObj);

  # need to check if this process is executed
  if(is.null(mbSetObj$analSet$rnaseq$resTable)){
    return();
  };
        
  #taxonomic class will be replaced with gene id in SDP
  if(mbSetObj$module.type == "sdp"){
    mbSetObj$analSet$rnaseq.taxalvl<-mbSetObj$dataSet$gene.id;
  }
        
  descr <- c("\\subsection{RNAseq methods}\n",
             "MicrobiomeAnalyst supports two powerful statistical methods including EdgeR and DESEQ2",
             "for performing differential abundance analysis. Both these methods were originally developed",
             "for RNA-Seq count data. However, these methods also seem to perform well or better than other statistical algorithms",
             "for metagenomic datasets.\\footnote{McMurdie P.J. \\textit{Waste not, want not: why rarefying microbiome data is inadmissible.}",
             "PLoS Computational Biology 2014.} They differ in their approach of data normalization and the algorithms used for evaluation",
             "of variance or dispersion. EdgeR utilizes RLE (Relative log expression) as a default normalization and assumes a ",
             "Negative Binomial model for count distributions. DESeq2 variance estimations are based on modeling the counts to",
             "Negative Binomial Generalized Linear Model. Usually, edgeR is more powerful (detecting more differential features) but",
             "also yields higher false positives. Whereas, DESeq2 is more robust and powerful in identifying the differential ",
             "features (i.e., lower false positives). However, it is more computationally demanding and take a considerable",
             "amount of time for big data sets.\n",
             "Features are considered to  be significant based on their adjusted p-value. The default is",
             "\\texttt{adj.p-value cutoff} = 0.05.",
             "\n\n",
             paste("Table", table.count<<-table.count+1,"shows the important features identified by","\\texttt{", mbSetObj$analSet$rnaseq.meth, "}method at", "\\texttt{", mbSetObj$analSet$rnaseq.taxalvl, "}level" ));

  cat(descr, file=rnwFile, append=TRUE);
  cat("\n\n", file=rnwFile, append=TRUE);
  cmdhist<-c("<<echo=false, results=tex>>=",
             "GetSigTable.RNASeq(mbSet)",
             "@");
        
  cat(cmdhist, file=rnwFile, append=TRUE, sep="\n");
  cat("\\clearpage", file=rnwFile, append=TRUE);
}

# LEFSE
CreateLEFSEdoc <- function(mbSetObj){
  
  mbSetObj <- .get.mbSetObj(mbSetObj);
    
  # need to check if this process is executed
  if(is.null(mbSetObj$analSet$lefse$resTable)){
    return();
  };
    
  #taxonomic class will be replaced with gene id in SDP
  if(mbSetObj$module.type == "sdp"){
    mbSetObj$analSet$lefse.taxalvl<-mbSetObj$dataSet$gene.id;
  }

  descr <- c("\\subsection{LDA Effect Size (LEfSe)}\n",
             "This method is specifically designed for biomarker discovery and explanation in high-dimensional metagenomic data.",
             "\\footnote{N Segata \\textit{Metagenomic biomarker discovery and explanation.}",
             "Nature Methods 2013.}",
             "It incorporates statistical significance with biological consistency (effect size) estimation.",
             "It performs non-parametric factorial Kruskal-Wallis (KW) sum-rank test to identify features with significant",
             "differential abundance with regard to experimental factor or class of interest, followed by Linear Discriminant",
             "Analysis (LDA) to calculate the effect size of each differentially abundant features. The result consists of all the",
             "features, the logarithmic value of the maximum mean among all the groups or classes, and if the features are",
             "differentially significant, the group with the highest mean and the logarithmic LDA score (Effect Size).\n",
             "Features are considered to  be significant based on their adjusted p-value. The default is",
             "\\texttt{adj.p-value cutoff} = 0.05.",
             "\n\n",
             paste("Table", table.count<<-table.count+1,"shows the important features identified by LEfSe at","\\texttt{", mbSetObj$analSet$lefse.taxalvl, "}level"));

  cat(descr, file=rnwFile, append=TRUE);
  cat("\n\n", file=rnwFile, append=TRUE);
  cmdhist<-c("<<echo=false, results=tex>>=",
             "GetSigTable.LEFSE(mbSet)",
             "@");
  cat(cmdhist, file=rnwFile, append=TRUE, sep="\n");
  cat("\\clearpage", file=rnwFile, append=TRUE);  

  cmdhist<-c("\\begin{figure}[htp]",
             "\\begin{center}",
             paste("\\includegraphics[width=1.0\\textwidth]{", mbSetObj$analSet$lefse_plot,"}", sep=""),
             paste("\\caption{Graphical summary", 
                   " at ", "\\texttt{", mbSetObj$analSet$lefse.taxalvl, "} level ",
                   "in group ", "\\texttt{", mbSetObj$analSet$meta, "}", " }", sep=""),      
             "\\end{center}",
             paste("\\label{", mbSetObj$analSet$lefse_plot,"}", sep=""),
             "\\end{figure}");
  
  cat(cmdhist, file=rnwFile, append=TRUE);
  cat("\\clearpage", file=rnwFile, append=TRUE);
}

# random forests
CreateRFdoc <- function(mbSetObj){
  
  mbSetObj <- .get.mbSetObj(mbSetObj);

  # need to check if this process is executed
  if(is.null(mbSetObj$analSet$rf)){
    return();
  }

  descr <- c("\\subsection{Random Forest (RF) }\n",
             "Random Forest is a supervised learning algorithm suitable for high dimensional data analysis.",
             "It uses an ensemble of classification trees, each of which is grown by random feature",
             "selection from a bootstrap sample at each branch. Class prediction is based on the",
             "majority vote of the ensemble. RF also provides other useful information such as OOB",
             "(out-of-bag) error and variable importance measure. During tree construction, about",
             "one-third of the instances are left out of the bootstrap sample. This OOB data",
             "is then used as test sample to obtain an unbiased estimate of the classification",
             "error (OOB error). Variable importance is evaluated by measuring the increase of the",
             "OOB error when it is permuted. The outlier measures are based on the proximities during tree construction.",
             "\n\n",
             "RF analysis is performed using the \\texttt{randomForest} package\\footnote{Andy Liaw and",
             "Matthew Wiener. \\textit{Classification and Regression by randomForest}, 2002, R News}.",
             paste("Table", table.count<<-table.count+1,"shows the confusion matrix of random forest."),
             paste("Figure", fig.count<<-fig.count+1,"shows the cumulative error rates of random forest analysis for given parameters.\n"),
             paste("Figure", fig.count<<-fig.count+1,"shows the important features ranked by random forest.\n"));

  cat(descr, file=rnwFile, append=TRUE);

  cmdhist<-c("\\begin{figure}[htp]",
             "\\begin{center}",
             paste("\\includegraphics[width=1.0\\textwidth]{", mbSetObj$imgSet$rf.cls,"}", sep=""),
             "\\caption{Cumulative error rates by Random Forest classification. The overall error rate is shown
             as the black line; the red and green lines represent the error rates for each class.}",
             "\\end{center}",
             paste("\\label{",mbSetObj$imgSet$rf.cls,"}", sep=""),
             "\\end{figure}");

  cat(cmdhist, file=rnwFile, append=TRUE, sep="\n");

  cmdhist<-c("<<echo=false, results=tex>>=",
             "GetRFConf.Table(mbSet)",
             "@",
             paste("The OOB error is ", GetRFOOB(mbSetObj)));
        
  cat(cmdhist, file=rnwFile, append=TRUE, sep="\n");

  cmdhist<-c( "\n\n",
              "\\begin{figure}[htp]",
              "\\begin{center}",
              paste("\\includegraphics[width=1.0\\textwidth]{", mbSetObj$imgSet$rf.imp,"}", sep=""),
              "\\caption{Significant features identified by Random Forest. The features are ranked by the mean
              decrease in classification accuracy when they are permuted.}",
              "\\end{center}",
              paste("\\label{",mbSetObj$imgSet$rf.imp,"}", sep=""),
              "\\end{figure}");
        
  cat(cmdhist, file=rnwFile, append=TRUE, sep="\n");
  cat("\\clearpage", file=rnwFile, append=TRUE);
}

# Functional Prediction
CreateFUNCPREDdoc <- function(mbSetObj){
  
  mbSetObj <- .get.mbSetObj(mbSetObj);
  
  # need to check if this process is executed
  if(is.null(mbSetObj$analSet$func.pred)){
    return();
  }
  
  descr <- c("\\subsection{Functional Prediction}\n",
             "16S rRNA sequencing data can also be utilized to predict the functional profiles",
             "of the corresponding microorganisms based on their phylogenetic distances or sequence",
             "similarities to those organisms whose complete genomes have been sequenced and annotated.",
             "Two most commonly available methods to perform such functional profiling are PICRUSt and ",
             "Tax4Fun. Along with the difference in the methodologies for prediction, the underlying database",
             "used for OTU annotation determines the usage of either of these approaches. PICRUSt is designed",
             "to work with Greengenes annotated features (OTUs) while Tax4Fun is based on SILVA annotated features",
             "for functional profiling. Due to such inconsistencies, both these approaches are limited and",
             "independent in their usage. MicrobiomeAnalyst allow flexibility in the annotation of OTUs, by",
             "asupporting both these approach.",
             "\n\n",
             "The result of both the approaches is a table of gene family abundances for each sample.",
             "Gene families are annotated using KEGG database (KO). Further, this predicted metagenomes",
             "can be organized into higher-level functional categories and also mapped to metabolic",
             "pathways using Shotgun Data Profiling (SDP) module in MicrobiomeAnalyst.",
             paste("Figure", fig.count<<-fig.count+1,"represents the boxplot for count distribution profile of predicted metagenomes."));
  
  cat(descr, file=rnwFile, append=TRUE);
  
  cmdhist<-c(
    "\\begin{figure}[htp]",
    "\\begin{center}",
    paste("\\includegraphics[width=1.0\\textwidth]{", mbSetObj$imgSet$func.pred,"}", sep=""),
    "\\caption{The count distribution (log-scaled) profile of metagenomes (KO) across each sample predicted using","\\texttt{", mbSetObj$analSet$func.meth, "}.}",
    "\\end{center}",
    paste("\\label{",mbSetObj$imgSet$func.pred,"}", sep=""),
    "\\end{figure}"
  );
  
  cat(cmdhist, file=rnwFile, append=TRUE, sep="\n");
  cat("\\clearpage", file=rnwFile, append=TRUE);
}

###############################################
## Shotgun Data profiling
##########################################

# Functional Profiling
CreateFUNCPROFdoc <- function(mbSetObj){
  
  mbSetObj <- .get.mbSetObj(mbSetObj);
  
  # need to check if this process is executed
  if(is.null(mbSetObj$analSet$func.prof)){
    return();
  }
  
  descr <- c("\\subsection{Functional Profiling}\n",
             "This method allows the mapping of gene abundance profile (output from shotgun metagenomics or predicted",
             "from PICRUSt/Tax4Fun) to functional or metabolic pathways for better biological understanding and interpretation.",
             "MicrobiomeAnalyst supports two high-quality resources KEGG and COG for providing organism-independent higher-level", 
             "annotation for genes. Such that functional profiling can be done at four higher-level functional categories including",
             "KEGG metabolisms, KEGG pathways, KEGG modules and COG functional categories.\n",
             "The most important parameter while performing functional profiling is how to compute abundance of higher functional",
             "categories since one gene can be assigned to multiple higher-level functional groups. MicrobiomeAnalyst offers three",
             "approaches to deal with this one to many mapping issue including Total hits, Total hits normalized by category size",
             "and Sum of the weighted hits. In Total hits method, the abundance of each category is calculated by summing up all",
             "the hits belong to that category. While, this calculated category abundance is further divided by category size in",
             "case of Total hits normalized by category size approach. Lastly, in Sum of the weighted hits method the weight of",
             "hits distributed equally to those categories which they belong.).",
             "\n\n",
             "The resultant higher-level functional profiles can be visualized as a stacked area plot across all the",
             "samples. The samples are organized by experimental factors to help visualize patterns of variations",
             "across different conditions.",
             paste("Figure", fig.count<<-fig.count+1,"represents the result of functional profiling as a stacked area plot."));
  
  cat(descr, file=rnwFile, append=TRUE,sep="\n");
  
  cmdhist<-c(
    "\\begin{figure}[htp]",
    "\\begin{center}",
    paste("\\includegraphics[width=1.0\\textwidth]{", mbSetObj$imgSet$func.prof,"}", sep=""),
    "\\caption{Functional profiling of gene abundance at ","\\texttt{", mbSetObj$analSet$func.lvl, "} level.}",
    "\\end{center}",
    paste("\\label{",mbSetObj$imgSet$func.prof,"}", sep=""),
    "\\end{figure}"
  );
  
  cat(cmdhist, file=rnwFile, append=TRUE, sep="\n");
  cat("\\clearpage", file=rnwFile, append=TRUE);
}

# create PCA doc
CreatePCAdoc <- function(mbSetObj){
  
  mbSetObj <- .get.mbSetObj(mbSetObj);
  
  # need to check if this process is executed
  if(is.null(mbSetObj$analSet$pca)){
    return();
  }
  
  descr <- c("\\subsection{Principal Component Analysis (PCA)}\n",
             "PCA is an unsupervised method aiming to find the directions that best",
             "explain the variance in a data set (X) without referring to class labels (Y).",
             "The data are summarized into much fewer variables called \\textit{scores} which",
             "are weighted average of the original variables. The weighting profiles are called",
             "\\textit{loadings}. The PCA analysis is performed using the \\texttt{prcomp} package.",
             "The calculation is based on singular value decomposition.\n",
             paste("Figure", fig.count<<-fig.count+1,"shows the 2-D scores plot between selected PCs\n")
             );
  
  cat(descr, file=rnwFile, append=TRUE);
  cmdhist<-c(
    "\\begin{figure}[htp]",
    "\\begin{center}",
    paste("\\includegraphics[width=1.0\\textwidth]{", mbSetObj$imgSet$pca,"}", sep=""),
    "\\caption{Scores plot between the selected PCs. The explained variances are shown in brackets.}",
    "\\end{center}",
    paste("\\label{",mbSetObj$imgSet$pca,"}", sep=""),
    "\\end{figure}"
  );
  cat(cmdhist, file=rnwFile, append=TRUE, sep="\n");
  cat("\\clearpage", file=rnwFile, append=TRUE, sep="\n");
}

###############################################
## Taxon Set enrichment analysis report
##########################################

# write .Rnw file template
CreateTaxaEnrichRnwReport<-function(mbSetObj, usrName){
  CreateHeader(usrName);
  CreateEnrichIntr();
  CreateEnrichOverview();
  CreateEnrichInputDoc(mbSetObj);
  CreateEnrichProcessDoc();
  CreateEnrichORAdoc();
  CreateFooter(mbSetObj);
}

CreateEnrichIntr<-function(){
  descr <- c("\\section{Background}\n",
             "Metagenome-wide association studies have already found strong association of microbes with host health and disease.", 
             "It also identifies a large number of microbes differentially regulated in various conditions.",
             "However, computational methods for analyzing such differentially regulated microbes from microbiome study are limited.",
             "TSEA or Taxon Set Enrichment Analysis is a way to identify biologically or ecologically meaningful",
             "patterns by analyzing them with context to pre-defined taxon set",
             "(microbes sharing some common trait) from a given list of significant features or microbes.",
             "approaches. In conventional approach, microbes are evaluated individually for their significance under conditions",
             "of study. Those microbes that have passed certain significance level are then combined to",
             "see if any meaningful patterns can be discerned. In contrast, TSEA directly investigates if",
             "a set of functionally related microbes without the need to preselect compounds based on",
             "some arbitary cut-off threshold. It has the potential to identify subtle but consistent changes",
             "among a group of related microbes, which may go undetected with the conventional approaches.",
             "\n\n",
             "Essentially, TSEA is a microbiome version of the popular GSEA (Gene Set Enrichment Analysis)",
             "software with its own collection of taxon set libraries as well as an implementation of",
             "user-friendly web-interfaces. GSEA is widely used in genomics data analysis and has proven to",
             "be a powerful alternative to conventional approaches. For more information, please refer to",
             "the original paper by Subramanian A, and a nice review paper by Nam D, Kim SY \\footnote{Dougu Nam Seon-Young Kim",
             "\\textit{Gene-set approach for expression pattern analysis}, Breif. in Bioinformatics 2008}.\n"
  );
  cat(descr, file=rnwFile, append=TRUE);
}

CreateEnrichOverview<-function(){
  descr <- c("\\section{TSEA Overview}\n",
             "Taxon set enrichment analysis consists of 4 steps - data input, data processing,",
             "data analysis, and results download. Based on the taxonomic resolution of microbes",
             "to analyze, three types of taxon sets are created and supported in MicrobiomeAnalyst",
             "Different taxon sets are selected based on different input types.",
             "Users can also perform taxon name mapping to higher taxonomic level and between a variety",
             "of microbes names and major database identifiers."
  );
  cat(descr, file=rnwFile, append=TRUE);
}

# create data input doc
CreateEnrichInputDoc<-function(mbSetObj){
  
  mbSetObj <- .get.mbSetObj(mbSetObj);
  
  descr <- c("\\section{Data Input}\n",
             "There are three enrichment analysis algorithms offered by TSEA. Accordingly, three",
             "different types of data inputs are required by these three approaches:\n"
  );
  
  cat(descr, file=rnwFile, append=TRUE);
  descr <- c(
    "\\begin{itemize}",
    "\\item{A list of microbes names characterized at any possible taxonomic level - entered as a one column data",
    "(\\textit{Mixed-level taxa});}",
    "\\item{A list of microbes names characterized at any species level - entered as a one column data",
    "(\\textit{Species-level taxa});}",
    "\\item{A list of microbes names (Binomial Nomenclature Name/GOLD ID/NCBI Taxonomy ID) characterized at
            any strain level - entered as a one column data",
    "(\\textit{Strain-level taxa)})}",
    "\\end{itemize}",
    "\n\n");
  
  cat(descr, file=rnwFile, append=TRUE, sep="\n");
  
  descr <- c("\\section{Selection of Taxon Set Library}\n",
             "Depending upon type of input list, Taxon set library will be selected.",
             "There are four built-in libraries offered by TSEA:\n"
  );
  cat(descr, file=rnwFile, append=TRUE, sep="\n");
  
  descr <- c(
    "\\begin{itemize}",
    "\\item{Mixed-level Taxon sets associated with human genetic variations - used with mixed-level taxa",
    "(\\textit{currently contains 1520 entries});}",
    "\\item{Mixed-level taxon sets associated with human diseases - used with mixed-level taxa",
    "(\\textit{currently contains 39 entries});}",
    "\\item{Species-level Taxon sets associated with human physiology, development, life styles etc. - used with species-level taxa)",
    "(\\textit{currently contains 170 entries})}",
    "\\item{Strain-level Taxon sets based on their shared phenotypic traits or ecological niches - used with strain-level taxa",
    "(\\textit{currently contains 100 entries})}",
    "\\end{itemize}",
    "\n\n");
  cat(descr, file=rnwFile, append=TRUE, sep="\n");
  
  if(mbSetObj$dataSet$q.type == "speciestaxa"){
    descr <- c("You have provided a list of microbes annotated at species-level. Species-level",
               "taxon sets will be used for performing enrichment analysis. \n\n");
  }else if(mbSetObj$dataSet$q.type == "straintaxa"){
    descr <- c("You have provided a list of microbes annotated at strain-level. Strain-level",
               "taxon sets will be used for performing enrichment analysis. \n\n");
  }else{
    if(mbSetObj$dataSet$tset.type=="dis") {
        descr <- c("You have provided a list of microbes annotated at mixed-level. Mixed-level",
                   "taxon sets associated with Human diseases will be used for performing enrichment analysis. \n\n");
    } else {
         descr <- c("You have provided a list of microbes annotated at mixed-level. Mixed-level",
                   "taxon sets associated with Human genetic variations will be used for performing enrichment analysis. \n\n");
    }
  }
  
  cat(descr, file=rnwFile, append=TRUE, sep="\n");
}

CreateEnrichProcessDoc<-function(){
  descr <- c("\\section{Data Processing}\n",
             "The first step is to match the user's entered microbes name with the microbes contained in the Taxon set library.",
             "All the microbes name should be in universally accepted format (Binomial Nomenclature).",
             "TSEA also provides conversion between microbe names and other database identifiers such",
             "as NCBI Taxonomy ID and GOLD ID. The unmapped name will be indicated by a \\textit{- or empty cell} and",
             "will be removed from further analysis.",
             "\n\n");
  cat(descr, file=rnwFile, append=TRUE, sep="\n");
  
  descr<-c("<<echo=false, results=tex>>=",
           "GetMapTable(mbSet)",
           "@",
           "\\clearpage\n\n");
  cat(descr, file=rnwFile, append=TRUE, sep="\n");
}


CreateEnrichORAdoc<-function(){

  descr <- c("\\section{Enrichment Analysis}\n",
    "Over Representation Analysis (ORA) is performed when",
    "a list of taxa or microbes is provided. The list of microbes can be obtained through",
    "differential abundance testing, or from biomarker analysis or from a clustering algorithm",
    "performed using MDP",
    "to investigate if some biologically meaningful patterns can be identified.",
    "\n\n",
    "ORA was implemented using the \\textit{hypergeometric test} to evaluate whether a particular",
    "Taxon set is represented more than expected by chance within the given compound list.",
    "One-tailed p values are provided after adjusting for multiple testing. \\textbf{Table 2} below",
    "provides the detail about enriched taxon set.");

  cat(descr, file=rnwFile, append=TRUE);
  cat("\n\n", file=rnwFile, append=TRUE); 
  descr<-c("<<echo=false, results=tex>>=",
           "GetORATable(mbSet)",
           "@",
           "\\clearpage\n\n");
  cat(descr, file=rnwFile, append=TRUE, sep="\n");
}

CreateFooter<-function(mbSetObj){
  
  mbSetObj <- .get.mbSetObj(mbSetObj);
  
  if(mbSetObj$module.type == "sdp"){
    descr <- c("\\section{Other Features}\n",
             "Please be advised that association analysis",
             "with its corresponding results are not included in this report.",
             "\n\n");
    cat(descr, file=rnwFile, append=TRUE); 
  }
  
  end <- c("\\vspace{5 mm}\n--------------------------------\n\n",
           "The report was generated on \\Sexpr{date()} with \\Sexpr{print(version$version.string)}.\n",
           "\\end{document}\n\n");
  cat(end, file=rnwFile, append=TRUE);
}


###############################################
## Projection with Public Data
##########################################

# write .Rnw file template
CreatePPDRnwReport<-function(mbSetObj, usrName){
    CreateHeader(usrName);
    CreatePPDIntr();
    CreatePPDOverview();
    CreateIOdoc();
    CreatePPDAnalDoc();
    CreatePPDResultDoc();
    CreateFooter(mbSetObj);
}

CreatePPDIntr<-function(){
        descr <- c("\\section{Background}\n",
                   "The 16S rRNA marker gene survey remain one of the most widely used method for microbiome",
                   "datasets, with large number of datasets publicly available. In principle, these studies",
                   "can be compared, as they based on the same bacterial 16S ribosomal RNA (rRNA) gene target.",
                   "However, some technical differences may arise when comparing two such different studies.",
                   "These technical differences include experimental protocols between laboratories; the way",
                   "samples are collected and stored, DNA extraction methods, PCR primers used for amplicons",
                   "generation, the targeted hyper-variable region of 16S rRNA, and the sequencing instruments",
                   "and technologies used. Despite these concerns, several cross-study comparisons of human",
                   "microbiota have revealed that biological differences outweigh the variations produced due",
                   "to these technical factors.\\footnote{Lozupone CA. \\textit{Meta-analyses of studies",
                   "of the human microbiota}, Genome Res. 2013}",
                   "\n\n",
                   "MicrobiomeAnalyst now allows researchers to perform a meta-analysis to reveal larger pictures",
                   "or novel insights beyond a single study through Projection with Public Data (PPD) module.",
                   "Users can perform visual analytics on their 16S rRNA data within the background of",
                   "compatible public datasets. Such meta-analysis will improve statistical power or biological understanding.\n");
        cat(descr, file=rnwFile, append=TRUE);
}

CreatePPDOverview<-function(){
  descr <- c("\\section{PPD Overview}\n",
             "Projection with Public Data consists of 4 main steps - data input and processing,",
             "public or reference dataset selection, data analysis and result visualization.");
  cat(descr, file=rnwFile, append=TRUE);
}

CreatePPDAnalDoc<-function(mbSetObj){
  
  mbSetObj <- .get.mbSetObj(mbSetObj);
  
  descr <- c("\\section{Selection of Reference Dataset}\n",
             "MicrobiomeAnalyst supports large-scale well-annotated 16S public datasets collected from Qiita.",
             "These datasets are internally processed using QIIME pipeline and used Greengenes identifier for", 
             "OTUs annotation.",
             "In this step, users are asked to select a public dataset for comparison and meta analysis with",
             "their own data. Currently, user can choose from 34",
             "public datasets for performing such meta-analysis in MicrobiomeAnalyst.\n");
        
  cat(descr, file=rnwFile, append=TRUE, sep="\n");

  descr <- c("\\subsection{Public dataset Library}\n",
             "The datasets are organized into 4 main categories by biome of interest, with a total of 34 datsets :\n",
             "\\begin{itemize}",
             "\\item{Human",
             "(\\textit{currently contains 24 datasets from 4 body sites (Gut, Skin, Oral, Vagina)});}",
             "\\item{Cow",
             "(\\textit{currently contains 2 datasets});}",
             "\\item{Rat",
             "(\\textit{currently contains 3 datasets});}",
             "\\item{Environment",
             "(\\textit{currently contains 5 datasets})}",
             "\\end{itemize}",
             "\n\n",
             mbSetObj$dataSet$lib.msg,
             "\n");
        
  cat(descr, file=rnwFile, append=TRUE, sep="\n");

  descr <- c("\\subsection{Meta analysis}\n",
             "User data is merged with selected reference dataset based on common OTUs or feature names.An internal mapping is performed for OTUs",
             "annotated other than a Greengenes identifier. At least 20 percent overlapping of OTUs is required between the user and reference data to have",
             "a meaningful comparison between the datasets. Further, samples with shallow sequencing depth (library size < 100) have been removed",
             "from further processing. Afterwards, all the samples are rarefied at even sequencing depth to bring them at the same scale", 
             "for assessment.",
             "\n\n",
             "Beta-diversity is mostly affected by those abundant taxa that are shared among the samples.",
             "Normalization tends to have minimal effect on data clustering patterns. These considerations",
             "have been adapted in MicrobiomeAnalyst to save computing time, in which the default PCoA is based on",
             "the Bray-Curtis index distance computed from the top 20% most abundant taxa. However,",
             "the user can also perform such analysis using complete dataset and several different distance",
             "methods to discover clustering patterns.",
             "\n");
  cat(descr, file=rnwFile, append=TRUE, sep="\n");
}

CreatePPDResultDoc<-function(mbSetObj){
  
  mbSetObj <- .get.mbSetObj(mbSetObj);
  
  descr <- c("\\section{Meta Analysis Result}\n",
             "The results are represented in an interactive 3D (website) and 2D (report) PCoA plot with node colors based on different",
             "experimental factors and node shapes indicating different studies or datasets. Users can",
             "intuitively rotate (mouse dragging), zoom in and out (mouse scrolling), or directly click",
             "on any node or sample of interest to view its composition at any taxonomic level.\n",
             "Community-level differences or similarities between the user and reference data is interrogated by calculating the amount of diversity",
             "that is shared between samples (beta-diversity), followed by clustering using unsupervised multivariate statistical methods such as",
             "Principal Coordinate Analysis (PCoA). These techniques often provide clear associations between subject characteristics and", 
             "overall diversity.",
             paste("Figure", fig.count<<-fig.count+1,"shows the 2-D PCoA plot where user uploaded samples are represented as triangular nodes."),
             "\n\n",
             mbSetObj$analSet$topo.msg,
             "\n");
        
  cat(descr, file=rnwFile, append=TRUE, sep="\n");
        
  fig <- c("\\begin{figure}[htp]",
           "\\begin{center}",
           paste("\\includegraphics[width=1.0\\textwidth]{",mbSetObj$imgSet$ppd.2d,"}",sep=""),
           "\\caption{2-D PCoA plot}",
           "\\end{center}",
           paste("\\label{",mbSetObj$imgSet$ppd.2d,"}", sep=""),
           "\\end{figure}",
           "\\clearpage\n\n");
        
  cat(fig, file=rnwFile, append=TRUE, sep="\n");
}
