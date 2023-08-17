##################################################
## R script for MicrobiomeAnalyst
## Description: predict KO abundance from 16S data (picrust(gg_12 and gg_13), tax4fun 1.0 and 2.0)
## Authors: Jeff Xia, jeff.xia@mcgill.ca; Yao Lu, yao.lu5@mail.mcgill.ca
###################################################

my.16sfun.anot<-function(mbSetObj, type, pipeline,ggversion) {
   print(c(type,pipeline,ggversion))
  mbSetObj <- .get.mbSetObj(mbSetObj);
  mbSetObj$dataSet$type <- type;
  merge.otu <- qs::qread("data.orig");
 
  merge.otu <- apply(merge.otu, 2, as.numeric);
  current.msg<<- "null"
  err.vec <<- "null"
 if(pipeline %in% c("Ref99NR","Ref100NR")){
    func.meth<-"Tax4Fun2";
    database_mode <- pipeline
    path_to_reference_data <- get.fun.lib.path(func.meth);
    
    if(nrow(merge.otu)<2000){
      seqInput<- as.list(rownames(mbSetObj$dataSet$data.orig))
    }else{
      merge.otu <- qs::qread("filt.data.orig")
     ASV_ID_mapping <- data.frame(mbSetObj$dataSet$master.nms)
      seqInput<- as.list(ASV_ID_mapping$OriginalID[match(rownames(merge.otu),ASV_ID_mapping$NewID)])
    }
    
  
    seqinr::write.fasta(seqInput,names=1:length(seqInput), nbchar=80,file.out = "user_otu_input.fasta")

    user_otu_table<- data.frame(cbind(ID=1:nrow(merge.otu),merge.otu),check.names = F)
    write.csv(user_otu_table,"user_otu_table.csv",row.names = F)
    runRefBlast(path_to_otus="user_otu_input.fasta", 
                               path_to_reference_data, 
                               path_to_temp_folder = "tax4fun2_prediction", 
                               database_mode, 
                               use_force = T, num_threads = 2, 
                               include_user_data = F, path_to_user_data = '', 
                               name_of_user_data = "")   

    makeFunctionalPrediction(path_to_otu_table='user_otu_table.csv',
           path_to_reference_data, 
           path_to_temp_folder = "tax4fun2_prediction", database_mode, 
           normalize_by_copy_number = T, min_identity_to_reference = 97,
           include_user_data = F, path_to_user_data = '', 
           name_of_user_data = '', use_uproc = T, normalize_pathways = F)  
              
 
  
    result <- read.csv("tax4fun2_prediction/KOpred.csv",check.names = F,row.names=1)
    result_path<- read.csv("tax4fun2_prediction/pathpred.csv",check.names = F,row.names=1)
    result<-as.data.frame(round(1000000*result),check.names=FALSE);  # get to integers
    result_path<-as.data.frame(round(1000000*result_path),check.names=FALSE);  # get to integers

}else{

    #getting whole taxa labels back
    if(type=="SILVA"){
     rownames(merge.otu) <- unname(mbSetObj$dataSet$comp_taxnm);
        func.meth<-"Tax4Fun";
        folderReferenceData <- get.fun.lib.path(func.meth);

        #load_tax4fun();
        if(pipeline=="qi_silva"){
            ModSilvaIds <- gsub("uncultured archaeon","",rownames(merge.otu));
            ModSilvaIds <- gsub("uncultured organism","",ModSilvaIds);
            ModSilvaIds <- gsub("uncultured bacterium","",ModSilvaIds);
            ModSilvaIds <- gsub("uncultured crenarchaeote","",ModSilvaIds);
            ModSilvaIds <- gsub("uncultured euryarchaeote","",ModSilvaIds);
            ModSilvaIds <- gsub("; ",";",ModSilvaIds);
            rownames(merge.otu)<-ModSilvaIds;
            merge.otu <- rowsum(as.data.frame(merge.otu),ModSilvaIds);
        }
        rownames(merge.otu) <- gsub(";  ",";",rownames(merge.otu));
        rownames(merge.otu) <- gsub("; ",";",rownames(merge.otu));
        rownames(merge.otu) <- gsub(" ","_",rownames(merge.otu))
        idx=which(!(grepl(";$",rownames(merge.otu))))
        rownames(merge.otu)[idx] <- paste0(rownames(merge.otu)[idx],";");
        rownames(merge.otu) <- gsub("Bacteroidota","Bacteroidetes",rownames(merge.otu));
        rownames(merge.otu) <- gsub("Enterobacterales","Enterobacteriales",rownames(merge.otu));
        rownames(merge.otu) <- stringr::str_trim(rownames(merge.otu) ,side="both")
        data<-list(sampleNames=colnames(merge.otu),otuTable=merge.otu);

         Tax4FunOutput <- Tax4Fun(data, folderReferenceData, fctProfiling = TRUE, refProfile = "UProC", shortReadMode = TRUE, normCopyNo = TRUE);
      
         if(length(Tax4FunOutput)==1){
          return(0)
         }
        result2<-as.data.frame(Tax4FunOutput[1],check.names=FALSE);

        #removing unnecessary data from column name.
        colnames(result2)<-sub("Tax4FunProfile.","",colnames(result2));
        colnames(result2)<-substr(colnames(result2),1,6);
        result<-t(result2);
        result<-as.data.frame(round(1000000*result),check.names=FALSE);  # get to integers
    } else {
      if(ggversion=="gg13"){
        func.meth<-"picrust13";
        picrust_path <- get.fun.lib.path(func.meth);
        if(type == "Greengenes"){
            a<-rownames(merge.otu);
            rownames(merge.otu)<-NULL;
            merge.otu<-data.frame(merge.otu,check.names=FALSE);
            merge.otu <- cbind(a,merge.otu)
            colnames(merge.otu)[1]<-"X.OTU_IDs";
      
            otu.dic <<- qs::qread(paste0(picrust_path, "/greengenes_tax.qs"));

            merge.otu[,1]<-as.character(merge.otu[,1]);
            merge.otu[,1]<-otu.dic[match(merge.otu$X.OTU_IDs,otu.dic$Greengenes),1];
            merge.otu<-merge.otu[!is.na(merge.otu$X.OTU_IDs),];
            nm<-rownames(merge.otu)<-merge.otu[,1];
            merge.otu<-merge.otu[,-1];
            merge.otu<-apply(merge.otu, 2, as.numeric);
            rownames(merge.otu)<-nm;
        }
    
        query<-as.data.frame(merge.otu,check.names=FALSE);
        samplenm<-colnames(query);
    
        #normalizing by 16S copy number
        copyno <- qs::qread(paste0(picrust_path, "/16S_copyno.qs"));
        result2<-merge(query,copyno, by ="row.names");
        index1<-match(samplenm,colnames(result2), nomatch = NA_integer_, incomparables = NULL);
        result2[index1]<-result2[index1]/result2[['X16S_rRNA_Count']];
        result2[index1]<-round(result2[index1],2);
        rownames(result2)<-result2[,1];
        result2<-result2[,-1];
    
        # update picrust database (gg_13_5), 
        #need to fetch and merge results from 14 parts of picrust to get around memory issue (from 10.8G => 780M)
        res <- data.frame(matrix(0, nrow=nrow(result2), ncol= 6909),check.names=FALSE);
        row.names(res) <- row.names(result2);
    
        pc.lib <- qs::qread(paste0(picrust_path, "/picrust_part1.qs"));
    
        row.hits <- match(row.names(result2), rownames(pc.lib));
        res[, 1:500] <- pc.lib[row.hits,];
        colnames(res)[1:500] <- colnames(pc.lib);
    
        pc.lib <- qs::qread(paste0(picrust_path, "/picrust_part2.qs"));
    
        res[, 501:1000] <- pc.lib[row.hits,];
        colnames(res)[501:1000] <- colnames(pc.lib);
     
        pc.lib <- qs::qread(paste0(picrust_path, "/picrust_part3.qs"));
    
        res[, 1001:1500] <- pc.lib[row.hits,];
        colnames(res)[1001:1500] <- colnames(pc.lib);
    
        pc.lib <- qs::qread(paste0(picrust_path, "/picrust_part4.qs"));
    
        res[, 1501:2000] <- pc.lib[row.hits,];
        colnames(res)[1501:2000] <- colnames(pc.lib);
    
        pc.lib <- qs::qread(paste0(picrust_path, "/picrust_part5.qs"));
    
        res[, 2001:2500] <- pc.lib[row.hits,];
        colnames(res)[2001:2500] <- colnames(pc.lib);

        pc.lib <- qs::qread(paste0(picrust_path, "/picrust_part6.qs"));
    
        row.hits <- match(row.names(result2), rownames(pc.lib));
        res[, 2501:3000] <- pc.lib[row.hits,];
        colnames(res)[2501:3000] <- colnames(pc.lib);
    
        pc.lib <- qs::qread(paste0(picrust_path, "/picrust_part7.qs"));
    
        res[, 3001:3500] <- pc.lib[row.hits,];
        colnames(res)[3001:3500] <- colnames(pc.lib);
    
        pc.lib <- qs::qread(paste0(picrust_path, "/picrust_part8.qs"));
    
        res[, 3501:4000] <- pc.lib[row.hits,];
        colnames(res)[3501:4000] <- colnames(pc.lib);
    
        pc.lib <- qs::qread(paste0(picrust_path, "/picrust_part9.qs"));
    
        res[, 4001:4500] <- pc.lib[row.hits,];
        colnames(res)[4001:4500] <- colnames(pc.lib);
    
        pc.lib <- qs::qread(paste0(picrust_path, "/picrust_part10.qs"));
    
        res[, 4501:5000] <- pc.lib[row.hits,];
        colnames(res)[4501:5000] <- colnames(pc.lib);

        pc.lib <- qs::qread(paste0(picrust_path, "/picrust_part11.qs"));
    
        row.hits <- match(row.names(result2), rownames(pc.lib));
        res[, 5001:5500] <- pc.lib[row.hits,];
        colnames(res)[5001:5500] <- colnames(pc.lib);
    
        pc.lib <- qs::qread(paste0(picrust_path, "/picrust_part12.qs"));
    
        res[, 5501:6000] <- pc.lib[row.hits,];
        colnames(res)[5501:6000] <- colnames(pc.lib);
    
        pc.lib <- qs::qread(paste0(picrust_path, "/picrust_part13.qs"));
    
        res[, 6001:6500] <- pc.lib[row.hits,];
        colnames(res)[6001:6500] <- colnames(pc.lib);
    
        pc.lib <- qs::qread(paste0(picrust_path, "/picrust_part14.qs"));
    
        res[, 6501:6910] <- pc.lib[row.hits,];
        colnames(res)[6501:6909] <- colnames(pc.lib);
    
       
        ko.nms <- colnames(res);
        res<-merge(result2,res, by ="row.names");
    
      }else if(ggversion=="gg12"){
        func.meth<-"picrust12";
        picrust_path <- get.fun.lib.path(func.meth);
        if(type == "Greengenes"){
            a<-rownames(merge.otu);
            rownames(merge.otu)<-NULL;
            merge.otu<-data.frame(merge.otu,check.names=FALSE);
            merge.otu <- cbind(a,merge.otu)
            colnames(merge.otu)[1]<-"X.OTU_IDs";
      
            otu.dic <<- qs::qread(paste0(picrust_path, "/greengenes_taxmap.qs"));

            merge.otu[,1]<-as.character(merge.otu[,1]);
            merge.otu[,1]<-otu.dic[match(merge.otu$X.OTU_IDs,otu.dic$Greengenes),1];
            merge.otu<-merge.otu[!is.na(merge.otu$X.OTU_IDs),];
            nm<-rownames(merge.otu)<-merge.otu[,1];
            merge.otu<-merge.otu[,-1];
            merge.otu<-apply(merge.otu, 2, as.numeric);
            rownames(merge.otu)<-nm;
        }
    
        query<-as.data.frame(merge.otu,check.names=FALSE);
        samplenm<-colnames(query);
    
        #normalizing by 16S copy number
        copyno <- qs::qread(paste0(picrust_path, "/16S_copyno.qs"));
        result2<-merge(query,copyno, by ="row.names");
        index1<-match(samplenm,colnames(result2), nomatch = NA_integer_, incomparables = NULL);
        result2[index1]<-result2[index1]/result2[['X16S_rRNA_Count']];
        result2[index1]<-round(result2[index1],2);
        rownames(result2)<-result2[,1];
        result2<-result2[,-1];
    
        # need to fetch and merge results from 5 parts of picrust to get around memory issue (from 2G => 400M)
        res <- data.frame(matrix(0, nrow=nrow(result2), ncol= 6885),check.names=FALSE);
        row.names(res) <- row.names(result2);
    
        pc.lib <- qs::qread(paste0(picrust_path, "/picrust_part1.qs"));
    
        row.hits <- match(row.names(result2), rownames(pc.lib));
        res[, 1:1377] <- pc.lib[row.hits,];
        colnames(res)[1:1377] <- colnames(pc.lib);
    
        pc.lib <- qs::qread(paste0(picrust_path, "/picrust_part2.qs"));
    
        res[, 1378:2754] <- pc.lib[row.hits,];
        colnames(res)[1378:2754] <- colnames(pc.lib);
    
        pc.lib <- qs::qread(paste0(picrust_path, "/picrust_part3.qs"));
    
        res[, 2755:4131] <- pc.lib[row.hits,];
        colnames(res)[2755:4131] <- colnames(pc.lib);
    
        pc.lib <- qs::qread(paste0(picrust_path, "/picrust_part4.qs"));
    
        res[, 4132:5508] <- pc.lib[row.hits,];
        colnames(res)[4132:5508] <- colnames(pc.lib);
    
        pc.lib <- qs::qread(paste0(picrust_path, "/picrust_part5.qs"));
    
        res[, 5509:6885] <- pc.lib[row.hits,];
        colnames(res)[5509:6885] <- colnames(pc.lib);
    
        ko.nms <- colnames(res);
        res<-merge(result2,res, by ="row.names");
   
      }

        index2<-match(samplenm,colnames(res), nomatch = NA_integer_, incomparables = NULL);
        index3<-match(ko.nms,colnames(res), nomatch = NA_integer_, incomparables = NULL);
    
        myList <- vector('list', length(index2));
    
        for (i in 1:length(index2)) {
            myList[[i]]<-data.frame(colSums(res[index3]*res[,index2[i]]),check.names=FALSE);
        }
    
        res <- NULL;
        gc();
    
        MyMerge <- function(x, y){
            df<- merge(x, y, by= "row.names", all.x= F, all.y= F);
            rownames(df) <- df$Row.names
            df$Row.names <- NULL
            return(df)
        }
    
        result<- Reduce(MyMerge, myList);
        #orignal class label
        colnames(result)<-samplenm;
        
    }
  }

  result <- round(result, digits =0);
  if(func.meth == "Tax4Fun2"){
    kos <- matrix(rownames(result))
    colnames(kos) <- "#NAME"
    tax4fun2_feats <- cbind(kos, result)
    fast.write(tax4fun2_feats, file="functionalprof_tax4fun2.csv", row.names = FALSE);
   
    result_path <- round(result_path, digits =0);
    pathkegg <- matrix(rownames(result_path))
    colnames(pathkegg) <- "#PATHWAY"
    tax4fun2_path <- cbind(pathkegg, result_path)
    fast.write( tax4fun2_path, file="pathwayprof_tax4fun2.csv", row.names = FALSE);
   
    result_path <- result_path[!apply(result_path==0,1,all), ]; 
    mbSetObj$analSet$pathway.pred<- result_path;
    result <- result[!apply(result==0,1,all), ]; #filtering zero counts across all
    mbSetObj$analSet$func.pred.tax4fun2<-result
  }else if(func.meth == "Tax4Fun"){
    kos <- matrix(rownames(result))
    colnames(kos) <- "#NAME"
    tax4_feats <- cbind(kos, result)
    fast.write(tax4_feats, file="functionalprof_tax4fun.csv", row.names = FALSE);
    result <- result[!apply(result==0,1,all), ]; #filtering zero counts across all
     mbSetObj$analSet$func.pred.tax4fun<-result
  }else if(func.meth == "picrust12"){
    kos <- matrix(rownames(result))
    colnames(kos) <- "#NAME"
    picrust_feats <- cbind(kos, result)
    fast.write(picrust_feats, file="functionalprof_picrust_gg12.csv", row.names=FALSE);
     result <- result[!apply(result==0,1,all), ]; #filtering zero counts across all
    mbSetObj$analSet$func.pred.gg12<-result
   }else if(func.meth == "picrust13"){
    kos <- matrix(rownames(result))
    colnames(kos) <- "#NAME"
    picrust_feats <- cbind(kos, result)
    fast.write(picrust_feats, file="functionalprof_picrust_gg13.csv", row.names=FALSE);
    result <- result[!apply(result==0,1,all), ]; #filtering zero counts across all
    mbSetObj$analSet$func.pred.gg13<-result

   }

 
  
  # create meta-data file
  meta <- as.matrix(mbSetObj$dataSet$sample_data)
  samplenames <- matrix(rownames(meta))
  colnames(samplenames) <- "#NAME"
  meta_all <- cbind(samplenames, meta)  
  fast.write(meta_all, file="metadata.csv", row.names = FALSE)
  
  # save as RDS for memory saving
 ;
  mbSetObj$analSet$func.meth<-func.meth;
  
  return(.set.mbSetObj(mbSetObj));
  
}


##### adapted from tax4fun2 package  https://doi.org/10.1186/s40793-020-00358-7

# Function to reduce blast results
blastTableReducer <-  function(path_to_blast_file = '') {
  if(file.size(path_to_blast_file) == 0) stop('Blast file empty!')
  
  id1 = ""
  file_in = file(description = path_to_blast_file, open = "r")
  file_out = file(description = paste0(path_to_blast_file, ".tmp"), open = "w")
  while (TRUE){
    line = readLines(con = file_in, n = 1)
    if (length(line) == 0) break
    
    id2 = strsplit(x = line, split = "\t", fixed = T)[[1]][1]
    if(id1 != id2){
      id1 = id2
      write(x = line, file = file_out, append = T)
    }
  }
  close(file_in)
  close(file_out)
  
  # Remove old file and rename tmp file
  file.remove(path_to_blast_file)
  file.rename(from = paste0(path_to_blast_file, ".tmp"), to = path_to_blast_file)
}


runRefBlast <-  function(path_to_otus="user_otus.fasta", 
                         path_to_reference_data = "Tax4Fun2_ReferenceData", 
                         path_to_temp_folder = "tax4fun2_prediction", 
                         database_mode = 'Ref100NR', 
                         use_force = T, num_threads = 1, 
                         include_user_data = F, path_to_user_data = '', 
                         name_of_user_data = ""){
  blast_bin = file.path(path_to_reference_data, "blast_bin/blastn")
  if (tolower(Sys.info()[["sysname"]]) != "windows") system(paste("chmod +x", blast_bin))
  if (tolower(Sys.info()[["sysname"]]) == "windows") blast_bin = file.path(path_to_reference_data, "blast_bin/blastn.exe")
  res = system(command = paste(blast_bin, "-help"), intern = T)
  if(length(res) == 0)
  {
    blast_bin = "blastn"
    res = system(command = paste(blast_bin, "-help"), intern = T)
    if(length(res) == 0) stop("blastn not found! Consider to use the buildDependencies() command.")
  }
  
  makeblastdb_bin = file.path(path_to_reference_data, "blast_bin/makeblastdb")
  if (tolower(Sys.info()[["sysname"]]) != "windows") system(paste("chmod +x", makeblastdb_bin))
  if (tolower(Sys.info()[["sysname"]]) == "windows") makeblastdb_bin = file.path(path_to_reference_data, "blast_bin/makeblastdb.exe")
  res = system(command = paste(makeblastdb_bin, "-help"), intern = T)
  if(length(res) == 0){
    makeblastdb_bin = "makeblastdb"
    res = system(command = paste(makeblastdb_bin, "-help"), intern = T)
    if(length(res) == 0) stop("makeblastdb not found! Consider to use the buildDependencies() command.")
  }
  
  
  # Check for the presence of the otu file and the reference data
  if(!file.exists(path_to_otus)) stop('Otu file not found!')
  if(!dir.exists(path_to_reference_data)) stop('Reference data not found!')
  
  # Choose which refernence data set is used
  if(database_mode == 'Ref99NR') {
    path_to_ref_db = file.path(path_to_reference_data, 'Ref99NR/Ref99NR.fasta')
  } else if (database_mode == 'Ref100NR') {
    path_to_ref_db = file.path(path_to_reference_data, 'Ref100NR/Ref100NR.fasta')
  } else {
    stop('Database mode unknown! valid choces are Ref99NR, Ref100NR')
  }
  if(!file.exists(path_to_ref_db)) stop('Reference database not found!')
  
  # Check for the presence of the tmp folder. If present, stop or delete, if not, create

  if(dir.exists(path_to_temp_folder)){
    if(use_force) unlink(x = path_to_temp_folder, recursive = T)
    if(!use_force) stop('Temporay folder exist! Use \'use_force = T\' to overwrite')
  }
  dir.create(path_to_temp_folder)
  
  # Write to log file
  path_to_log_file = file.path(path_to_temp_folder, 'logfile1.txt')
  write(x = "RefBlast", file = path_to_log_file, append = F)
  write(x = database_mode, file = path_to_log_file, append = T)
  write(x = path_to_otus, file = path_to_log_file, append = T)
  if(include_user_data) write(x = 'User data will be included', file = path_to_log_file, append = T)
  write(x = date(), file = path_to_log_file, append = T)
  
  message('Copy and generate database')
  
  cmd = paste(makeblastdb_bin, '-dbtype nucl -in', path_to_ref_db)
  if (tolower(Sys.info()[["sysname"]]) == "windows") system(cmd, show.output.on.console = F)
  if (tolower(Sys.info()[["sysname"]]) != "windows") system(cmd, ignore.stdout = T, ignore.stderr = T)
  message('Reference blast started')
  cmd = paste(blast_bin, '-db', path_to_ref_db, '-query', path_to_otus, '-evalue 1e-20 -max_target_seqs 1000000 -outfmt 6 -out', file.path(path_to_temp_folder, 'ref_blast.txt'), '-num_threads', num_threads)
  if (tolower(Sys.info()[["sysname"]]) == "windows") system(cmd, show.output.on.console = F)
  if (tolower(Sys.info()[["sysname"]]) != "windows") system(cmd, ignore.stdout = T, ignore.stderr = T)
  blastTableReducer(file.path(path_to_temp_folder, 'ref_blast.txt'))
  message('Reference blast finished')
  message('Cleanup finished')
  
  if(include_user_data){
    message('Copy and generate user database')
    path_to_user_db = file.path(path_to_user_data, name_of_user_data, 'USER.fasta')
    file.copy(from = path_to_user_db, to = file.path(path_to_temp_folder, 'user_blast.fasta'))
    cmd = paste(makeblastdb_bin, '-dbtype nucl -in', file.path(path_to_temp_folder, 'user_blast.fasta'))
    if (tolower(Sys.info()[["sysname"]]) == "windows") system(cmd, show.output.on.console = F)
    if (tolower(Sys.info()[["sysname"]]) != "windows") system(cmd, ignore.stdout = T, ignore.stderr = T)
    message('User blast started')
    cmd = paste(blast_bin, '-db', file.path(path_to_temp_folder, '/user_blast.fasta'), '-query', path_to_otus, '-evalue 1e-20 -max_target_seqs 1000000 -outfmt 6 -out', file.path(path_to_temp_folder, 'user_blast.txt'), '-num_threads', num_threads)
    if (tolower(Sys.info()[["sysname"]]) == "windows") system(cmd, show.output.on.console = F)
    if (tolower(Sys.info()[["sysname"]]) != "windows") system(cmd, ignore.stdout = T, ignore.stderr = T)
    blastTableReducer(file.path(path_to_temp_folder, 'user_blast.txt'))
    message('User blast finished')
    unlink(file.path(path_to_temp_folder, 'user_blast.fasta*'))
    message('Cleanup finished')
  }
  
}



makeFunctionalPrediction <- function(path_to_otu_table='user_otus.fasta',
                                    path_to_reference_data = "Tax4Fun2_ReferenceData", 
                                    path_to_temp_folder = "tax4fun2_prediction", database_mode = 'Ref99NR', 
                                    normalize_by_copy_number = T, min_identity_to_reference = 97,
                                    include_user_data = F, path_to_user_data = '', 
                                    name_of_user_data = '', use_uproc = T, normalize_pathways = F){
  # Read old log file to control the database mode
  path_to_log_file = file.path(path_to_temp_folder, 'logfile1.txt')
  log_file = read.delim(path_to_log_file, header = F)
  if(database_mode != as.character(log_file$V1[2])) stop('Logfile indicates other database mode')
  if(as.character(log_file$V1[4]) == "User data will be included" & include_user_data == FALSE) warning("Your log file indicates that you used user data in the first step\nSet include_user_data to TRUE to incorporate!")
  if(as.character(log_file$V1[4]) != "User data will be included" & include_user_data == TRUE)
  {
    warning("Your log file indicates that you did not used user data in the first step\nSet include_user_data will be set to FALSE!")
    include_user_data = FALSE
  }
  
  # Set path to reference files
  if(database_mode == 'Ref99NR') {
    path_to_ref_profiles = file.path(path_to_reference_data, 'Ref99NR')
  } else if (database_mode == 'Ref100NR') {
    path_to_ref_profiles = file.path(path_to_reference_data, 'Ref100NR')
  } else {
    stop('Database mode unknown! valid choces are Ref99NR and Ref100NR')
  }
  
  if(min_identity_to_reference < 1) min_identity_to_reference = min_identity_to_reference * 100
  if(min_identity_to_reference < 90) warning("Minimum identity of less than 90% will likly results in inaccurate predictions!")
  message(paste0("Using minimum idenity cutoff of ", min_identity_to_reference, "% to nearest neighbor"))
  # Write to new log file
  path_to_log_file = file.path(path_to_temp_folder, 'logfile2.txt')
  write(x = "Tax4fun2 beta", file = path_to_log_file, append = F)
  write(x = date(), file = path_to_log_file, append = T)
  
  # Reading and reducing the blast file
  ref_blast_result = read.delim(file.path(path_to_temp_folder, 'ref_blast.txt'), h = F)
  ref_blast_result_reduced = ref_blast_result[which(ref_blast_result$V3 >= min_identity_to_reference), 1:2]
  
  if(include_user_data) {
    ref_blast_result_reduced_v1 = as.character(ref_blast_result_reduced$V1)
    ref_blast_result_reduced_v2 = as.character(ref_blast_result_reduced$V2)
    user_blast_result = read.delim(paste(path_to_temp_folder, '/user_blast.txt', sep = ''), h = F)
    user_blast_result_reduced = user_blast_result[which(user_blast_result$V3 >= min_identity_to_reference), 1:2]
    for(i in 1:nrow(user_blast_result_reduced)) {
      user_id = as.character(user_blast_result_reduced$V1)[i]
      if(user_id %in% ref_blast_result_reduced_v1){
        j = which(ref_blast_result_reduced$V1 == user_id)
        ref_blast_result_reduced_v2[j] = as.character(user_blast_result_reduced$V2)[i]
      } else {
        ref_blast_result_reduced_v1 = c(ref_blast_result_reduced_v1, user_id)
        ref_blast_result_reduced_v2 = c(ref_blast_result_reduced_v2, as.character(user_blast_result_reduced$V2)[i])
      }
    }
    ref_blast_result_reduced = data.frame(V1 = ref_blast_result_reduced_v1, V2 = ref_blast_result_reduced_v2)
  }
  #ref_blast_result_reduced
  
  # Reading and filtering the otu table
  otu_table = read.csv(path_to_otu_table,,check.names = F)
  otu_table_reduced = merge(x = ref_blast_result_reduced, y = otu_table, by.x = 'V1', by.y = names(otu_table)[1])[,-1]
  otu_table_reduced_aggregated = aggregate(x = otu_table_reduced[,-1], by = list(otu_table_reduced[,1]), sum)
  
  # Write unknown fraction to log file
  if((ncol(otu_table) - 1) == 1){
    unknown_fraction1 = as.data.frame(round(1 - sum(ifelse(otu_table_reduced[,-1]>0,1,0)) / sum(ifelse(otu_table[,-1]>0,1,0)), digits = 5))
    write(x = 'Unknown fraction (amount of otus unused in the prediction) for each sample:', file = path_to_log_file, append = T)
    write.table(x = unknown_fraction1, file = path_to_log_file, append = T, quote = F, sep = ': ', row.names = T, col.names = F)
    unknown_fraction2 = as.data.frame(round(1 - sum(otu_table_reduced_aggregated[,-1]) / sum(otu_table[,-1]), digits = 5))
    write(x = 'Unknown fraction (amount of sequences unused in the prediction) for each sample:', file = path_to_log_file, append = T)
    write.table(x = unknown_fraction2, file = path_to_log_file, append = T, quote = F, sep = ': ', row.names = T, col.names = F)
  } else {
    unknown_fraction1 = as.data.frame(round(1 - colSums(ifelse(otu_table_reduced[,-1]>0,1,0)) / colSums(ifelse(otu_table[,-1]>0,1,0)), digits = 5))
    write(x = 'Unknown fraction (amount of otus unused in the prediction) for each sample:', file = path_to_log_file, append = T)
    write.table(x = unknown_fraction1, file = path_to_log_file, append = T, quote = F, sep = ': ', row.names = T, col.names = F)
    unknown_fraction2 = as.data.frame(round(1 - colSums(otu_table_reduced_aggregated[,-1]) / colSums(otu_table[,-1]), digits = 5))
    write(x = 'Unknown fraction (amount of sequences unused in the prediction) for each sample:', file = path_to_log_file, append = T)
    write.table(x = unknown_fraction2, file = path_to_log_file, append = T, quote = F, sep = ': ', row.names = T, col.names = F)
  }
  
  # normalize or not normalize is the question
  n = 1
  if(use_uproc) n = 3
  if(normalize_by_copy_number) n = n + 1
  
  # Generate reference profile
  message('Generating reference profile')
  reference_profile = NULL
  for(reference_id in otu_table_reduced_aggregated$Group.1){
    #print(reference_id)
    if(grepl(pattern = 'user', x = reference_id)){
      path_to_profile = file.path(path_to_user_data, name_of_user_data)
      reference_id = gsub("_[0-9]*", "", reference_id)
      reference_file_path = file.path(path_to_profile, paste0(reference_id, '.tbl'))
    } else{
      path_to_profile = path_to_ref_profiles
      reference_file_path = file.path(path_to_profile, paste0(reference_id, '.tbl.gz'))
    }
    reference_file = read.delim(file = reference_file_path)
    reference_profile = rbind(reference_profile, as.numeric(reference_file[,n]))
  }
  dim(reference_profile)
  
  ko_list = read.delim(file.path(path_to_reference_data, 'KEGG/ko.txt'))
  
  # Calculate functional profiles sample-wise
  message('Generating functional profile for:')
  functional_prediction = NULL
  for(sample in 2:ncol(otu_table_reduced_aggregated)){
    message(names(otu_table_reduced_aggregated[sample]))
    functional_prediction_sample = reference_profile * as.numeric(otu_table_reduced_aggregated[,sample])
    functional_prediction_sample = colMeans(functional_prediction_sample)
    functional_prediction_sample = functional_prediction_sample / sum(functional_prediction_sample)
    if(is.na(sum(functional_prediction_sample))) functional_prediction_sample[1:nrow(ko_list)] = 0
    functional_prediction = cbind(functional_prediction, functional_prediction_sample)
  }
  colnames(functional_prediction) = names(otu_table)[2:ncol(otu_table_reduced_aggregated)]
  functional_prediction_final = data.frame(KO = ko_list$ko, functional_prediction, description = ko_list$description,check.names = F)
  if(ncol(functional_prediction) >= 2) keep = which(rowSums(functional_prediction) > 0)
  if(ncol(functional_prediction) == 1) keep = which(functional_prediction > 0)
  if (length(keep) == 0) stop("No functional prediction possible!\nEither no nearest neighbor found or your table is empty!")
  functional_prediction_final = functional_prediction_final[keep,]
 rownames(  functional_prediction_final)<-  functional_prediction_final[,1]
 tax4fun_KO_annotation<-unique(functional_prediction_final[,c("KO","description")])
functional_prediction_final$description <- functional_prediction_final$KO<-NULL
 fast.write(functional_prediction_final, file = 'tax4fun2_prediction/KOpred.csv')
   fast.write(tax4fun_KO_annotation, file = 'tax4fun_KO_annotation.csv', row.names = F)

  # Converting the KO profile to a profile of KEGG pathways
  message('Converting functions to pathways')
  ko2ptw = read.delim(file.path(path_to_reference_data, 'KEGG/ko2ptw.txt'))
  if(normalize_pathways) functional_prediction_norm = functional_prediction / ko_list$pathway_count
  pathway_prediction = aggregate(x = functional_prediction[ko2ptw$nrow,], by = list(ko2ptw$ptw), sum)
  if(ncol(pathway_prediction) >= 3){
    col_sums = colSums(pathway_prediction[,-1])
    col_sums[col_sums == 0] = 1
    pathway_prediction[,-1] = t(t(pathway_prediction[,-1]) / col_sums)
    keep = which(rowSums(pathway_prediction[,-1]) > 0)
  } else {
    pathway_prediction[,-1] = t(t(pathway_prediction[,-1]) / sum(pathway_prediction[,-1]))
    keep = which(pathway_prediction[,2] > 0)
  }
  if(sum(pathway_prediction[,-1]) == 0) stop("Conversion to pathway failed!")
  names(pathway_prediction) = names(otu_table)
  names(pathway_prediction)[1] = 'pathway'
  
  ptw_desc = read.delim(paste(path_to_reference_data, '/KEGG/ptw.txt', sep = ''))
  pathway_prediction_final = merge(pathway_prediction, ptw_desc)[keep,]
   rownames(pathway_prediction_final)<-pathway_prediction_final[,1]
    tax4fun_pathway_annotation<-unique(pathway_prediction_final[,c("pathway","level1","level2","level3")])
    fast.write(tax4fun_pathway_annotation, file = 'tax4fun_pathway_annotation.csv', row.names = F)
    pathway_prediction_final$pathway <- pathway_prediction_final$level1 <- pathway_prediction_final$level2 <- pathway_prediction_final$level3 <- NULL
    fast.write(pathway_prediction_final, file = 'tax4fun2_prediction/pathpred.csv')

}

##modify tax4fun to remove all 0 columns

Tax4Fun <- function(Tax4FunInput,folderReferenceData, fctProfiling=TRUE,refProfile="UProC",shortReadMode=TRUE,normCopyNo=TRUE){
  Tax4FunReferenceData <- importTax4FunReferenceData(folderReferenceData)
  ### clean the taxonomy names for match
  cleanedNms <- DeepCleanTaxaNames(rownames(Tax4FunInput$otuTable))

  #Intersect Mapping SILVA to KEGG and user OTU table
rank <- colnames(cleanedNms)[length(colnames(cleanedNms))]
commonOTUs <- intersect(Tax4FunReferenceData$SilvaTaxmat[,rank],cleanedNms[,rank])
if(length(commonOTUs)==0){
AddErrMsg("No feature in your data has been detected in Tax4Fun prediction database!");
 return(0);
}
indexInput <- match(commonOTUs,cleanedNms[,rank])
indexSILVAToKEGG <- match(commonOTUs,Tax4FunReferenceData$SilvaTaxmat[,rank])
subsetOTUTables <- as.matrix(Tax4FunInput$otuTable[indexInput,])
keepidx <-  colSums(subsetOTUTables)>0
subsetOTUTables <- subsetOTUTables[,keepidx]
subsetSILVAToKEGG <- Tax4FunReferenceData$SilvaToKEGGMappingMat[indexSILVAToKEGG,]
subsetSILVAToKEGG <- as.data.frame(as.matrix(subsetSILVAToKEGG))
   
  if(length(which(keepidx==F))>0){
    rmidx = which(keepidx==F)
    rmsmp = colnames(subsetOTUTables)[rmidx]
    if(length(rmidx)==1){
      current.msg<<-paste0("Sample ", rmsmp, " was removed during the processing due to the abundance of all the matched taxa are 0 in this sample.")
    }else{
      rmsmp <- paste(rmsmp,collapse = ", ")
      current.msg<<-paste0("Samples ", rmsmp, " were removed during the processing due to the abundance of all the matched taxa are 0 in these samples.")
    }
  }
 
  #Calculate taxonomic profile for KEGG organisms
  if(fctProfiling){
    FctCat <- Tax4FunReferenceData$KEGGKOInformation
    if(refProfile=="UProC"){
      if(shortReadMode){
        RefProfile <- Tax4FunReferenceData$FctAbundancesKEGGBacArchUProCShort
      }else{
        RefProfile <- Tax4FunReferenceData$FctAbundancesKEGGBacArchUProCLong
      }
    }else if(refProfile=="PAUDA"){
      if(shortReadMode){
        RefProfile <- Tax4FunReferenceData$FctAbundancesKEGGBacArchPAUDAShort
      }else{
        RefProfile <- Tax4FunReferenceData$FctAbundancesKEGGBacArchPAUDALong
      }
    }else{
      print("Invalid functional profiling method. Using default functional profiling method.")
      if(shortReadMode){
        RefProfile <- Tax4FunReferenceData$FctAbundancesKEGGBacArchUProCShort
      }else{
        RefProfile <- Tax4FunReferenceData$FctAbundancesKEGGBacArchUProCLong
      }
    }
    
  }else{
    RefProfile <- Tax4FunReferenceData$PathwayAbundancesKEGGBacArch
    FctCat <- Tax4FunReferenceData$PathwayInformationKEGGBacArch
  }
  
  Tax4FunProfile <- matrix(data=0,nrow=ncol(subsetOTUTables),ncol=nrow(RefProfile))
  for(i in 1:ncol(subsetOTUTables)){
    KEGGTaxProfile <- t(subsetSILVAToKEGG) %*% subsetOTUTables[,i]
    #Normalize by KEGG copy number
    if(normCopyNo){
      NormKEGGTaxProfile <- KEGGTaxProfile/Tax4FunReferenceData$KEGGBacArchCopyNumbers
    }else{
      NormKEGGTaxProfile <- as.data.frame(KEGGTaxProfile)
    }
    NormKEGGTaxProfile <- NormKEGGTaxProfile/sum(NormKEGGTaxProfile)
    #Calculate TaxFun profile
    Tax4FunProfile[i,] <- as.matrix(RefProfile) %*% NormKEGGTaxProfile$V1
  }
  
  if(fctProfiling){
    ###functional profiling
    FctCat <- FctCat[which(!colSums(Tax4FunProfile)==0),c(1,2)]
    Tax4FunProfile <- Tax4FunProfile[,which(!colSums(Tax4FunProfile)==0)]
    Tax4FunProfile <- as.matrix(Tax4FunProfile)
    
    if(ncol(Tax4FunProfile)>1){
      colnames(Tax4FunProfile) <- paste(FctCat$V2,FctCat$V1,sep="; ")
      rownames(Tax4FunProfile) <- colnames(subsetOTUTables)
    }else{
      Tax4FunProfile <- as.matrix(t(Tax4FunProfile))
      colnames(Tax4FunProfile) <- paste(FctCat$V2,FctCat$V1,sep="; ")
      rownames(Tax4FunProfile) <- Tax4FunInput$sampleNames
    }
  }else{
    ###metabolic profiling
    FctCat <- FctCat[which(!colSums(Tax4FunProfile)==0),c(1,2)]
    Tax4FunProfile <- Tax4FunProfile[,which(!colSums(Tax4FunProfile)==0)]
    Tax4FunProfile <- as.matrix(Tax4FunProfile)
    
    if(ncol(Tax4FunProfile)>1){
      colnames(Tax4FunProfile) <- paste(FctCat$V3,FctCat$V4,sep="; ")
      rownames(Tax4FunProfile) <- Tax4FunInput$sampleNames
    }else{
      Tax4FunProfile <- as.matrix(t(Tax4FunProfile))
      colnames(Tax4FunProfile) <- paste(FctCat$V3,FctCat$V4,sep="; ")
      rownames(Tax4FunProfile) <- Tax4FunInput$sampleNames
    }
  }
  
  #colnames(Tax4FunProfile) <- paste(FctCat$V2,FctCat$V1,sep="; ")
  #rownames(Tax4FunProfile) <- Tax4FunInput$sampleNames
  FTU <- 1-colSums(subsetOTUTables)/colSums(Tax4FunInput$otuTable[,keepidx])
  names(FTU) <- colnames(subsetOTUTables)
  Tax4FunProfile <- list(Tax4FunProfile=Tax4FunProfile, FTU=FTU,fctProfiling=fctProfiling,refProfile=refProfile,shortReadMode=shortReadMode)
  class(Tax4FunProfile) <- "Tax4Fun"
  return(Tax4FunProfile)
} 


DeepCleanTaxaNames <- function(nms){
  feat_nm =data.frame(nms,check.names=FALSE);
  names(feat_nm)<-"Rank";
  taxonomy<-suppressWarnings(splitstackshape::cSplit(feat_nm,"Rank",";"));
  taxmat= data.frame(matrix(NA, ncol = 7, nrow = nrow(taxonomy)),check.names=FALSE);
  colnames(taxmat) <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species");
  taxmat[,1:ncol(taxonomy)]<-taxonomy;
  taxmat<-taxmat[colSums(!is.na(taxmat))>0];
  if("Species" %in% colnames( taxmat)){
    sps = data.frame(taxmat[,c('Genus','Species')])
    if(grepl("^g__",sps[1,"Genus"])){
      sps$Species = gsub("^s__","",sps$Species)
      if(all(!grepl(" ",sps$Species)) & all(!grepl("_",sps$Species))){
        sps$Genus =gsub("^g__","",sps$Genus)
        idxsp = which(!is.na(sps$Genus) & !is.na(sps$Species) & sps$Genus!="" & sps$Species!="")
        sps$Species[idxsp] = paste0("s__",sps$Genus[idxsp],"_",sps$Species[idxsp])
      }
    }else if(grepl("^g_",sps[1,"Genus"])){
      sps$Species = gsub("^s_","",sps$Species)
      if(all(!grepl(" ",sps$Species)) & all(!grepl("_",sps$Species))){
        sps$Genus =gsub("^g_","",sps$Genus)
        idxsp = which(!is.na(sps$Genus) & !is.na(sps$Species) & sps$Genus!="" & sps$Species!="")
        sps$Species[idxsp] = paste0("s__",sps$Genus[idxsp],"_",sps$Species[idxsp])
      }
    }else if(grepl("^g_0__",sps[1,"Genus"])){
      sps$Species = gsub("^s_0__","",sps$Species)
      if(all(!grepl(" ",sps$Species)) & all(!grepl("_",sps$Species))){
        sps$Genus =gsub("^g_0__","",sps$Genus)
        idxsp = which(!is.na(sps$Genus) & !is.na(sps$Species) & sps$Genus!="" & sps$Species!="")
        sps$Species[idxsp] = paste0("s_0__",sps$Genus[idxsp],"_",sps$Species[idxsp])
      }
    }else{
      if(all(!grepl(" ",sps$Species)) & all(!grepl("_",sps$Species))){
        idxsp = which(!is.na(sps$Genus) & !is.na(sps$Species) & sps$Genus!="" & sps$Species!="")
        sps$Species[idxsp] = paste0(sps$Genus[idxsp],"_",sps$Species[idxsp])
      }
    }
    taxmat[,'Species'] <- sps$Species
  }
  
  taxmat <- apply(taxmat,2,function(x) gsub("'","",x))
  taxmat <- apply(taxmat,2,function(x) gsub("\\[|\\]","",x))
  
  return(taxmat)

}

importTax4FunReferenceData <- function(folder){
  
  if(substr(folder,nchar(folder),nchar(folder))=="/"){
    pathReferenceData <- folder
  }else{
    pathReferenceData <- paste(folder,"/",sep="")
  }
    
  referenceData <- list()

  tmpReferenceData <- readRDS(paste(pathReferenceData,"PathwayAbundancesKEGGBacArch.RData",sep=""))
  referenceData$PathwayAbundancesKEGGBacArch <- tmpReferenceData
  tmpReferenceData <- readRDS(paste(pathReferenceData,"PathwayInformationKEGGBacArch.RData",sep=""))
  referenceData$PathwayInformationKEGGBacArch <- tmpReferenceData
  
  tmpReferenceData <- readRDS(paste(pathReferenceData,"KEGGKOInformation.RData",sep=""))
  referenceData$KEGGKOInformation <- tmpReferenceData

  tmpReferenceData <- readRDS(paste(pathReferenceData,"KEGGBacArchCopyNumbers.RData",sep=""))
  referenceData$KEGGBacArchCopyNumbers <- tmpReferenceData
  
  tmpReferenceData <- readRDS(paste(pathReferenceData,"FctAbundancesKEGGBacArchPAUDALong.RData",sep=""))
  referenceData$FctAbundancesKEGGBacArchPAUDALong <- tmpReferenceData
  tmpReferenceData <- readRDS(paste(pathReferenceData,"FctAbundancesKEGGBacArchPAUDAShort.RData",sep=""))
  referenceData$FctAbundancesKEGGBacArchPAUDAShort <- tmpReferenceData
  tmpReferenceData <- readRDS(paste(pathReferenceData,"FctAbundancesKEGGBacArchUProCLong.RData",sep=""))
  referenceData$FctAbundancesKEGGBacArchUProCLong <- tmpReferenceData
  tmpReferenceData <- readRDS(paste(pathReferenceData,"FctAbundancesKEGGBacArchUProCShort.RData",sep=""))
  referenceData$FctAbundancesKEGGBacArchUProCShort <- tmpReferenceData
  
  
  tmpReferenceData <- readRDS(paste(pathReferenceData,"SilvaToKEGGMappingMat.RData",sep=""))
  referenceData$SilvaToKEGGMappingMat <- tmpReferenceData
  tmpReferenceData <- readRDS(paste(pathReferenceData,"SilvaIDs.RData",sep=""))
  referenceData$SilvaIDs <- tmpReferenceData
  tmpReferenceData <- readRDS(paste(pathReferenceData,"SilvaTaxmat.RData",sep=""))
  referenceData$SilvaTaxmat <- tmpReferenceData

  return(referenceData)
}

# helper functions for local testing
get.fun.lib.path <- function(type){

 type <- tolower(type);

 # if(.on.public.web){
  if(file.exists("/home/glassfish/resources/MicrobiomeAnalyst")){ # on server
    if(type == "tax4fun"){
        return("/home/glassfish/resources/MicrobiomeAnalyst/tax4fun/SILVA123");
    }else if(type == "tax4fun2"){
        return("/home/glassfish/resources/MicrobiomeAnalyst/tax4fun2");
    }else if(type == "picrust12"){
        return("/home/glassfish/resources/MicrobiomeAnalyst/picrust12");
    }else{
        return("/home/glassfish/resources/MicrobiomeAnalyst/picrust13");
    }
  }
  # xia laptop
  if(file.exists("/Users/jeffxia/Dropbox/resources")){
    if(type == "tax4fun"){
        return("/Users/jeffxia/Dropbox/resources/MicrobiomeAnalyst/tax4fun/SILVA123");
    }else if(type == "tax4fun2"){
        return("/Users/jeffxia/Dropbox/resources/MicrobiomeAnalyst/tax4fun2");
    }else if(type == "picrust12"){
        return("/Users/jeffxia/Dropbox/resources/MicrobiomeAnalyst/picrust12");
    }else{
        return("/Users/jeffxia/Dropbox/resources/MicrobiomeAnalyst/picrust13");
    }
  }

  # add your path here

  # yao laptop
  if(file.exists("/Users/lzy/Documents")){
    if(type == "tax4fun"){
        return("/Users/lzy/Documents/MicrobiomeAnalystR-master/function_database/SILVA123");
    }else if(type == "tax4fun2"){
        return("/Users/lzy/Documents/MicrobiomeAnalystR-master/function_database/Tax4Fun2");
    }else if(type == "picrust12"){
        return("/Users/lzy/Documents/MicrobiomeAnalystR-master/function_database/picrust12");
    }else{
        return("/Users/lzy/Documents/MicrobiomeAnalystR-master/function_database/picrust13");
    }
  }

}
