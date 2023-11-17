
my.sparcc.net <- function(mbSetObj=NULL, corr.net.name, networkType="static", netLayout="kk", netTextSize){
  load_igraph()
  mbSetObj <- .get.mbSetObj(mbSetObj);
  abundOpt <- mbSetObj$analSet$abund.opt;
  edge.list <- mbSetObj$analSet$network_cor;
  edge.list <- edge.list[,c("Taxon1", "Taxon2"), drop=FALSE];
  g <- graph_from_data_frame(edge.list, directed = FALSE, vertices = NULL);
  E(g)$weight <- abs(mbSetObj$analSet$network_cor[, "Correlation"]);
  E(g)$correlation <- mbSetObj$analSet$network_cor[, "Correlation"];
  colnames(edge.list) <- c("Source", "Target");
  nodes <- unique(c(edge.list[,1], edge.list[,2]));
  node.list <- data.frame(Id=nodes, Name=nodes,check.names=FALSE); 
  abundance_table <- mbSetObj$analSet$boxdatacor[,-length(colnames(mbSetObj$analSet$boxdatacor))]
  ab <- apply(abundance_table,2,function(x){sum(x)/length(x)})  
  taxa <- mbSetObj$analSet$filt.taxa.table;
  
  colorOpt <- mbSetObj$analSet$corr_color_opt;
  # annotation
  nms <- V(g)$name;
  inx <- !nms %in% taxa[,length(colnames(taxa))]
  toDelete = nms[inx]
  g <- delete_vertices(g, toDelete);
  nms <- V(g)$name;
  taxa <- taxa[match(nms, taxa[,length(colnames(taxa))]),] 
  if(all(c(!colorOpt %in% colnames(taxa), colorOpt != "expr"))){
    
    current.msg <<- "Invalid taxa is selected for coloring (must be same or higher taxonomy level of selected taxa used for correlation calculation)"
    return(0);
  }
  hit.inx <- match(nms, node.list[,1]);
  lbls <- node.list[hit.inx, 2];
  
  # setup shape (gene circle, other squares)
  shapes <- rep("circle", length(nms));
  itypes <- rep("circle", length(nms));
  seeds <- rep("circle", length(nms)); 
  # get edge data
  edge.mat <- get.edgelist(g);
  edge.mat1 <- data.frame(edge.mat,check.names=FALSE);
  edge.mat1$color <- ComputeColorGradientCorr(E(g)$correlation);
  edge.mat1 <- as.matrix(edge.mat1)

  edge.mat <- cbind(id=1:nrow(edge.mat), source=edge.mat[,1], target=edge.mat[,2], color = edge.mat1[,3], weight=E(g)$weight, correlation = E(g)$correlation);
  # now get coords
  pos.xy <- PerformLayOut(g);
 
  # get the note data
  node.btw <- as.numeric(betweenness(g));
  node.eig <- eigen_centrality(g);
  node.eig <- as.numeric(node.eig$vector);
  node.tra <- transitivity(g,type=c("local"))
  node.dgr <- as.numeric(igraph::degree(g)); 
  node.exp <- as.numeric(get.vertex.attribute(g, name="abundance", index = V(g)));
  # node size to abundance values
  if(vcount(g)>500){
    min.size = 1;
  }else if(vcount(g)>200){
    min.size = 2;
  }else{
    min.size = 2;
  } 
  abundance <- unname(ab[nms])
  node.sizes <- as.numeric(rescale2NewRange((log(abundance+1))^2, min.size, 15));

  # update node color based on betweenness
  require("RColorBrewer");
  topo.val <- log(node.btw+1);
  exp_table <- mbSetObj$analSet$diff_table;
  colVecNms <- rep("NA", length(topo.val))
  nms2=nms
 
  if(colorOpt == "expr"){
    nms1 = strsplit(nms, "__")
    if(length(nms1[[2]])>1){
      l = unlist(lapply(nms1, function(x) unlist(x[2])));
      inx = is.na(l)
      if(length(which(inx == T))>0){
        l[inx]= paste0(nms1[[which(inx == T)]][1],"__")
      }
      nms2=l
    }else{
      nms2 = nms
    }
    exp_table = exp_table[match(nms2, exp_table$tax_name), ]
    topo.colsb <- topo.colsb1 <-ComputeColorGradientCorr(exp_table$median_diff);
    topo.colsw <- topo.colsw1 <-ComputeColorGradientCorr(exp_table$median_diff);
    
    smpl <- data.frame(sample_data(mbSetObj$dataSet$proc.phyobj),check.names=FALSE);
    subsmpl <- smpl[mbSetObj$dataSet$selected.grps,,drop=FALSE]
    if(mbSetObj$dataSet$cor.method == "sparcc"){
      abund_data <- qs::qread("sparcc_data.qs")
    }else if(mbSetObj$dataSet$cor.method %in% c("secom_p1","secom_p2","secom_dist")){
      abund_data <- qs::qread("secom_data.qs")
    }else{
      abund_data <- t(as.matrix(mbSetObj$analSet$netcorr_data))
    }

    grp.vec = factor(subsmpl[,mbSetObj$dataSet$meta])
    filt_data = abund_data[,mbSetObj$dataSet$selected.grps]
    filt_data = as.data.frame(t(filt_data),check.names=FALSE);
    meta_vec = rep(grp.vec, ncol(filt_data))
    r = stack(filt_data)
    r$meta = meta_vec
    
    smp.nms = as.character(unique(r$ind))
    meta.nms = as.character(unique(r$meta))
    
    abundColor = generateColorArr(length(meta.nms));
    
    propList = list()
    for(i in 1:length(smp.nms)){
      smp.nm = as.character(smp.nms[i])
      for(j in 1:length(meta.nms)){
        meta.nm = as.character(unique(r$meta)[j])
        tmp = r[which(r$ind == smp.nm),]
        tmp = tmp[which(tmp$meta == meta.nm),]
        if(abundOpt == "meanlog"){
          vals = log(tmp$values, 10)
          average = mean(vals)
        }else if(abundOpt == "medianlog"){
          vals = log(tmp$values, 10)
          average = median(vals)
        }else if(abundOpt == "mean"){
          total = sum(tmp$values)
          average = total/length(tmp$values)
        }else{
          average=median(tmp$values)
        }
        if(length(propList[[smp.nm]]) == 0){
          propList[[smp.nm]] = list()
        }
        propList[[smp.nm]][[meta.nm]]=average
      }
    }
    
  }else if(colorOpt == "metadata"){

  }else{

    color_var <- as.character(unique(taxa[,colorOpt]));
    colVec=generateColorArr(length(color_var));
    names(colVec) = unique(taxa[,colorOpt])
    colVecNms<-names(colVec[as.character(taxa[,colorOpt])])
    topo.colsb  <- unname(colVec[as.character(taxa[,colorOpt])])
    topo.colsw  <- topo.colsb
    propList = "NA"
    abundColor = "NA"
  }
  # color based on expression
  bad.inx <- is.na(node.exp) | node.exp==0;
  
  if(!all(bad.inx)){
    exp.val <- node.exp;
    node.colsb.exp <- ComputeColorGradientCorr(exp.val); 
    node.colsw.exp <- ComputeColorGradientCorr(exp.val);
    node.colsb.exp[bad.inx] <- "#d3d3d3"; 
    node.colsw.exp[bad.inx] <- "#c6c6c6"; 
  }else{
    node.colsb.exp <- rep("#D3D3D3",length(node.dgr)); 
    node.colsw.exp <- rep("#C6C6C6",length(node.dgr)); 
  }
  
  network_prop = list();
  for(i in 1:length(node.sizes)){
    network_prop[[i]]  <- list(
      eigen = node.eig[i],
      transitivity = node.tra[i]
    )
  }
  
  type <- rep(FALSE,length(node.dgr));
  node_attr <- list.vertex.attributes(g);
  node_attr <- node_attr[which(node_attr!="names")] 
  
  # now create the json object
  nodes <- vector(mode="list");

  for(i in 1:length(node.sizes)){
    nodes[[i]] <- list(
      id=nms[i], 
      idx=i,
      label=nms2[i],
      size=node.sizes[i], 
      tsize=node.sizes[i],
      type="gene",
      itype=itypes[i],
      taxon=colVecNms[i],
      color=topo.colsb[i],
      colorb=topo.colsb[i],
      colorw=topo.colsw[i],
      topocolb=topo.colsb[i],
      topocolw=topo.colsw[i],
      expcolb=node.colsb.exp[i],
      expcolw=node.colsw.exp[i],
      x=pos.xy[i,1],
      y= pos.xy[i,2],
      user =network_prop[[i]],
      attributes=list(
        degree=node.dgr[i], 
        between=node.btw[i],
        expr = exp_table$mean_diff[i],
        eigen = node.eig[i],
        transitivity = node.tra[i]
      )
    );
  }

  if(!.on.public.web){
    if(networkType == "static"){
      
      # static plot
      load_ggraph()
      
      p <- ggraph(g, layout = netLayout) + theme_void() +
        geom_edge_fan(color="gray20", width=0.5, alpha=0.5) +
        geom_node_point(color=topo.colsb, size=node.sizes, alpha = 0.95) +
        geom_node_text(aes(label = V(g)$name), size = netTextSize, repel=TRUE, nudge_y = 0.05, nudge_x = 0.05, check_overlap = TRUE) +
        # 10% space above/to the side of x
        scale_y_continuous(expand = expansion(mult = c(.1, .1))) +
        scale_x_continuous(expand = expansion(mult = c(.1, .1)))
      
      filename <- paste0(gsub(".json", "", corr.net.name), ".png")
      ggsave(p, file=filename, width = 7, height = 7)
      
    }else{ 
      
      # interactive plot
      load_visNetwork()
      
      data <- toVisNetworkData(g)
      data[["nodes"]]$color <- topo.colsb
      data[["nodes"]]$size <- node.sizes
      network <- visNetwork(nodes = data$nodes, edges = data$edges, 
                            idToLabel = TRUE, height = "900px", width = "100%") %>%
        visEdges(color = list(color = "lightgrey", highlight = "red"))
      filename <- paste0(gsub(".json", "", corr.net.name), ".html")
      visSave(network, file = filename)
    }
  }
  # save node table
  nd.tbl <- data.frame(Id=nms, Label=lbls, Degree=node.dgr, Betweenness=round(node.btw,2),check.names=FALSE);
  # order 
  ord.inx <- order(nd.tbl[,3], nd.tbl[,4], decreasing = TRUE)
  nd.tbl <- nd.tbl[ord.inx, ];
  fast.write(nd.tbl, file="node_table.csv", row.names=FALSE);
  
  if(.on.public.web){
    # covert to json
    taxNodes <- list(tax=colVecNms,colors=topo.colsb);
    netData <- list(nodes=nodes,edges=edge.mat, abund=propList, abundColor=abundColor, taxa=taxa, taxaNodes = taxNodes)
    sink(corr.net.name);
    cat(RJSONIO::toJSON(netData));
    sink();
  }
  return(.set.mbSetObj(mbSetObj));
}
