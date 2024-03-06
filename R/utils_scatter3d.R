##################################################
## R script for MicrobiomeAnalyst
## Description: Compute json file for 3D scatter plot
## Author: Guangyan Zhou, guangyan.zhou@mail.mcgill.ca
##################################################

my.json.scatter <- function(mbSetObj=NA, filenm, containsLoading=F){
  library(igraph);
  ## 
  mbSetObj <- .get.mbSetObj(mbSetObj);
  res <- mbSetObj$analSet$pca$score
  nodes <- vector(mode="list");
  names <- res$name;
  if(ncol(res$xyz) > nrow(res$xyz)){
    orig.pos.xyz <- t(res$xyz);
  }else{
    orig.pos.xyz <- res$xyz;
  }
  
  ticksX <- pretty(range(orig.pos.xyz[,1]*1.2),10, bound=F);
  ticksY <- pretty(range(orig.pos.xyz[,2]*1.2),10, bound=F);
  ticksZ <- pretty(range(orig.pos.xyz[,3]*1.2),10, bound=F);
  
  #add two nodes, 1 with all min values and another with all max values. For scaling purposes
  minNode <-  c(min(ticksX), min(ticksY), min(ticksZ));
  maxNode <-  c(max(ticksX), max(ticksY), max(ticksZ));
  
  # Add the new rows to the data frame
  orig.pos.xyz.mod <- rbind(orig.pos.xyz, minNode);
  orig.pos.xyz.mod <- rbind(orig.pos.xyz.mod, maxNode);
  
  #scaling
  pos.xyz <- orig.pos.xyz.mod;
  pos.xyz[,1] <- scale_range(orig.pos.xyz.mod[,1], -1,1);
  pos.xyz[,2] <- scale_range(orig.pos.xyz.mod[,2], -1,1);
  pos.xyz[,3] <- scale_range(orig.pos.xyz.mod[,3], -1,1);
  
  #remove last two rows
  pos.xyz <- pos.xyz[1:(nrow(pos.xyz) - 2), ]
  
  metadf <- res$facA
  if(!is.factor(metadf)){
    metadf <- as.factor(metadf);
  }
  col = vector();
  
  meta.vec <- metadf
  meta.vec.num = as.integer(metadf)
  
  
  if(all(sapply(meta.vec, function(x) !is.na(as.numeric(x))))){
    sorted_unique_meta_vec <- sort(unique(as.numeric(levels(meta.vec))))
    col.s <- generateColorGradient(sorted_unique_meta_vec)
    # Create a named vector where names are sorted unique meta.vec values and values are corresponding colors
    col <- generateColorGradient(meta.vec)
    # Prepare legend data
    legendData <- list(label = sorted_unique_meta_vec, color = col.s)
    
  }else{
    col.s <- gg_color_hue(length(levels(metadf)))
    for(i in 1:length(meta.vec.num)){
      col[i] = col.s[meta.vec.num[i]];
    }
    legendData <- list(label=unique(meta.vec),color=col.s)
    
  }
  
  
  if("facB" %in% names(res)){
    meta.vec2 <- res$facB
    phyobjdata <- qs::qread("merged.data.qs");
    metadf <- as.matrix(sample_data(phyobjdata));
    metadf <- as.data.frame(metadf);
    shape <- vector();
    meta.vec.num <- as.integer(as.factor(meta.vec2))
    shape.s <- c("circle", "triangle", "square", "diamond")[1:length(unique(meta.vec2))];
    
    for(i in 1:length(meta.vec.num)){
      shape[i] = shape.s[meta.vec.num[i]];
    }
    legendData2 <- list(label=unique(meta.vec2),shape=shape.s);
  }
  
  nodeSize = 18;
  
  for(i in 1:length(names)){
    nodes[[i]] <- list(
      id=names[i],
      label=names[i],
      size=nodeSize,
      meta=meta.vec[i],
      fx = unname(pos.xyz[i,1])*1000,
      fy = unname(pos.xyz[i,2])*1000,
      fz = unname(pos.xyz[i,3])*1000,
      origX = unname(orig.pos.xyz[i,1]),
      origY = unname(orig.pos.xyz[i,2]),
      origZ = unname(orig.pos.xyz[i,3]),
      colorb=col[i]
    );
    
    if("facB" %in% names(res)){
      nodes[[i]][["meta2"]] <- meta.vec2[i]
      nodes[[i]][["shape"]] <- shape[i]
      
    }
  }
  
  
  ticks <- list(x=ticksX, y=ticksY, z=ticksZ);
  library(RJSONIO)
  
  if(!containsLoading){
    netData <- list(nodes=nodes,
                    edges="NA",
                    meta=metadf, 
                    loading="NA",
                    axis=res$axis, 
                    ticks=ticks,
                    metaCol = legendData);
  }else{
    res2 <- qs::qread("loading3d.qs");
    
    if(ncol(res2$xyz) > nrow(res2$xyz)){
      orig.load.xyz <- t(res2$xyz);
    }else{
      orig.load.xyz <- res2$xyz;
    }
    
    ticksX <- pretty(range(orig.load.xyz[,1]),10);
    ticksY <- pretty(range(orig.load.xyz[,2]),10);
    ticksZ <- pretty(range(orig.load.xyz[,3]),10);
    
    #add two nodes, 1 with all min values and another with all max values. For scaling purposes
    minNode <-  c(min(ticksX), min(ticksY), min(ticksZ));
    maxNode <-  c(max(ticksX), max(ticksY), max(ticksZ));
    
    # Add the new rows to the data frame
    orig.load.xyz.mod <- rbind(orig.load.xyz, minNode)
    orig.load.xyz.mod <- rbind(orig.load.xyz.mod, maxNode)
    
    load.xyz <- orig.load.xyz.mod;
    
    load.xyz[,1] <- scale_range(orig.load.xyz.mod[,1], -1,1);
    load.xyz[,2] <- scale_range(orig.load.xyz.mod[,2], -1,1);
    load.xyz[,3] <- scale_range(orig.load.xyz.mod[,3], -1,1);
    
    #remove last two rows
    load.xyz <- load.xyz[1:(nrow(load.xyz) - 2), ];
    
    names <- res2$name;
    if("entrez" %in% names(res2)){
      ids <- res2$entrez;
    }else{
      ids <- res2$name;
    }
    colres <- rgba_to_hex_opacity(res2$cols);
    colorb <- colres[[1]];
    opacity_array <- colres[[2]];
    nodes2 <- vector(mode="list");
    
    for(i in 1:length(names)){
      nodes2[[i]] <- list(
        id=ids[i],
        label=names[i],
        size=24,
        opacity=opacity_array[i],
        fx = unname(load.xyz[i,1])*1000,
        fy = unname(load.xyz[i,2])*1000,
        fz = unname(load.xyz[i,3])*1000,
        origX = unname(orig.load.xyz[i,1]),
        origY = unname(orig.load.xyz[i,2]),
        origZ = unname(orig.load.xyz[i,3]),
        seedArr = "notSeed",
        colorb=colorb[i],
        colorw=colorb[i]
      );
      
    }
    
    ticksLoading <- list(x=ticksX, y=ticksY, z=ticksZ);
    
    netData <- list(omicstype="", 
                    nodes=nodes,
                    meta=metadf, 
                    loading=nodes2, 
                    misc="", ticks=ticks,
                    ticksLoading=ticksLoading,
                    axis=res$axis,
                    axisLoading=res2$axis, 
                    metaCol = legendData);
    #res2$pos.xyz <- load.xyz;
    #qs::qsave(res2, "pca3d_loadings.qs");
  }
  
  if("facB" %in% names(res)){
    netData$metaShape <- legendData2;
  }
  
  rownames(pos.xyz) <- res$name;
  mbSetObj$analSet$pos.xyz <- pos.xyz;
  qs::qsave(pos.xyz,"pos.xyz.qs");
  
  
  sink(filenm);
  cat(toJSON(netData));
  sink();
  return(.set.mbSetObj(mbSetObj));
}



# Define a function to convert RGBA to Hex and opacity values
rgba_to_hex_opacity <- function(rgba_array) {
  
  # Define an empty vector to store the hex color values
  hex_values <- vector("character", length = length(rgba_array))
  
  # Define an empty vector to store the opacity values
  opacity_values <- vector("numeric", length = length(rgba_array))
  
  # Loop through each element in the input array
  for (i in seq_along(rgba_array)) {
    rgba <- rgba_array[i]
    rgba <- gsub("rgba\\(", "",rgba);
    rgba <- gsub("\\)", "",rgba);
    # Convert the RGBA value to a list of numeric values
    rgba <- strsplit(rgba, ",")[[1]]
    
    rgba <- as.numeric(rgba)
    # Convert the RGB values to hexadecimal notation
    hex <- rgb(rgba[1], rgba[2], rgba[3], maxColorValue = 255)
    hex_values[i] <- hex
    
    # Extract the opacity value from the RGBA string
    opacity_values[i] <- rgba[4]
  }
  
  # Return a list containing the hex color values and opacity values
  return(list(hex_values = hex_values, opacity_values = opacity_values))
}


scale_range <- function(x, new_min = 0, new_max = 1) {
  range <- pretty(x,10);
  old_min <- min(range);
  old_max <- max(range);
  
  (x - old_min) / (old_max - old_min) * (new_max - new_min) + new_min
}


ComputeEncasing <- function(filenm, type, names.vec, level=0.95, omics="NA"){
  
  
  level <- as.numeric(level)
  names = strsplit(names.vec, "; ")[[1]]
  
  pos.xyz <- qs::qread("pos.xyz.qs");
  #print(head(pos.xyz));
  inx = rownames(pos.xyz) %in% names;
  coords = as.matrix(pos.xyz[inx,c(1:3)])
  mesh = list()
  if(type == "alpha"){
    library(alphashape3d)
    library(rgl)
    sh=ashape3d(coords, 1.0, pert = FALSE, eps = 1e-09);
    mesh[[1]] = as.mesh3d(sh, triangles=T);
  }else if(type == "ellipse"){
    library(rgl);
    pos=cov(coords, y = NULL, use = "everything");
    mesh[[1]] = ellipse3d(x=as.matrix(pos), level=level);
  }else{
    library(ks);
    res=kde(coords);
    r = plot(res, cont=level*100, display="rgl");
    sc = scene3d();
    mesh = sc$objects;
  }
  library(RJSONIO);
  sink(filenm);
  cat(toJSON(mesh));
  sink();
  return(filenm);
}


gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}

generateColorGradient <- function(values, low="green", high="red") {
  # Normalize values to the range from 0 to 1
  values <- as.numeric(values);
  normalized_values <- (values - min(values)) / (max(values) - min(values))
  
  # Generate colors
  colors <- grDevices::colorRampPalette(c(low, high))(length(values))
  
  # Map normalized values to colors
  return(colors[as.numeric(cut(normalized_values, breaks = length(values), labels = FALSE))])
}

