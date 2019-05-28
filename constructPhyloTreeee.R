
constructPhyloTree <- function(constrcutSeqTabRes = constrcutSeqTabRes,
                         ...){
  #tic("constructPhyloTree")
  require(DECIPHER);
  require(phangorn);
  
  #get data;
  seqtab.nochim = constructSeqTabRes$seqtab.nochim;
  
  print("get sequences");
  seqs <- getSequences(seqtab.nochim);
  names(seqs) <- seqs;
  print("do alignment");
  alignment <- AlignSeqs(DNAStringSet(seqs), 
                         anchor = NA,
                         verbose = TRUE);
  print("construct tree");
  phangAlign <- phyDat(as(alignment, "matrix"), type = "DNA")
  dm <- dist.ml(phangAlign)
  treeNJ <- NJ(dm) # Note, tip order != sequence order
  fit = pml(treeNJ, data = phangAlign)
  fitGTR <- update(fit, k = 4, inv = 0.2)
  fitGTR <- optim.pml(fitGTR, model = "GTR", 
                      optInv = TRUE, 
                      optGamma = TRUE,
                      rearrangement = "stochastic", 
                      control = pml.control(trace = 0));
  ifelse(!dir.exists("phylogenetic_tree"),
         dir.create("phylogenetic_tree"),
         FALSE);
  write.tree(fitGTR$tree,
             file = file.path("phylogenetic_tree", "phylotree.tre"));
  #detach("package:phangorn", unload=TRUE);
  #print(toc());
}

