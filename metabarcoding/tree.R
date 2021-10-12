## ---- eval = FALSE-------------------------------------------------------------------------------
## #seqtab<- readRDS("seqtab_V4.rds")
## 
## seqs <- getSequences(seqtab)
## 
## # This propagates to the tip labels of the tree.
## # At this stage ASV labels are full ASV sequence
## names(seqs) <- seqs
## alignment <- AlignSeqs(DNAStringSet(seqs),
##                        anchor = NA,
##                        iterations = 5,
##                        refinements = 5)
## 
## print("Computing pairwise distances from ASVs...")
## Sys.time()
## phang.align <- phyDat(as(alignment, "matrix"), type = "DNA")
## dm <- dist.ml(phang.align)
## treeNJ <- NJ(dm) # Note, tip order is not sequence order
## fit = pml(treeNJ, data = phang.align)
## print("...Done.")
## Sys.time()
## 
## print("Fit GTR model...")
## Sys.time()
## fitGTR <- update(fit, k = 4, inv = 0.2)
## print("...Done.")
## Sys.time()
## 
## print("Computing likelihood of tree...")
## Sys.time()
## fitGTR <- optim.pml(fitGTR,
##                     model = "GTR",
##                     optInv = TRUE,
##                     optGamma = TRUE,
##                     rearrangement = "stochastic",
##                     control = pml.control(trace = 0))
## print("...Done".)
## Sys.time()