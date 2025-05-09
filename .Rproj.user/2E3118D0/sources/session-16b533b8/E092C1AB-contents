# copyright: Wei Sheng (wsheng@stu.njmu.edu.cn)

#' Scale gene expression for scGSEA
#'
#' This function normalizes single cell gene expression levels acorss cells for each gene. 
#'
#' @param obj either a `Seurat` or a `SingleCellExperiment` object
#' @param gene_expr_frac fliltering out genes if the non-zero fraction less than this number
#' @param kNN the k nearest neighbors for data smoothing/denoising. 
#' @return the obj with added data in slot `assay` of `Seurat` and `metadata` of `SingleCellExperiment`
#' @author Wei Sheng (wsheng at stu.njmu.edu.cn)
#' @examples obj <- scGSEA.scale.data(obj, kNN=51)
#' \dontrun{
#' }
#' @export
scGSEA.scale.data <- function(obj, gene_expr_frac=0.1, kNN=21) { 
  obj.class <- class(obj)
  stopifnot(obj.class %in% c("SingleCellExperiment", "Seurat"))
  if(obj.class == "Seurat") { 
    assay_name <- "RNA"
    assay_obj <- obj[[assay_name]]
    assay_class <- class(assay_obj)[1]
    if (assay_class == "Assay5") {
        print(assay_class)
        data <- as.matrix(GetAssayData(obj, assay = assay_name, layer = "data"))
      } else if (assay_class == "Assay") {
        print(assay_class)
        data <- as.matrix(GetAssayData(obj, assay = assay_name, slot = "data"))
      } else {
        stop("Unsupported assay class.")
      }
    #data <- as.matrix(obj@assays$RNA@data)

    gene.exp.idx <- apply(data>0, 1, sum) > gene_expr_frac * ncol(data) 
    gene.exp <- data[gene.exp.idx,]
    if(kNN>0) {
      if(all(dim(obj@graphs$RNA_snn) == ncol(data))) {
        gene.exp.sm <- do.call(cbind, sapply(1:ncol(gene.exp), function(i) {
          kNN.w <- tail(sort(obj@graphs$RNA_snn[i,]), n=(1+kNN))
          kNN.names <- names(kNN.w)
          gene.exp[, kNN.names] %*% kNN.w / sum(kNN.w)
        }, simplify=F))
        colnames(gene.exp.sm) <- colnames(gene.exp)
        gene.exp <- gene.exp.sm
      }
      else {
        warning("dim of obj@graphs$RNA_snn does not match the cell number in the object, kNN is set to 0!")
      }
    }
    gene.exp.cell_norm <- t(apply(gene.exp, 1, function(x) (x-mean(x))/sqrt(var(x))) )
    # obj@assays$scGSEA <- obj[["RNA"]]
    # obj@assays$scGSEA@data <- gene.exp.cell_norm
    obj[["scGSEA"]] <- CreateAssayObject(data = gene.exp.cell_norm)
  }
  if(obj.class == "SingleCellExperiment") {
    if(kNN>0) {
      warning("kNN>0 is not compatible with sce obj, kNN is set to 0!")
    }
    data <- as.matrix(logcounts(obj))
    gene.exp.idx <- apply(data>0, 1, sum) > gene_expr_frac * ncol(data) ## Para: expressed in the fraction of cell
    gene.exp <- data[gene.exp.idx,]
    gene.exp.cell_norm <- t(apply(gene.exp, 1, function(x) (x-mean(x))/sqrt(var(x))) )
    metadata(obj)[["scGSEA.scale.data"]] <- gene.exp.cell_norm 
  }
  return(obj)
}

#' Get gene list from a `SeqGeneSet` object
#' 
#' This function extracts a list of genes in geneSet(s) specificed by GeneSetName. 
#' 
#' @param gs.obj a `SeqGeneSet` object
#' @param GeneSetName Name(s) of geneSet(s)
#' @return a vector of gene names
#' @examples 
#' \dontrun{
#' }
#' @author Wei Sheng (wsheng at stu.njmu.edu.cn) 
#' @export
getGSgenelist <- function(gs.obj, GeneSetName) {
  stopifnot( is(gs.obj, "SeqGeneSet") & GeneSetName %in% gs.obj@GSNames )
  gs.obj@geneList [ gs.obj@GS[[ which( gs.obj@GSNames %in%  GeneSetName ) ]] ]
}

#' Load and store gene sets for scGSEA
#' 
#' This function loads gene sets and creats a `SeqGeneSet` object
#'
#' @param geneset.file the file path and name of GeneSet in .gmt format 
#' @param sc.obj either a `Seurat` or a `SingleCellExperiment` object
#' @param geneID.type the gene ID type, currently support gene symbol and Ensembl
#' @param use.HVG logical, using highly variable genes only or not  
#' @param genesetsize.min minimum number of genes in a gene set, gene set with smaller than this number of genes will be excluded
#' @param genesetsize.max maximum number of genes in a gene set, gene set with greater than this number of genes will be excluded
#' @return a `SeqGeneSet` object
#' @examples 
#' \dontrun{
#' }
#' @author Wei Sheng (wsheng at stu.njmu.edu.cn) 
#' @export
scLoadGS <- function(geneset.file, sc.obj, geneID.type = c("gene.symbol", "ensembl"), use.HVG=FALSE,
                     genesetsize.min = 5, genesetsize.max = 1000)  {
  
  geneID.type <- match.arg(geneID.type, c("gene.symbol", "ensembl"))
  sc.obj.class <- class(sc.obj)
  stopifnot(sc.obj.class %in% c("SingleCellExperiment", "Seurat"))
  
  if(use.HVG) {
    if((sc.obj.class %in% "Seurat") && (length(VariableFeatures(sc.obj))!=0)) {
      genes <- intersect(VariableFeatures(sc.obj), rownames(sc.obj@assays$scGSEA@data))       
    } else {
      warning("HVG information is not available in sce obj, use.HVG is set FALSE!")
    }
  }
  else { 
    if(sc.obj.class %in% "Seurat") { 
      genes <- rownames(sc.obj@assays$scGSEA@data)
    } else {
      genes <- rownames(metadata(sc.obj)[["scGSEA.scale.data"]])
    }
  }
  gs <- loadGenesets(geneset.file, genes, geneID.type=geneID.type, singleCell = TRUE, 
                   genesetsize.min=genesetsize.min, genesetsize.max=genesetsize.max)
  stopifnot(length(gs@GSNames) > 0)
  gs

}

#' Calculate ES score
#' 
#' This function is used to calcuate ES score for each cell and each gene set 
#' 
#' @param sc.obj either a `Seurat` or a `SingleCellExperiment` object
#' @param gs.obj a `SeqGeneSet` object
#' @param weighted.type the weighted type in GSEA
#' @return a `SeqGeneSet` object with data added to slot sc.ES 
#' @examples 
#' \dontrun{
#' }
#' @author Wei Sheng (wsheng at stu.njmu.edu.cn) 
#' @export
scCalES <- function(sc.obj, gs.obj, weighted.type = 0) {
  if( is(sc.obj, "Seurat") ) {
    ES <- log10(1+apply(sc.obj@assays$scGSEA@data[gs.obj@geneList, ], 2, function(x) 
      apply(calES(gs.obj, x, weighted.type = weighted.type), 1, max)))
  }
  if( is(sc.obj, "SingleCellExperiment") ) {
    ES <- log10(1+apply(metadata(sc.obj)[["scGSEA.scale.data"]][gs.obj@geneList, ], 2, function(x) 
      apply(calES(gs.obj, x, weighted.type = weighted.type), 1, max)))
  }
  if(length(gs.obj@GSNames) == 1) { 
    ES <- matrix(ES, nrow=1)
    rownames(ES) <- gs.obj@GSNames
  } else { 
    ES <- data.frame(ES, row.names = gs.obj@GSNames)
  }
  if( is(sc.obj, "Seurat") ) {
    colnames(ES) <- colnames(sc.obj@assays$scGSEA@data)
  } 
  if( is(sc.obj, "SingleCellExperiment") ) {
    colnames(ES) <- colnames(metadata(sc.obj)[["scGSEA.scale.data"]])
  } 
  gs.obj@sc.ES <- as.matrix(ES)
  return(gs.obj)
}

#' ES score normalizaton 
#' 
#' This function calculates ES scores for a permutated data sets, 
#' which is then used for normalization of the observed ES scores 
#' 
#' @param sc.obj either a `Seurat` or a `SingleCellExperiment` object
#' @param gs.obj a `SeqGeneSet` object
#' @param weighted.type the weighted type in GSEA
#' @param perm.time the number of permutations, with default value the number of cells 
#' @param seed the seed for random sampling  
#' @return a `SeqGeneSet` object with data added to slot sc.ES.perm and sc.normFlag reset as TRUE
#' @examples 
#' \dontrun{
#' }
#' @author Wei Sheng (wsheng at stu.njmu.edu.cn) 
#' @export
scESnorm <- function(sc.obj, gs.obj, weighted.type = 0, perm.time=ncol(sc.obj@scale.data), seed=1234) {
  if(gs.obj@sc.normFlag) {return(gs.obj)}
  if(nrow(gs.obj@sc.ES) == 0) { 
    scCalES(sc.obj, gs.obj, weighted.type) 
  }
  ### calculate ES on shuffled dataset
  set.seed(seed)
  if( is(sc.obj, "Seurat") ) {
    gene.exp.shuffled <- t(apply(sc.obj@assays$scGSEA@data[gs.obj@geneList, ], 1, function(x) sample(x)) )
  }
  if( is(sc.obj, "SingleCellExperiment") ) {
    gene.exp.shuffled <- t(apply(metadata(sc.obj)[["scGSEA.scale.data"]][gs.obj@geneList, ], 1, function(x) sample(x)) )
  }
  ES.shuffled <- log10(1+apply(gene.exp.shuffled, 2, function(x) 
    apply(calES(gs.obj, x, weighted.type = weighted.type), 1, max)))
  if(length(gs.obj@GSNames) == 1) {
    ES.shuffled <- matrix(ES.shuffled, nrow=1)
    rownames(ES.shuffled) <- gs.obj@GSNames
  } else {
    ES.shuffled <- data.frame(ES.shuffled, row.names = gs.obj@GSNames)
  }
  ES.shuffled.mean <- apply(ES.shuffled, 1, mean)
  ### normalize ES and ES.shuffled
  NES <- gs.obj@sc.ES / ES.shuffled.mean
  NES.shuffled <- ES.shuffled / ES.shuffled.mean
  gs.obj@sc.ES <- as.matrix(NES)
  gs.obj@sc.ES.perm <- as.matrix(NES.shuffled)
  gs.obj@sc.normFlag = TRUE
  return(gs.obj)
}

#' Add scGESA results to single cell object
#' 
#' This functions adds scGESA results to single cell object
#' 
#' @return a `SeqGeneSet` object with data added to slot sc.ES.perm and sc.normFlag reset as TRUE
#' @examples 
#' \dontrun{
#' }
#' @author Wei Sheng (wsheng at stu.njmu.edu.cn) 
#' @export
add.scGESA.to.obj <- function(sc.obj, gs.obj) {
  stopifnot(gs.obj@sc.normFlag)
  if( is(sc.obj, "Seurat") ) {
    stopifnot(all(colnames(gs.obj@sc.ES) == colnames(sc.obj@assays$scGSEA@data)))
    stopifnot(all(colnames(gs.obj@sc.ES) == rownames(sc.obj@meta.data)))
    #sc.obj <- SetDimReduction(object=sc.obj, reduction.type = 'scGSEA', slot = "cell.embeddings", new.data = as.matrix(t(gs.obj@sc.ES)))
    #sc.obj <- SetDimReduction(object=sc.obj, reduction.type = 'scGSEA', slot = "key", new.data = "scGSEA_")
    #scGSEA <- CreateDimReducObject(embeddings = as.matrix(t(gs.obj@sc.ES)), key="scGSEA_", assay = "RNA")
    sc.obj@meta.data <- data.frame(sc.obj@meta.data, t(gs.obj@sc.ES))
  }
  if( is(sc.obj, "SingleCellExperiment") ) {
    stopifnot(all(colnames(gs.obj@sc.ES) == colnames(counts(sc.obj))))
    stopifnot(all(colnames(gs.obj@sc.ES) == rownames(colData(sc.obj))))
    reducedDim(sc.obj, "scGSEA") <- t(gs.obj@sc.ES)
    colData(sc.obj) <- DataFrame(colData(sc.obj), t(gs.obj@sc.ES))
  }
  sc.obj
}


#' Calculate ES significance
#' 
#' This function calculates ES significance, including p-value, FDR, and FWER
#' 
#' @return a `SeqGeneSet` object with data added to slots sc.pval, sc.FDR, and sc.FWER
#' @examples 
#' \dontrun{
#' }
#' @author Wei Sheng (wsheng at stu.njmu.edu.cn) 
#' @export
scCalSignif <- function(gs.obj) {
  stopifnot(gs.obj@sc.normFlag)
  NES <- gs.obj@sc.ES
  NES.shuffled <- gs.obj@sc.ES.perm
  NES.shuffled.gs_max <- apply(NES.shuffled, 2, max)
  shuffle.time <- ncol(NES.shuffled)
  NES.pval <- apply(NES, 2, function(NES_v) sapply(1:length(NES_v), function(i) sum(NES_v[i] <= NES.shuffled[i,]) / shuffle.time ))
  NES.FWER <- apply(NES, 2, function(NES_v) sapply(NES_v, function(NES_i) sum(NES_i <= NES.shuffled.gs_max) / shuffle.time ))
  NES.FDR <- apply(NES, 2, function(NES_v) sapply(NES_v, function(NES_i) {
    sum( apply(NES_i <= NES.shuffled, 2, sum) ) / sum( apply(NES_i <= NES, 2, sum) ) }))
  NES.FDR[NES.FDR>1] <- 1
  gs.obj@sc.pval <- as.matrix(NES.pval)
  gs.obj@sc.FDR <- as.matrix(NES.FDR)
  gs.obj@sc.FWER <- as.matrix(NES.FWER)
  return(gs.obj)
}


