#' Title
#'
#' @param ChunkName 
#' @param type 
#' @param Universe 
#' @param OrgDb 
#' @param Levels 
#'
#' @returns
#' @export
#'
#' @examples
GenerateLocalGoCode<-function(
    ChunkName="Test",
    type=c("BP","CC","MF"),
    Universe="rownames(merged_seurat)",
    OrgDb="org.Mm.eg.db",
    Levels=3){
  type <- match.arg(type)
  out<-c()
  if(type=="BP"){out<-c(out,paste0(paste0(rep("#",Levels),collapse = "")," Biological Process\n\n"))}
  if(type=="CC"){out<-c(out,paste0(paste0(rep("#",Levels),collapse = "")," Cellular Component\n\n"))}
  if(type=="MF"){out<-c(out,paste0(paste0(rep("#",Levels),collapse = "")," Molecular Function\n\n"))}
  out<-c(out,paste0("```{r ",ChunkName," ",type,"}\n"))
  out<-c(out,"if(length(LocalResults$Gene)>0){\n")
  out<-c(out,paste0("LocalGO <- enrichGO(gene = LocalResults[,c('gene')],
                    universe = ",Universe,",
                    keyType = 'SYMBOL',
                    OrgDb = ",OrgDb,",
                    ont = '",type,"',
                    pAdjustMethod = 'bonferroni',
                    pvalueCutoff  = 0.01,
                    qvalueCutoff  = 0.05,
                    readable = TRUE)\n"))
  out<-c(out,"if(!is.null(LocalGO)){LocalGO<-LocalGO@results\n")
  out<-c(out,"LocalGO<-LocalGO[LocalGO[,'p.adjust']<0.05,]\n")
  out<-c(out,"LocalGO[,'pvalue']<-formatC(LocalGO[,'pvalue'], format = 'e', digits = 2)\n")
  out<-c(out,"LocalGO[,'p.adjust']<-formatC(LocalGO[,'p.adjust'], format = 'e', digits = 2)\n")
  out<-c(out,"LocalGO[!is.na(LocalGO[,'qvalue']),'qvalue']<-formatC(LocalGO[!is.na(LocalGO[,'qvalue']),'qvalue'], format = 'e', digits = 2)\n")
  out<-c(out,"LocalGO[,'geneID']<-gsub('/',' ', LocalGO[,'geneID'])\n")
  out<-c(out,"LocalGO<-LocalGO[,c('ID', 'Description', 'Count', 'GeneRatio', 'BgRatio', 'pvalue', 'p.adjust', 'qvalue', 'geneID')]\n")
  out<-c(out,"datatable(LocalGO, extensions = c('Buttons','Responsive'), filter=list(position='top'), options=list(dom='Bfrtip', buttons=c('copy','csv','excel')),rownames=FALSE)\n")
  out<-c(out,"}}\n")
  out<-c(out,"```\n\n")
  out<-paste(out,collapse = "",sep="")
  return(out)
}

#' Title
#'
#' @param ChunkName 
#' @param Universe 
#' @param OrgDb 
#' @param Levels 
#'
#' @returns
#' @export
#'
#' @examples
GenerateAllGOsCode<-function(
    ChunkName="Test",
    Universe="rownames(merged_seurat)",
    OrgDb="org.Mm.eg.db",
    Levels=3){
  out<-c()
  out<-c(out,paste0(paste0(rep(":",Levels+1),collapse = "")," .panel-tabset\n\n"))  
  out<-c(out,GenerateLocalGoCode(
    ChunkName=ChunkName,
    type="BP",
    Universe=Universe,
    OrgDb=OrgDb,
    Levels=Levels))
  out<-c(out,GenerateLocalGoCode(
    ChunkName=ChunkName,
    type="CC",
    Universe=Universe,
    OrgDb=OrgDb,
    Levels=Levels))
  out<-c(out,GenerateLocalGoCode(
    ChunkName=ChunkName,
    type="MF",
    Universe=Universe,
    OrgDb=OrgDb,
    Levels=Levels))
  out<-c(out,paste0(paste0(rep(":",Levels+1),collapse = ""),"\n\n"))  
  out<-paste(out,collapse = "",sep="")
  return(out)
}

#' Title
#'
#' @param ChunkName 
#' @param ResultsDF 
#' @param SeuratObject 
#' @param OrgDb 
#' @param Levels 
#'
#' @returns
#' @export
#'
#' @examples
GenerateLocalResultsAllGOs<-function(
    ChunkName="Test",
    ResultsDF="merged_seurat.markers[merged_seurat.markers$cluster=='2',]",
    SeuratObject="merged_seurat",
    OrgDb="org.Mm.eg.db",
    Levels=3){
  out<-c()
  out<-c(out,paste0(paste0(rep(":",Levels+1),collapse = "")," .panel-tabset\n\n"))  
  out<-c(out,paste0(paste0(rep("#",Levels),collapse = "")," Differentially expressed\n\n"))
  out<-c(out,paste0("```{",ChunkName," DEG}\n"))
  out<-c(out,paste0("LocalResults<-",ResultsDF,"\n"))
  out<-c(out,"LocalResults<-LocalResults[LocalResults$p_val_adj<0.05, ]\n")
  out<-c(out,"LocalResults<-LocalResults[order(LocalResults$avg_log2FC,decreasing=TRUE), ]\n")
  out<-c(out,"LocalMarker<-LocalResults$gene[1]\n")
  out<-c(out,paste0("DimPlot(",SeuratObject,",label = TRUE)+FeaturePlot(merged_seurat,features = LocalMarker, cols = c('black', 'gold'))\n"))
  out<-c(out,"LocalMarker[,'avg_log2FC']<-formatC(LocalMarker[,'avg_log2FC'], format = 'e', digits = 2)\n")
  out<-c(out,"LocalMarker[,'p_val']<-formatC(LocalMarker[,'p_val'], format = 'e', digits = 2)\n")
  out<-c(out,"LocalMarker[,'p_val_adj']<-formatC(LocalMarker[,'p_val_adj'], format = 'e', digits = 2)\n")
  out<-c(out,"DT::datatable(LocalResults[,c('gene','pct.1','pct.2','avg_log2FC','p_val','p_val_adj')],colnames = c('Gene','pct in cluster','pct out of cluster','avg log2FC','p val','p val adj'),rownames = FALSE, extensions = 'Buttons', options = list( dom = 'Bfrtip', buttons = c('copy', 'csv', 'excel')),escape=FALSE)\n")
  out<-c(out,"```\n")
  out<-c(out,GenerateAllGOsCode(
    ChunkName=paste(ChunkName,"DEG"),
    Universe=paste0("rownames(",SeuratObject,")"),
    OrgDb=OrgDb,
    Levels=Levels+1))
  
  
  out<-c(out,paste0(paste0(rep("#",Levels),collapse = "")," UpRegulated\n\n"))
  out<-c(out,paste0("```{",ChunkName," UP}\n"))
  out<-c(out,paste0("LocalResults<-",ResultsDF,"\n"))
  out<-c(out,"LocalResults<-LocalResults[LocalResults$avg_log2FC>0, ]\n")
  out<-c(out,"LocalResults<-LocalResults[LocalResults$p_val_adj<0.05, ]\n")
  out<-c(out,"LocalResults<-LocalResults[order(LocalResults$avg_log2FC,decreasing=TRUE), ]\n")
  out<-c(out,"LocalMarker<-LocalResults$gene[1]\n")
  out<-c(out,paste0("DimPlot(",SeuratObject,",label = TRUE)+FeaturePlot(merged_seurat,features = LocalMarker, cols = c('black', 'gold'))\n"))
  out<-c(out,"LocalMarker[,'avg_log2FC']<-formatC(LocalMarker[,'avg_log2FC'], format = 'e', digits = 2)\n")
  out<-c(out,"LocalMarker[,'p_val']<-formatC(LocalMarker[,'p_val'], format = 'e', digits = 2)\n")
  out<-c(out,"LocalMarker[,'p_val_adj']<-formatC(LocalMarker[,'p_val_adj'], format = 'e', digits = 2)\n")
  out<-c(out,"DT::datatable(LocalResults[,c('gene','pct.1','pct.2','avg_log2FC','p_val','p_val_adj')],colnames = c('Gene','pct in cluster','pct out of cluster','avg log2FC','p val','p val adj'),rownames = FALSE, extensions = 'Buttons', options = list( dom = 'Bfrtip', buttons = c('copy', 'csv', 'excel')),escape=FALSE)\n")
  out<-c(out,"```\n")
  out<-c(out,GenerateAllGOsCode(
    ChunkName=paste(ChunkName,"UP"),
    Universe=paste0("rownames(",SeuratObject,")"),
    OrgDb=OrgDb,
    Levels=Levels+1))
  
  out<-c(out,paste0(paste0(rep("#",Levels),collapse = "")," UpRegulated\n\n"))
  out<-c(out,paste0("```{",ChunkName," DOWN}\n"))
  out<-c(out,paste0("LocalResults<-",ResultsDF,"\n"))
  out<-c(out,"LocalResults<-LocalResults[LocalResults$avg_log2FC<0, ]\n")
  out<-c(out,"LocalResults<-LocalResults[LocalResults$p_val_adj<0.05, ]\n")
  out<-c(out,"LocalResults<-LocalResults[order(LocalResults$avg_log2FC,decreasing=FALSE), ]\n")
  out<-c(out,"LocalMarker<-LocalResults$gene[1]\n")
  out<-c(out,paste0("DimPlot(",SeuratObject,",label = TRUE)+FeaturePlot(merged_seurat,features = LocalMarker, cols = c('black', 'gold'))\n"))
  out<-c(out,"LocalMarker[,'avg_log2FC']<-formatC(LocalMarker[,'avg_log2FC'], format = 'e', digits = 2)\n")
  out<-c(out,"LocalMarker[,'p_val']<-formatC(LocalMarker[,'p_val'], format = 'e', digits = 2)\n")
  out<-c(out,"LocalMarker[,'p_val_adj']<-formatC(LocalMarker[,'p_val_adj'], format = 'e', digits = 2)\n")
  out<-c(out,"DT::datatable(LocalResults[,c('gene','pct.1','pct.2','avg_log2FC','p_val','p_val_adj')],colnames = c('Gene','pct in cluster','pct out of cluster','avg log2FC','p val','p val adj'),rownames = FALSE, extensions = 'Buttons', options = list( dom = 'Bfrtip', buttons = c('copy', 'csv', 'excel')),escape=FALSE)\n")
  out<-c(out,"```\n")
  out<-c(out,GenerateAllGOsCode(
    ChunkName=paste(ChunkName,"DOWN"),
    Universe=paste0("rownames(",SeuratObject,")"),
    OrgDb=OrgDb,
    Levels=Levels+1))
  out<-c(out,paste0(paste0(rep(":",Levels+1),collapse = ""),"\n\n"))  
  out<-paste(out,collapse = "",sep="")
  return(out)
}

