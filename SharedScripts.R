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
  out<-c(out,paste0("```{r ",ChunkName," ",type,",warning=FALSE}\n"))
  #out<-c(out,"if(length(LocalResults$Gene)>0){\n")
  out<-c(out,paste0("LocalGO <- enrichGO(gene = LocalResults[,c('gene')],
                    universe = ",Universe,",
                    keyType = 'SYMBOL',
                    OrgDb = ",OrgDb,",
                    ont = '",type,"',
                    pAdjustMethod = 'bonferroni',
                    pvalueCutoff  = 0.01,
                    qvalueCutoff  = 0.05,
                    readable = TRUE)\n"))
  out<-c(out,"if(!is.null(LocalGO)){\n")
  out<-c(out,"if(sum(LocalGO@result$qvalue<0.05 & LocalGO@result$p.adjust<0.01)>5){capture.output(LocalPW<-pairwise_termsim(LocalGO));if(length(LocalPW@termsim)>2){capture.output(T<-treeplot(LocalPW));print(T)}}\n")
  out<-c(out,"LocalGO<-LocalGO@result\n")
  out<-c(out,"LocalGO<-LocalGO[LocalGO[,'p.adjust']<0.05,]\n")
  out<-c(out,"LocalGO[,'pvalue']<-formatC(LocalGO[,'pvalue'], format = 'e', digits = 2)\n")
  out<-c(out,"LocalGO[,'p.adjust']<-formatC(LocalGO[,'p.adjust'], format = 'e', digits = 2)\n")
  out<-c(out,"LocalGO[!is.na(LocalGO[,'qvalue']),'qvalue']<-formatC(LocalGO[!is.na(LocalGO[,'qvalue']),'qvalue'], format = 'e', digits = 2)\n")
  out<-c(out,"LocalGO[,'geneID']<-gsub('/',' ', LocalGO[,'geneID'])\n")
  out<-c(out,"LocalGO<-LocalGO[,c('ID', 'Description', 'Count', 'GeneRatio', 'BgRatio', 'pvalue', 'p.adjust', 'qvalue', 'geneID')]\n")
  out<-c(out,"datatable(LocalGO, extensions = c('Buttons','Responsive'), filter=list(position='top'), options=list(dom='Bfrtip', buttons=c('copy','csv','excel')),rownames=FALSE)\n")
  out<-c(out,"}\n")
  #out<-c(out,"}\n")
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
  out<-c(out,paste0(paste0(rep(":",Levels+2),collapse = "")," {.panel-tabset}\n\n"))  
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
  out<-c(out,paste0(paste0(rep(":",Levels+2),collapse = ""),"\n\n"))  
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
    Universe="rownames(merged_seurat)",
    SeuratObject="merged_seurat",
    Identity="seurat_clusters",
    OrgDb="org.Mm.eg.db",
    Levels=3){
  if(Universe=="rownames(merged_seurat)" & SeuratObject!="merged_seurat") {Universe=paste0("rownames(",SeuratObject,")")}
  out<-c()
  out<-c(out,paste0(paste0(rep(":",Levels+1),collapse = "")," {.panel-tabset}\n\n"))  
  out<-c(out,paste0(paste0(rep("#",Levels),collapse = "")," Differentially expressed\n\n"))
  out<-c(out,paste0("```{r ",ChunkName," DEG}\n"))
  out<-c(out,"if(!'package:clusterProfiler' %in% search()){library(clusterProfiler)}\n")
  out<-c(out,"if(!'package:enrichplot' %in% search()){library(enrichplot)}\n")
  out<-c(out,"if(!'package:utils' %in% search()){library(utils)}\n")
  out<-c(out,paste0("LocalResults<-",ResultsDF,"\n"))
  out<-c(out,"if(sum(LocalResults$p_val_adj<0.05)>1){LocalResults<-LocalResults[LocalResults$p_val_adj<0.05, ]}\n")
  out<-c(out,"LocalResults<-LocalResults[order(LocalResults$avg_log2FC,decreasing=TRUE), ]\n")
  out<-c(out,"LocalMarker<-LocalResults\n")
  out<-c(out,"LocalMarker<-LocalMarker[order(LocalResults$avg_log2FC,decreasing=TRUE),]\n")
  out<-c(out,"if(sum(LocalMarker$pct.1>0.5)>1){LocalMarker<-LocalMarker[LocalMarker$pct.1>0.5,]}\n")
  out<-c(out,"LocalMarker<-rownames(LocalMarker)[1]\n")
  out<-c(out,"LocalMarker<-gsub('\\\\.[0-9]*$','',LocalMarker)\n")
  out<-c(out,paste0("Ident(",SeuratObject,")<-",Identity,"\n"))
  out<-c(out,paste0("DimPlot(",SeuratObject,",label = TRUE)+FeaturePlot(",SeuratObject,",features = LocalMarker, cols = c('black', 'gold'))\n"))
  out<-c(out,"LocalResults[,'avg_log2FC']<-formatC(LocalResults[,'avg_log2FC'], format = 'e', digits = 2)\n")
  out<-c(out,"LocalResults[,'p_val']<-formatC(LocalResults[,'p_val'], format = 'e', digits = 2)\n")
  out<-c(out,"LocalResults[,'p_val_adj']<-formatC(LocalResults[,'p_val_adj'], format = 'e', digits = 2)\n")
  out<-c(out,"DT::datatable(LocalResults[,c('gene','pct.1','pct.2','avg_log2FC','p_val','p_val_adj')],colnames = c('Gene','pct in cluster','pct out of cluster','avg log2FC','p val','p val adj'),rownames = FALSE, extensions = 'Buttons', options = list( dom = 'Bfrtip', buttons = c('copy', 'csv', 'excel')),escape=FALSE)\n")
  out<-c(out,"```\n")
  out<-c(out,GenerateAllGOsCode(
    ChunkName=paste(ChunkName,"DEG"),
    Universe=Universe,
    OrgDb=OrgDb,
    Levels=Levels+1))
  
  
  out<-c(out,paste0(paste0(rep("#",Levels),collapse = "")," UpRegulated\n\n"))
  out<-c(out,paste0("```{r ",ChunkName," UP}\n"))
  out<-c(out,paste0("LocalResults<-",ResultsDF,"\n"))
  out<-c(out,"LocalResults<-LocalResults[LocalResults$avg_log2FC>0, ]\n")
  out<-c(out,"if(sum(LocalResults$p_val_adj<0.05)>1){LocalResults<-LocalResults[LocalResults$p_val_adj<0.05, ]}\n")
  out<-c(out,"LocalResults<-LocalResults[order(LocalResults$avg_log2FC,decreasing=TRUE), ]\n")
  out<-c(out,"LocalMarker<-LocalResults\n")
  out<-c(out,"if(sum(LocalMarker$pct.1>0.5)>1){LocalMarker<-LocalMarker[LocalMarker$pct.1>0.5,]}\n")
  out<-c(out,"LocalMarker<-rownames(LocalMarker)[1]\n")
  out<-c(out,paste0("Ident(",SeuratObject,")<-",Identity,"\n"))
  out<-c(out,paste0("DimPlot(",SeuratObject,",label = TRUE)+FeaturePlot(",SeuratObject,",features = LocalMarker, cols = c('black', 'gold'))\n"))
  out<-c(out,"LocalResults[,'avg_log2FC']<-formatC(LocalResults[,'avg_log2FC'], format = 'e', digits = 2)\n")
  out<-c(out,"LocalResults[,'p_val']<-formatC(LocalResults[,'p_val'], format = 'e', digits = 2)\n")
  out<-c(out,"LocalResults[,'p_val_adj']<-formatC(LocalResults[,'p_val_adj'], format = 'e', digits = 2)\n")
  out<-c(out,"DT::datatable(LocalResults[,c('gene','pct.1','pct.2','avg_log2FC','p_val','p_val_adj')],colnames = c('Gene','pct in cluster','pct out of cluster','avg log2FC','p val','p val adj'),rownames = FALSE, extensions = 'Buttons', options = list( dom = 'Bfrtip', buttons = c('copy', 'csv', 'excel')),escape=FALSE)\n")
  out<-c(out,"```\n")
  out<-c(out,GenerateAllGOsCode(
    ChunkName=paste(ChunkName,"UP"),
    Universe=Universe,
    OrgDb=OrgDb,
    Levels=Levels+1))
  
  out<-c(out,paste0(paste0(rep("#",Levels),collapse = "")," DownRegulated\n\n"))
  out<-c(out,paste0("```{r ",ChunkName," DOWN}\n"))
  out<-c(out,paste0("LocalResults<-",ResultsDF,"\n"))
  out<-c(out,"if(!'gene' %in% colnames(LocalResults)){LocalResults$gene<-rownames(LocalResults)}\n")
  out<-c(out,"if(sum(LocalResults$avg_log2FC<0)>1){LocalResults<-LocalResults[LocalResults$avg_log2FC<0, ]}\n")
  out<-c(out,"if(sum(LocalResults$p_val_adj<0.05)>1){LocalResults<-LocalResults[LocalResults$p_val_adj<0.05, ]}\n")
  out<-c(out,"LocalResults<-LocalResults[order(LocalResults$avg_log2FC,decreasing=FALSE), ]\n")
  out<-c(out,"LocalMarker<-LocalResults\n")
  out<-c(out,"if(sum(LocalMarker$pct.2>0.5)>1){LocalMarker<-LocalMarker[LocalMarker$pct.2>0.5,]}\n")
  out<-c(out,"LocalMarker<-rownames(LocalMarker)[1]\n")
  out<-c(out,paste0("Ident(",SeuratObject,")<-",Identity,"\n"))
  out<-c(out,paste0("DimPlot(",SeuratObject,",label = TRUE)+FeaturePlot(",SeuratObject,",features = LocalMarker, cols = c('black', 'gold'))\n"))
  out<-c(out,"LocalResults[,'avg_log2FC']<-formatC(LocalResults[,'avg_log2FC'], format = 'e', digits = 2)\n")
  out<-c(out,"LocalResults[,'p_val']<-formatC(LocalResults[,'p_val'], format = 'e', digits = 2)\n")
  out<-c(out,"LocalResults[,'p_val_adj']<-formatC(LocalResults[,'p_val_adj'], format = 'e', digits = 2)\n")
  out<-c(out,"DT::datatable(LocalResults[,c('gene','pct.1','pct.2','avg_log2FC','p_val','p_val_adj')],colnames = c('Gene','pct in cluster','pct out of cluster','avg log2FC','p val','p val adj'),rownames = FALSE, extensions = 'Buttons', options = list( dom = 'Bfrtip', buttons = c('copy', 'csv', 'excel')),escape=FALSE)\n")
  out<-c(out,"```\n")
  out<-c(out,GenerateAllGOsCode(
    ChunkName=paste(ChunkName,"DOWN"),
    Universe=Universe,
    OrgDb=OrgDb,
    Levels=Levels+1))
  out<-c(out,paste0(paste0(rep(":",Levels+1),collapse = ""),"\n\n"))  
  out<-paste(out,collapse = "",sep="")
  return(out)
}

