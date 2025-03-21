library(shiny)
library(AnnotationDbi)
library(ComplexHeatmap)
library(InteractiveComplexHeatmap)
library(org.Mm.eg.db)
library(dplyr)
library(tidyr)
library(tibble)
library(data.table)
library(circlize)
library(shinyWidgets)
library(paletteer)
library(dendsort)

# randomids<-sample(unique(data2plot2$GENE),size =100 )
# randomgenes<-sample(unique(data2plot2$SYMBOL),size =100 )
# 
# tcounts2<-fread("Salmon_Gene_table_TPM.csv")
# 
# res <- AnnotationDbi::select(org.Mm.eg.db,
#                              keys = tcounts2$V1,
#                              keytype = "ENSEMBL",
#                              columns = c("ENSEMBL", "SYMBOL"))
# res$SYMBOL[is.na(res$SYMBOL)]<-res$ENSEMBL[is.na(res$SYMBOL)]
# 
# 
# 
# data2plot2<-tcounts2 %>% rename("V1"="GENE") %>% pivot_longer(-GENE, names_to = "Sample", values_to = "TPMs") %>%
#   mutate(LPI_annotation=case_when(
#     Sample == "Alpha-2nd_PP-Pro-II" ~ "Alpha-2nd cell/PP-Pro-II",
#     Sample == "EP3" ~ "EP3",
#     Sample == "EP2" ~ "EP2",
#     Sample == "Trunk-Duct-EP1" ~ "Trunk/Duct/EP1",
#     Sample == "Alpha-1st_eraly_and_late" ~ "Alpha-1st eraly & late",
#     Sample == "EP4-early" ~ "EP4-early",
#     Sample == "Trunk" ~ "Trunk",
#     Sample == "MP-late_and_Tip" ~ "MP-late + Tip",
#     Sample == "EP4-late" ~ "EP4-late",
#     Sample == "Beta-early" ~ "Beta-early",
#     Sample == "Alpha-2nd_PP-Pro-I" ~ "Alpha-2nd cell/PP-Pro-I",
#     Sample == "Delta-early" ~ "Delta-early",
#     Sample == "Beta-late" ~ "Beta-late",
#     Sample == "Beta" ~ "Beta",
#     Sample == "Alpha-2nd_cell" ~ "Alpha-2nd cell",
#     Sample == "Acinar" ~ "Acinar",
#     Sample == "PP" ~ "PP",
#     Sample == "Others" ~ "Others",
#     Sample == "Delta-late" ~ "Delta-late",
#     Sample == "Immune" ~ "Immune",
#     Sample == "Pre-alpha-1st" ~ "Pre-alpha-1st cell",
#     Sample == "MP-early" ~ "MP-early"))%>%
#   left_join(res, by=c("GENE"="ENSEMBL"),relationship = "many-to-many") %>%
#   mutate(LPI_annotation=factor(LPI_annotation,levels=orderCT))

#data2plot2 %>% group_by(GENE) %>% mutate(M=mean(TPMs), S=sd(TPMs)) %>% ungroup() %>% mutate(Z.score=ifelse(S==0,0,(TPMs-M)/S)) %>% dplyr::select(GENE,SYMBOL,LPI_annotation,TPMs,Z.score) %>% write.table("TPM_tibble.csv",quote=F,sep=",", row.names = F)



# Load your gene expression data
# Replace this with your actual data loading code

ui = fluidPage(
  titlePanel("Gene Expression Heatmap"),
  sidebarLayout(
    sidebarPanel(
      actionButton("useDemo1", "Demo genes"),
      textAreaInput("gene_list", "Enter gene names (one per line, ENSEMBLID or GeneNames are accepted):", rows = 10),
      actionButton("plot_heatmap", "Generate Heatmap"),
      numericInput("threshold", "TPM threshold for binary heatmap", min=1,1)
    ),
    mainPanel(
      tabsetPanel(
        id="tabset",
        tabPanel("TPMs",InteractiveComplexHeatmapOutput("heatmap_0")),
        tabPanel("Z-score",InteractiveComplexHeatmapOutput("heatmap_1")),
        tabPanel("binary",InteractiveComplexHeatmapOutput("heatmap_2"))
      ))
  )
)



server = function(input, output, session) {
  default_genes <- fread("data/Default_genes.tab", header = F) %>% .$V1
  
  
  observe({
    if (input$useDemo1) {
      updateTextInput(session, "gene_list", value = paste(default_genes,collapse = "\n"))
    }
  })
  
  celltypecolors<-c("MP-early"="#f8ac30","MP-late + Tip"="#f0853c",
                    "Trunk"="#eb5c3f","Trunk/Duct/EP1"="#77f07f","Acinar"="#dfdd19",
                    "Pre-alpha-1st cell"="#0099FF", "Alpha-1st eraly & late"="#0019ff","EP2"="#3ec995","EP3"="#41a0ae",
                    "EP4-early"="#36669c","EP4-late"="#3a2f6b","Alpha-2nd cell/PP-Pro-I"="#70AF63","Epsilon"="#4f75a3",
                    "Alpha-2nd cell/PP-Pro-II"="#3A6032","Alpha-2nd cell"="#4E8243","PP"="#e0322b",
                    "Beta-early"="#A163AF","Beta-late"="#583260","Beta"="#774382",
                    "Delta"="#53a980","Delta-early"="#a9d4bf","Delta-late"="#295440",
                    "Immune"="#b10808","Others"="gray50")
  
  orderCT<-c("MP-early",
             "MP-late + Tip","Trunk","Trunk/Duct/EP1",
             "EP2","EP3","EP4-early","EP4-late","Beta-early","Beta-late","Beta",
             "Alpha-2nd cell/PP-Pro-I","Alpha-2nd cell/PP-Pro-II","Alpha-2nd cell",
             "PP","Epsilon","Delta-early","Delta-late","Delta","Pre-alpha-1st cell",
             "Alpha-1st eraly & late","Acinar","Immune",
             "Others")
  
  data2plot2<-read.table("data/TPM_tibble.csv",sep = ",",header = T)
  colorder<-orderCT[orderCT %in% unique(data2plot2$LPI_annotation)]
  
  col_fun_tpms = colorRamp2(seq(from=0,to=100,length.out = 10), paletteer_d("ggsci::deep_orange_material"))
  col_fun_z = colorRamp2(seq(from=-2,to=2,length.out = 11), paletteer_d("Redmonder::dPBIPuGn"))
  col_fun_bin = colorRamp2(c(0,1), c("white", "gray30"))
  
  observeEvent(input$plot_heatmap, {
    genes2look <- unlist(strsplit(input$gene_list, "\n"))
    temp<-data2plot2 %>% filter(GENE %in% genes2look |SYMBOL %in% genes2look ) %>%
      group_by(SYMBOL, LPI_annotation) %>% 
      reframe(TPM=mean(TPMs), Z.score=median(Z.score))
    mat_raw<-temp%>%dplyr::select(SYMBOL, LPI_annotation,TPM)%>%
      pivot_wider(id_cols = "SYMBOL", values_from = "TPM", names_from = "LPI_annotation") %>%
      column_to_rownames("SYMBOL") %>% as.matrix()
    mat<-temp%>%dplyr::select(SYMBOL, LPI_annotation,Z.score)%>%
      pivot_wider(id_cols = "SYMBOL", values_from = "Z.score", names_from = "LPI_annotation") %>%
      column_to_rownames("SYMBOL") %>% as.matrix()
    
    mat<-mat[,colorder]
    mat_raw<-mat_raw[,colorder]
    
    mat2<-(mat_raw>=input$threshold)*1
    
    
    ha = HeatmapAnnotation("Cell Type" = colorder,col = list("Cell Type"=celltypecolors),
                           show_legend = F)
    
    cortab<-cor(t(mat),method = "pearson")
    cortab[is.na(cortab)]<-1
    roworder<-dendsort(hclust(as.dist(1-cortab)))
    
    hr<-Heatmap(mat_raw, name = "Raw TPMs", col = col_fun_tpms,cluster_columns = F,
                top_annotation = ha,cluster_rows = roworder)
    hr<-draw(hr)
    makeInteractiveComplexHeatmap(input, output, session, hr, "heatmap_0")
    hm<-Heatmap(mat, name = "Z-score", col = col_fun_z,cluster_columns = F,
                top_annotation = ha,cluster_rows = roworder)
    hm<-draw(hm)
    makeInteractiveComplexHeatmap(input, output, session, hm, "heatmap_1")
    hm2<-Heatmap(mat2, name = "Binary TPM", col =c("white","gray30"),cluster_columns = F,
                 top_annotation = ha,cluster_rows = roworder, show_heatmap_legend = F)
    hm2<-draw(hm2)
    makeInteractiveComplexHeatmap(input, output, session, hm2, "heatmap_2")
  })
}

shinyApp(ui, server)
