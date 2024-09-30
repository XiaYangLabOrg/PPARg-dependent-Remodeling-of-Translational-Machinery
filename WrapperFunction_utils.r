library(ggplot2)
library(varhandle)
library(WriteXLS)
library(reshape2)
library(stringr)


convertMouseGeneList <- function(x){
  
  require("biomaRt")
  human = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
  mouse = useMart("ensembl", dataset = "mmusculus_gene_ensembl")
  
  genesV2 = getLDS(attributes = c("mgi_symbol"), filters = "mgi_symbol", values = x , 
                   mart = mouse, attributesL = c("hgnc_symbol"), martL = human, uniqueRows=T)
  
  # Print the first 6 genes found to the screen
  print(head(genesV2))
  return(genesV2)
}


convertDfGeneColumnMouseHuman <- function(df, toSpecies="human", forPathway=FALSE){
  MT_gene_names = c("MT-ND1"="ND1",  "MT-ND2"="ND2",  "MT-CO1"="COX1",  "MT-CO2"="COX2",  "MT-ATP8"="ATP8", 
                    "MT-ATP6"="ATP6", "MT-CO3"="COX3","MT-ND3" ="ND3", "MT-ND4L"="ND4L", "MT-ND4"="ND4",
                    "MT-ND5"="ND5",  "MT-ND6"="ND6",  "MT-CYB"="CYB")
  if(toSpecies=="mouse"){
    convertedtoMouse <- convertMouseGeneList(df$GENE)
    df$MOUSE = convertedtoMouse$MGI.symbol[match(df$GENE, convertedtoMouse$HGNC.symbol)]
    df$MOUSE[which(is.na(df$MOUSE))] <- tolower(df$GENE[which(is.na(df$MOUSE))])
    df$GENE <- df$MOUSE
    df$GENE <- paste0(toupper(x=substr(df$GENE, start = 1, stop = 1)),tolower(substring(df$GENE, first = 2)))
    # correct for mt- genes
    MT_mouse_gene_names = MT_gene_names
    new_names = c()
    for(name in names(MT_gene_names)){
      new_names = append(new_names, 
                         paste0(tolower(unlist(strsplit(name, split = "-"))[1]), "-",
                                paste0(toupper(x=substr(unlist(strsplit(name, split = "-"))[2], start = 1, stop = 1)),
                                       tolower(substring(unlist(strsplit(name, split = "-"))[2], first = 2)))))
    }
    names(MT_mouse_gene_names) = new_names
    
    for(gene in 1:nrow(df)){
      if(sum(MT_mouse_gene_names==df$gene[gene])>0){
        df$GENE[gene] = names(MT_mouse_gene_names)[MT_mouse_gene_names==df$gene[gene]]
      }
      else{
        next
      }
    }
    return(df)
  }
  else if(toSpecies=="human"){
    convertedToHuman <- convertMouseGeneList(df$GENE)
    df$HUMAN <- convertedToHuman$HGNC.symbol[match(df$GENE, convertedToHuman$MGI.symbol)]
    df$HUMAN[which(is.na(df$HUMAN))] <- toupper(df$GENE[which(is.na(df$HUMAN))])
    df$GENE <- df$HUMAN
    df$HUMAN <- NULL
    if(sum(is.na(df$GENE))>0) cat("Warning: NAs created.\n")
    if(forPathway){ # convert MT- genes
      for(gene in 1:nrow(df)){
        if(grepl("^MT-",df$GENE[gene])){
          df$GENE[gene] = MT_gene_names[df$GENE[gene]]
        }
        else{
          next
        }
      }
    }
    return(df)
  }
  else cat("Put a valid toSpecies value - 'human' or 'mouse'\n")
}


addInfoAndTrim <- function(DEG_df, 
                           p_cutoff=0.01, 
                           multiple_comparisons=TRUE, 
                           convertToHuman=TRUE, 
                           lfc_cutoff=0.1){
  # filter DEG in DEG_df by p_val cutoff
  DEG_df = DEG_df[DEG_df$p_val<p_cutoff,]
  # filter DEG in DEG_df by avg_logFC cutoff (absolute value)
  DEG_df = DEG_df[abs(DEG_df$avg_logFC)>lfc_cutoff,]
  module = c()
  module_ct = c()
  module_comp = c()
  for(row in 1:nrow(DEG_df)){
    if(DEG_df$avg_logFC[row]>0){
      if(multiple_comparisons){
        module[row] = paste(DEG_df$Comparison[row], DEG_df$Cell_type[row], "UP", sep = ".")
        module_ct[row] = paste(DEG_df$Cell_type[row], "UP", sep = ".")
        module_comp[row] = paste(DEG_df$Comparison[row], "UP", sep = ".")
      }
      else{
        module[row] = paste(DEG_df$Cell_type[row], "UP", sep = "_")
      }
    }
    else{
      if(multiple_comparisons){
        module[row] = paste(DEG_df$Comparison[row], DEG_df$Cell_type[row], "DOWN", sep = ".")
        module_ct[row] = paste(DEG_df$Cell_type[row], "DOWN", sep = ".")
        module_comp[row] = paste(DEG_df$Comparison[row], "DOWN", sep = ".")
      }
      else{
        module[row] = paste(DEG_df$Cell_type[row], "DOWN", sep = "_")
      }
    }
  }
  DEG_df$MODULE = module
  DEG_df$MODULE_ct = module_ct
  DEG_df$MODULE_comp = module_comp
  DEG_df$MODULE_no_direct = paste(DEG_df$Cell_type, DEG_df$Comparison)
  DEG_df = DEG_df[order(DEG_df$MODULE),]
  
  # if(convertToHuman){
  #  temp <- convertDfGeneColumnMouseHuman(df = DEG_df, toSpecies = "human", forPathway = FALSE)
  #    DEG_df$HUMAN = temp$GENE
  # }
  DEG_df$HUMAN = toupper(DEG_df$GENE)
  return(DEG_df)
}


pathway_enrichment <- function(deg_list, 
                               FDR_threshold=NULL, 
                               pval_threshold=0.01, 
                               celltype_cluster,
                               identifier, 
                               pVal=TRUE){
  
  # pVal this is our method of switching off between FDR and pval threshold
  list_of_pathway_databases = list.files(path = "./resource/", pattern = "*.txt", full.names = TRUE)
  
  # only take specific module
  #deg_list = deg_list[deg_list$MODULE==celltype_cluster,]
  deg_list[["gene"]] <- deg_list$GENE
  deg_list[["gene"]] <- toupper(deg_list[["gene"]]) #changes the gene names to upper case for compatibility with database (redundant)
  deg_list[["module"]] <- rep(celltype_cluster, nrow(deg_list))
  deg_list <- deg_list[,c("module","gene")]
  unique_module1 <- unique(deg_list$module)
  module1_len <- length(unique(deg_list$module))
  x <- list()
  for(z in 1:length(list_of_pathway_databases)){
    pathway_database <- list_of_pathway_databases[z]
    database_name = unlist(strsplit(pathway_database,"/"))
    database_name = database_name[length(database_name)]
    database_name = unlist(strsplit(database_name, ".txt"))[1]
    print(database_name)
    
    Module2 <- tool.read(pathway_database)
    Unique_module2 <- unique(Module2$module)
    
    Module2_len <- length(unique(Module2$module))
    
    data_matrix_for_enrichment <- data.frame()
    List_initial <- 1
    # go through the different databases
    for(k in 1:Module2_len){
      data_matrix_for_enrichment[List_initial,1] <- unique_module1
      data_matrix_for_enrichment[List_initial,2] <- length(deg_list$gene[which(match(deg_list$module, unique_module1)>0)])
      
      Overlapped_genes <- intersect(deg_list$gene[which(match(deg_list$module, unique_module1)>0)],
                                    Module2$gene[which(match(Module2$module, Unique_module2[k])>0)])
      data_matrix_for_enrichment[List_initial,3] <- length(Overlapped_genes)
      data_matrix_for_enrichment[List_initial,4] <- length(Module2$gene[which(match(Module2$module, Unique_module2[k])>0)])
      
      data_matrix_for_enrichment[List_initial, 5] <- 20000
      data_matrix_for_enrichment[List_initial,6] <- Unique_module2[k]
      
      if(length(Overlapped_genes)){
        data_matrix_for_enrichment[List_initial,7] <- paste(Overlapped_genes, collapse = ",")
      } else{
        data_matrix_for_enrichment[List_initial,7] <- c("NULL")
      }
      data_matrix_for_enrichment[List_initial,8] <- database_name
      
      List_initial=List_initial + 1
    }
    ifelse(!dir.exists(file.path(paste0("./Csvs","_",identifier,"/pathway_enrichment"))), dir.create(file.path(paste0("./Csvs","_",identifier,"/pathway_enrichment")),recursive = T), FALSE)
    write.table(data_matrix_for_enrichment, file = paste0("./Csvs","_",identifier,"/pathway_enrichment/temp1.",celltype_cluster,".dat"), quote = FALSE, sep = "\t",
                row.names = FALSE, col.names = FALSE)
    
    record_mat <- read.table(paste0("./Csvs","_",identifier,"/pathway_enrichment/temp1.",celltype_cluster,".dat"))
    record_length <- dim(record_mat)
    enrichment_score <- data.frame()
    
    for(i in 1:record_length[1]){
      enrichment_score[i,1] <- record_mat[i,1]
      enrichment_score[i,2] <- phyper(record_mat[i,3], record_mat[i,4], record_mat[i,5]-record_mat[i,4], record_mat[i,2], lower.tail = FALSE)
      enrichment_score[i,3] <- record_mat[i,3]/record_mat[i,2]*record_mat[i,5]/record_mat[i,4]
      enrichment_score[i,4] <- 0
      enrichment_score[i,5]<-record_mat[i,8]
      enrichment_score[i,6]<-record_mat[i,6]
      enrichment_score[i,7]<-record_mat[i,3]
      enrichment_score[i,8]<-record_mat[i,7]
    }
    enrichment_score[,4]<-p.adjust(enrichment_score[,2], 'bonferroni')
    colnames(enrichment_score) <- c("Module","Pval","Enrichment","FDR","PathwaySource","Pathway","nOverlap","Overlap")
    enrichment_score <- enrichment_score[order(enrichment_score$FDR),]
    x[[database_name]] <- data.frame(enrichment_score)
    
    if(z==1){
      all_pathways_df <- enrichment_score
    }else{
      all_pathways_df <- rbind(all_pathways_df, enrichment_score)
    }
  }
  
  all_pathways_df <- all_pathways_df[order(all_pathways_df$FDR),]
  all_pathways_df <- all_pathways_df[which(all_pathways_df$FDR < 0.05),]
  # all_pathways_df <- all_pathways_df[which(all_pathways_df$Pval < 0.01),]
  all_pathways_df <- all_pathways_df[which(all_pathways_df$nOverlap > 3),]
  x[["Combined"]] <- data.frame(all_pathways_df)
  x[["Combined"]][["Disease_Model"]] <- rep(celltype_cluster, nrow(x[["Combined"]]))
  
  total = data.frame(stringsAsFactors = FALSE)
  total = rbind(x[[1]], x[[2]])
  total = rbind(total, x[[3]])
  total = rbind(total, x[[4]])
  total$ModuleSize = nrow(deg_list)
  ifelse(!dir.exists(file.path(paste0("./consolidated_pathways","_",identifier))), dir.create(file.path(paste0("./consolidated_pathways","_",identifier)),recursive = T), FALSE)
  write.table(total, paste0("./consolidated_pathways","_",identifier,"/", celltype_cluster, ".txt"), row.names = FALSE, quote = FALSE, sep = "\t")
  return(total)
  
}


tool.read <- function(file, vars=NULL) {
  if(is.null(file)) return(data.frame())
  if(file == "") return(data.frame())
  dat <- read.delim(file=file, header=TRUE,
                    na.strings=c("NA", "NULL", "null", ""),
                    colClasses="character", comment.char="",
                    stringsAsFactors=FALSE)
  if(is.null(vars) == FALSE) dat <- dat[,vars]
  dat <- na.omit(dat)
  return(dat)
}



# Function to transform pathway names
transform_pathway_name <- function(pathway) {
  # Remove the first term followed by an underscore
  pathway <- str_replace(pathway, "^[^_]+_", "")
  
  # Replace underscores with spaces and convert to lowercase
  pathway <- str_replace_all(pathway, "_", " ")
  pathway <- tolower(pathway)
  
  # Capitalize the first letter of the entire string
  pathway <- paste(toupper(substr(pathway, 1, 1)), substr(pathway, 2, nchar(pathway)), sep = "")
  
  return(pathway)
}
