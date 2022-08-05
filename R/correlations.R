##############################################################################
# Generate Jaccard Partitions                                                #
# Copyright (C) 2022                                                         #
#                                                                            #
# This program is free software: you can redistribute it and/or modify it    #
# under the terms of the GNU General Public License as published by the      #
# Free Software Foundation, either version 3 of the License, or (at your     #
# option) any later version. This program is distributed in the hope that    #
# it will be useful, but WITHOUT ANY WARRANTY; without even the implied      #
# warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the   #
# GNU General Public License for more details.                               #     
#                                                                            #
# Elaine Cecilia Gatto | Prof. Dr. Ricardo Cerri | Prof. Dr. Mauri Ferrandin #
# Federal University of Sao Carlos (UFSCar: https://www2.ufscar.br/) Campus  #
# Sao Carlos Computer Department (DC: https://site.dc.ufscar.br/)            # 
# Program of Post Graduation in Computer Science                             #
# (PPG-CC: http://ppgcc.dc.ufscar.br/) Bioinformatics and Machine Learning   #
# Group (BIOMAL: http://www.biomal.ufscar.br/)                               #
#                                                                            #
##############################################################################


########################################################################
# WORSKSPACE
########################################################################
FolderRoot = "~/jaccard"
FolderScripts = "~/jaccard/R"


############################################################################
# FUNCTION COMPUTE JACCARD                                                 #
#   Objective:                                                             #
#      Modeling correlations with index Jaccard                            #
#   Parameters                                                             #
#       ds: specific dataset information                                   #
#       resLS: specific dataset label space                                #
#       dataset_name: dataset name. It is used to save files.              #
#       namesLabels: label names                                           #
#       FolderHC: hclust and cutree folder path                            #
#       number_folds: number of folds created                              #
#   Return                                                                 #
#       correlation matrix                                                 #
############################################################################
computeJaccard <- function(ds,
                           dataset_name,
                           number_dataset, 
                           number_cores, 
                           number_folds, 
                           folderResults,
                           resLS,
                           namesLabels){
  
  diretorios = directories(dataset_name, folderResults)
  
    s = 1
    jaccardParalel <- foreach(s = 1:number_folds) %dopar%{
      
      cat("\n\nFold: ", s)
      
      #####################################################################
      FolderRoot = "~/jaccard"
      FolderScripts = "~/jaccard/R"
      
      #####################################################################
      # LOAD LIBRARIES
      setwd(FolderScripts)
      source("libraries.R")
      
      setwd(FolderScripts)
      source("utils.R")
      
      #####################################################################
      merge_matrix <- function(matrix_correlation){
        #cat("\n\tMerge matrix")
        matrix_correlation <- round(matrix_correlation,4)
        melt_mat_cor <- melt(matrix_correlation)
        return (melt_mat_cor)
        cat("\n")
        gc()
      }
      
      #####################################################################
      get_lower_tri<-function(matrix_correlation){
        #cat("\n\tGet Lower Tri")
        matrix_correlation_1[upper.tri(matrix_correlation)] <- NA
        return(matrix_correlation_1)
        cat("\n")
        gc()
      }
      
      #####################################################################
      get_upper_tri <- function(matrix_correlation){
        #cat("\n\tGet Upper Tri")
        matrix_correlation[lower.tri(matrix_correlation)]<- NA
        return(matrix_correlation)
        cat("\n")
        gc()
      }
      
      #####################################################################
      cut_matrix <- function(matrix_correlation, measure){
        #cat("\n\tCut Matrix")
        upper_tri <- get_upper_tri(matrix_correlation)
        melt_mat_cor <- melt(upper_tri, na.rm = TRUE) 
        return(melt_mat_cor)
        cat("\n")
        gc()
      }
      
      #####################################################################
      reorder_mat_cor <- function(matrix_correlation){
        #cat("\n\tReorder Matrix")
        dd <- as.dist((1-matrix_correlation)/2)
        hc <- hclust(dd)
        print(hc)
        matrix_correlation <- matrix_correlation[hc$order, hc$order]
        return(matrix_correlation)
        cat("\n")
        gc()
      }
      
      #####################################################################
      # created folder for the split
      FolderHCES = paste(diretorios$folderPartitions, "/Split-", s, sep="")
      if(dir.exists(FolderHCES)==TRUE){
        cat("\n")
      } else{
        dir.create(FolderHCES)  
      }
      setwd(FolderHCES)
      
      #cat("\nGET LABELS SPACE\n")
      classes = resLS$Classes[s]
      classes = data.frame(classes)
      classes = t(classes)
      
      #cat("\nCOMPUTES JACCARD\n")
      matrix_correlation = distance(classes, method = "jaccard", 
                                    use.row.names = TRUE)
      write.csv(matrix_correlation, "matrix_correlation.csv", 
                row.names = FALSE)
      nomesRotulos = colnames(matrix_correlation)
      
      #cat("\nCHECK NA\n")
      matrix_correlation_na = is.na(matrix_correlation)
      matrix_correlation_na_2 = replace(x = matrix_correlation, 
                                        list = is.na(matrix_correlation), 
                                        values = 0)
      
      #cat("\nGET COL NAMES\n")
      rownames(matrix_correlation) <- namesLabels
      matrix_correlation_2 = as.matrix(matrix_correlation)
      
      #cat("\nREORGANIZE\n")
      matrix_correlation_order <- reorder_mat_cor(matrix_correlation_2)
      upper_tri <- get_upper_tri(matrix_correlation_order)
      
      #cat("\nMELT MATRIX\n")
      melt_mat_cor <- melt(upper_tri, na.rm = TRUE)
      write.csv(melt_mat_cor, "melt_mat_cor.csv", row.names = FALSE)
    
      gc()
    }
  
  
  gc()
  cat("\n##################################################################")
  cat("\n# END OF COMPUTE JACCARD FUNCTION                                #")
  cat("\n##################################################################")
  cat("\n\n\n\n")
}


###########################################################################
# FUNCTION CUTREE HCLUST                                                   
#   Objective                                                              
#       Partitions the correlation matrix using a hierarchical clustering 
#    algorithm 
#   Parameters
#       ds: specific dataset information
#       resLS: specific dataset label space                            
#       dataset_name: dataset name. It is used to save files.          
#       namesLabels: label names                                       
#       FolderHClust: hclust and cutree folder path                    
#       number_folds: number of folds created                          
#   Return                                                             
#       partitions and graphics                                        
###########################################################################
CutreeHClust <- function(ds,
                         dataset_name,
                         number_dataset, 
                         number_cores, 
                         number_folds, 
                         folderResults,
                         namesLabels,
                         resLS){
  
  diretorios = directories(dataset_name, folderResults)
  
  # method	
  # the agglomeration method to be used. This should be
  # (an unambiguous abbreviation of) one of "ward.D", "ward.D2", "single", 
  # "complete", "average" (= UPGMA), "mcquitty" (= WPGMA), "median" (= WPGMC) 
  # or "centroid" (= UPGMC).
  
  s = 1
  cutreeParalel <- foreach(s = 1:number_folds) %dopar%{
    
    cat("\nFold: ", s) 
    
    FolderRoot = "~/jaccard"
    FolderScripts = "~/jaccard/R"
    
    # LOAD LIBRARIES
    setwd(FolderScripts)
    source("libraries.R")
    
    setwd(FolderScripts)
    source("utils.R")
    
    diretorios = directories(dataset_name, folderResults)
    
    #cat("\nCreate DataFrames\n")
    fold = c(0)
    average = c(0)
    complete = c(0)
    mcquitty = c(0)
    single = c(0)
    wardD = c(0)
    wardD2 = c(0)
    c4 = data.frame(fold, average, complete, mcquitty, single, wardD, wardD2)
    
    fold = c(0)
    metodo = c(0)
    coeficiente = c(0)
    coefHC = data.frame(fold, metodo, coeficiente)
    
    # methods
    metodos = c("average", "complete", "mcquitty", "single", "ward.D", "ward.D2")    
    
    
    #####################################################################
    # created folder for the split
    FolderHClust = paste(diretorios$folderPartitions, "/Split-", s, sep="")
    if(dir.exists(FolderHClust)==FALSE){dir.create(FolderHClust)}
    
    #cat("\nOpen matrix correlation\n")
    setwd(FolderHClust)
    matrix_correlation = data.frame(read.csv("matrix_correlation.csv"))
    rownames(matrix_correlation) <- namesLabels
    matrix_correlation_2 = as.matrix(matrix_correlation)    
    
    # for the first method to the last
    i = 1
    for(i in i:length(metodos)){
      
      cat("\n==============================")
      cat("\nMethod: ", metodos[i], "\n")
      
      #cat("\nCreates the folder to save information for this method\n")      
      FolderMethods = paste(FolderHClust, "/", metodos[i], sep="")
      if(dir.exists(FolderMethods)==FALSE){dir.create(FolderMethods)}
      
      #cat("\nCreates the folder to save graphics\n")            
      FolderGraphics = paste(FolderMethods, "/Graphics", sep="")
      if(dir.exists(FolderGraphics)==FALSE){dir.create(FolderGraphics)}
      
      #cat("\nCreates the folder to save clusters\n")                  
      FolderClusters = paste(FolderMethods, "/Clusters", sep="")
      if(dir.exists(FolderClusters)==FALSE){dir.create(FolderClusters)}
      
      #cat("\nDEND\n")                  
      Dend <- matrix_correlation_2 %>% as.dist %>% hclust(method = metodos[i]) %>% as.dendrogram
      
      #cat("\nOTTER DENDRO\n")                  
      OtterDendro = as.dendrogram(hclust(d = as.dist(matrix_correlation_2),
                                         method=metodos[i]))
      
      #cat("\nAsDist = as.dist(matrix_correlation)\n")                  
      AsDist = as.dist(matrix_correlation)
      
      #cat("\nAsDistMatrix = as.matrix(AsDist)\n")                  
      AsDistMatrix = as.matrix(AsDist)
      
      #cat("\nHC = hclust(AsDist, method=metodos[i])\n")                  
      HC = hclust(AsDist, method=metodos[i])
      
      #cat("\nDendro = as.dendrogram(HC)\n")                  
      Dendro = as.dendrogram(HC)
      
      #cat("\nCreates the folder to save clusters\n")                  
      DendData <- dendro_data(Dendro, type = "rectangle")
      
      #cat("\nSAVE COEFF\n")                  
      fold = s
      metodo = metodos[i]
      coeficiente = coef.hclust(HC)
      coefHC = rbind(coefHC, data.frame(fold, metodo, coeficiente))
      #print(coefHC)   
      
      #cat("\nGRAPHIC: RADIAL\n")
      setwd(FolderGraphics)
      pdf("radial.pdf", width = 10, height = 8)
      print(plot(as.phylo(HC), type = "radial", cex = 0.6, no.margin = TRUE))
      dev.off()
      cat("\n")
      
      #cat("\nGRAPHIC: FAN\n")
      pdf("fan.pdf", width = 10, height = 8)
      print(plot(as.phylo(HC), type = "fan", cex = 0.6, no.margin = TRUE))
      dev.off()
      cat("\n")
      
      #cat("\nGRAPHIC: UNROOT\n")      
      pdf("unroot.pdf", width = 10, height = 8)
      print(plot(as.phylo(HC), type = "unrooted", cex = 0.6, 
                 no.margin = TRUE))
      dev.off()
      cat("\n")
      
      #cat("\nGRAPHIC: CLADOGRAM\n")
      pdf("cladogram.pdf", width = 10, height = 8)
      print(plot(as.phylo(HC), type = "cladogram", cex = 0.6, 
                 no.margin = TRUE))
      dev.off()
      cat("\n")
      
      #cat("\nGRAPHIC: DENDRO\n")
      pdf("hc_plot.pdf", width = 10, height = 8)
      print(plot(Dendro))
      print(with(pvclust:::hc2axes(as.hclust(Dendro)), 
                 text(x.axis, y.axis, round(y.axis, 2),col = "red", 
                      adj = c(0.5, 1.5), cex = 0.5)))
      dev.off()
      cat("\n")     
      
      ############################################################### 
      #cat("\nClustering: from the first to the last label\n")                  	  
      clusters3 = data.frame(c(0))
      tables = data.frame()
      
      k = 1
      for(k in 1:ds$Labels){
        
        cat("\n===============")
        cat("\ncluster: ", k)   
        
        setwd(FolderClusters)
        
        #cat("\nCUTREE\n")        
        cutLabels = cutree(HC, k)
        clusters = data.frame(cutree(HC, k))
        names(clusters) = "grupo"
        rotulos = c(rownames(clusters))
        
        #cat("\nSAVE CUTREE\n")
        write.csv(clusters, paste("cluster_", k, ".csv", sep=""))		
        
        clusters3 = cbind(clusters3, clusters$grupo)
        names(clusters3)[k+1] = paste("particao-",k, sep="")
        
        k = k + 1 # increment label
        gc() # garbage collection
        
      } # fim dos clusters
      
      #cat("\nSAVE ALL CUTREE\n")        
      setwd(FolderMethods)
      clusters4 = clusters3[,-1]
      clusters5 = data.frame(rotulos, clusters4)
      write.csv(clusters5, paste(metodos[i], "-partitions.csv", sep=""), 
                row.names = FALSE)
      
      i = i + 1 # increment hclust method
      gc() # garbage collection
      
    } # fim do método
    
    #cat("\nSAVE ALL COEFICIENTES\n")        
    setwd(FolderHClust)
    write.csv(coefHC[-1,], paste("fold-", s, "-coeficientes.csv", sep=""), 
              row.names = FALSE)
    
    gc()
  }
  
  gc()
  cat("\n##################################################################")
  cat("\n# END OF THE FUNCTION HCLUST AND CUTREE                          #")
  cat("\n##################################################################")
  cat("\n\n\n\n")
}


#######################################################################
# FUNCTION BEST COEFFICIENT HCLUST                                 
#   Objective                                                      
#       lists the best HCLUST coefficients found for 10-folds      
#   Parameters                                                     
#       number_folds: number of folds created                      
#       FolderHClust: hclust and cutree folder path                
#   Return                                                         
#       csv with the best coefficients                             
###################################################################
bestCoefficient <- function(ds,
                            dataset_name,
                            number_dataset, 
                            number_cores, 
                            number_folds, 
                            folderResults,
                            resLS,
                            namesLabels){
  
  diretorios = directories(dataset_name, folderResults)
  
  todos = data.frame()
  primeiro = data.frame()
  segundo = data.frame()
  terceiro = data.frame()
  quarto = data.frame()
  quinto = data.frame()
  sexto = data.frame()
  
  # fold 1 to s
  s = 1
  while(s<=number_folds){  
    
    cat("\nFold: ", s)
    
    # enter the fold
    FolderSplit = paste(diretorios$folderPartitions, "/Split-", s, sep="")
    
    # get the coeffs
    setwd(FolderSplit)
    coef = data.frame(read.csv(paste("fold-", s, "-coeficientes.csv", sep="")))
    coef = coef[order(coef$coeficiente, decreasing = TRUE),]
    
    todos = rbind(todos, coef)
    
    um = coef[1,]
    primeiro = rbind(primeiro, um)
    
    dois = coef[2,]
    segundo = rbind(segundo, dois)
    
    tres = coef[3,]
    terceiro = rbind(terceiro, tres)
    
    quatro = coef[4,]
    quarto = rbind(quarto, quatro)
    
    cinco = coef[5,]
    quinto = rbind(quinto, cinco)
    
    seis = coef[6,]
    sexto = rbind(sexto, seis)
    
    s = s +1
    gc()
  }
  
  setwd(diretorios$folderReports)
  write.csv(todos, paste(dataset_name, "-All-Coefs.csv", sep=""), 
            row.names = FALSE)
  
  setwd(diretorios$folderReports)
  write.csv(primeiro, paste(dataset_name, "-1-Coefs.csv", sep=""), 
            row.names = FALSE)
  
  setwd(diretorios$folderReports)
  write.csv(segundo, paste(dataset_name, "-2-Coefs.csv", sep=""), 
            row.names = FALSE)
  
  setwd(diretorios$folderReports)
  write.csv(terceiro, paste(dataset_name, "-3-Coefs.csv", sep=""), 
            row.names = FALSE)
  
  setwd(diretorios$folderReports)
  write.csv(quarto, paste(dataset_name, "-4-Coefs.csv", sep=""), 
            row.names = FALSE)
  
  setwd(diretorios$folderReports)
  write.csv(quinto, paste(dataset_name, "-5-Coefs.csv", sep=""), 
            row.names = FALSE)
  
  setwd(diretorios$folderReports)
  write.csv(sexto, paste(dataset_name, "-6-Coefs.csv", sep=""), 
            row.names = FALSE)
  
  
  gc()
  cat("\n##################################################################################################")
  cat("\n# END OF THE FUNCTION HIGIEST COEFFICIENT                                                        #")
  cat("\n##################################################################################################")
  cat("\n\n\n\n")
}

#############################################################################
# Please, any errors, contact us: elainececiliagatto@gmail.com              #
# Thank you very much!                                                      #
#############################################################################
