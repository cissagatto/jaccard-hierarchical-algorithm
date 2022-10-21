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

#cores = seecol(pal_unikn_pref, 200)


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
    
    ########################################################################
    # TODOS OS MÉTODOS
    fold = c(0)
    average = c(0)
    complete = c(0)
    mcquitty = c(0)
    single = c(0)
    wardD = c(0)
    wardD2 = c(0)
    c4 = data.frame(fold, average, complete, mcquitty, single, wardD, wardD2)
    
    
    ########################################################################
    # APENAS OS COEFICIENTES
    fold = c(0)
    metodo = c(0)
    coeficiente = c(0)
    coefHC = data.frame(fold, metodo, coeficiente)
    
    
    #####################################################################
    # TODAS AS PARTIÇÕES COM RÓTULOS
    num.fold = c(0)
    num.part = c(0)
    num.group = c(0)
    metodo.name = c("")
    names.labels = c(0)
    AllPartitions = data.frame(num.fold, num.part, num.group, 
                               names.labels, metodo.name)
    
    
    #####################################################################
    # FREQUENCIA
    num.fold = c(0)
    metodo.name = c("")
    num.part = c(0)
    num.group = c(0)
    frequencia = c(0)
    final = data.frame(num.fold, num.part, num.group, frequencia, metodo.name)
    
    #####################################################################
    # RESUMO DA PARTIÇÃO
    num.fold = c(0)
    metodo.name = c("")
    num.part = c(0)
    num.group = c(0)
    resumePartitions = data.frame(num.fold, num.part, num.group, metodo.name)
    
    #####################################################################
    # methods
    metodos = c("average", "complete", "mcquitty", "single", "ward.D", "ward.D2")    
    
    
    #####################################################################
    # created folder for the split
    FolderHClust = paste(diretorios$folderPartitions, "/Split-", s, sep="")
    if(dir.exists(FolderHClust)==FALSE){dir.create(FolderHClust)}
    
    #FolderOS = paste(diretorios$folderOutputDataset, "/Split-", s, sep="")
    #if(dir.exists(FolderOS)==FALSE){dir.create(FolderOS)}
    
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
      
      FolderMOS = paste(diretorios$folderOutputDataset, "/", metodos[i], sep="")
      if(dir.exists(FolderMOS)==FALSE){dir.create(FolderMOS)}
      
      FolderSMOS = paste(FolderMOS, "/Split-", s, sep="")
      if(dir.exists(FolderSMOS)==FALSE){dir.create(FolderSMOS)}
      
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
      
      #####################################################################
      group = c(0)
      label = c(0)
      allPartitions2 = data.frame(label, group)
      
      fold = c(0)
      partition = c(0)
      num.groups = c(0)
      grupos_por_particao = data.frame(fold, partition, num.groups)
      
      k = 1
      for(k in 1:ds$Labels){
        cat("\ncluster: ", k)   
        
        FolderPOS = paste(FolderSMOS , "/Partition-", k, sep="")
        if(dir.exists(FolderPOS)==FALSE){dir.create(FolderPOS)}
        
        group = c(0)
        label = c("")
        clusters3 = cbind(group, label)
             
        cutLabels = cutree(HC, k)
        clusters = data.frame(cutree(HC, k))
        names(clusters) = "grupo"
        label = c(rownames(clusters))
        
        group = c(clusters$grupo)
        label = label
        clusters3 = data.frame(group, label)
        
        #cat("\nSAVE CUTREE\n")
        setwd(FolderClusters)
        write.csv(clusters3, paste("partition-", k, ".csv", sep=""), 
                  row.names = FALSE)		
        
        setwd(FolderPOS)
        write.csv(clusters3, paste("partition-", k, ".csv", sep=""), 
                  row.names = FALSE)		
        
        #cat("\nFrequencia")
        frequencia1 = count(clusters3, clusters3$group)
        names(frequencia1) = c("grupo", "frequencia")
        
        freq = frequencia1
        names(freq) = c("group", "totalLabels")
        setwd(FolderPOS)
        write.csv(freq, 
                  paste("fold-", s, "-labels-per-group-partition-",
                        k, ".csv", sep=""), 
                  row.names = FALSE)	
        
        fold = s
        partition = k
        num.groups = k
        teste = data.frame(fold, partition, num.groups)
        grupos_por_particao = rbind(grupos_por_particao, teste)
        
        
        num.fold = s
        metodo.name = metodos[i]
        num.part = k
        num.group = frequencia1$grupo
        frequencia = frequencia1$frequencia
        final = rbind(final, data.frame(num.fold, num.part, 
                                        num.group, frequencia, metodo.name))
        
        ############################################################################################################
        # cat("\nData frame")
        num.fold = s
        metodo.name = metodos[i]
        num.part = k
        num.group = k
        resumePartitions = rbind(resumePartitions, data.frame(num.fold, 
                                                              num.part, 
                                                              num.group,
                                                              metodo.name))
        
        ############################################################################################################       
        num.fold = s
        num.part = k
        num.group = clusters3$group
        metodo.name = metodos[i]
        names.labels = clusters3$label
        AllPartitions = rbind(AllPartitions, data.frame(num.fold, num.part, 
                                                         num.group,
                                                         names.labels,
                                                        metodo.name))
        
        ############################################################################################################
        nomesDosRotulos = clusters3$rotulos
        group = clusters3$group
        allPartitions2 = cbind(allPartitions2, group)
        b = k + 2
        names(allPartitions2)[b] = paste("partition-", k, sep="")
        
        k = k + 1 
        gc()
        
      } # fim do cluster
      
      allPartitions2 = allPartitions2[,c(-1,-2)]
      
      write.csv(allPartitions2, 
                paste(FolderMethods, "/fold-", s, "-", metodos[i], 
                      "-all-partitions.csv", sep=""), 
                row.names = FALSE)  
      
      write.csv(allPartitions2, 
                paste(FolderSMOS, "/fold-", s, 
                      "-all-partitions.csv", sep=""), 
                row.names = FALSE)  
      
      grupos_por_particao = grupos_por_particao[c(-1,-2),]
      w = nrow(grupos_por_particao)
      grupos_por_particao = grupos_por_particao[-w,]
      
      write.csv(grupos_por_particao, 
                paste(FolderSMOS, "/fold-", s, 
                      "-groups-per-partition.csv", sep=""), 
                row.names = FALSE) 
      
      write.csv(grupos_por_particao, 
                paste(FolderMethods, "/fold-", s, 
                      "-groups-per-partition.csv", sep=""), 
                row.names = FALSE) 
      
      print(system(paste("rm -r ", FolderSMOS, "/Partition-1", sep="")))
      print(system(paste("rm -r ", FolderSMOS, "/Partition-", ds$Labels, sep="")))
      
      
      i = i + 1 # increment hclust method
      gc() # garbage collection
      
    } # fim do método
    
    #cat("\nSAVE ALL COEFICIENTES\n")        
    setwd(FolderHClust)
    write.csv(coefHC[-1,], paste("fold-", s, "-coeficientes.csv", sep=""), 
              row.names = FALSE)
    
    write.csv(resumePartitions[-1,], 
              paste(FolderHClust, "/fold-", s, "-groups-per-partition.csv", sep=""), 
              row.names = FALSE)
    
    write.csv(final[-1,], 
              paste(FolderHClust, "/fold-", s, "-frequency.csv", sep=""), 
              row.names = FALSE)
    
    
    gc()
  } # fim do fold
  
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
  
  write.csv(todos, paste(diretorios$folderPartitions, "/", 
                         dataset_name, "-All-Coefs.csv", sep=""), 
            row.names = FALSE)
  
  write.csv(primeiro, paste(diretorios$folderPartitions, "/", 
                             dataset_name, "-1-Coefs.csv", sep=""), 
            row.names = FALSE)
  
  write.csv(segundo,paste(diretorios$folderPartitions, "/", 
                            dataset_name, "-2-Coefs.csv", sep=""), 
            row.names = FALSE)
  
  write.csv(terceiro, paste(diretorios$folderPartitions, "/", 
                             dataset_name, "-3-Coefs.csv", sep=""), 
            row.names = FALSE)
  
  write.csv(quarto, paste(diretorios$folderPartitions, "/", 
                           dataset_name, "-4-Coefs.csv", sep=""), 
            row.names = FALSE)
  
  write.csv(quinto, paste(diretorios$folderPartitions, "/", 
                          dataset_name, "-5-Coefs.csv", sep=""), 
            row.names = FALSE)
  
  write.csv(sexto, paste(diretorios$folderPartitions, "/", 
                          dataset_name, "-6-Coefs.csv", sep=""), 
            row.names = FALSE)
  
  a = data.frame(count(primeiro, primeiro$metodo))
  colnames(a) = c("metodo", "total")
  a = cbind(ordem = "1º.Lugar", a)
  
  b = data.frame(count(segundo, segundo$metodo))
  colnames(b) = c("metodo", "total")
  b = cbind(ordem = "2º.Lugar", b)
  
  c = data.frame(count(terceiro, terceiro$metodo))
  colnames(c) = c("metodo", "total")
  c = cbind(ordem = "3º.Lugar", c)
  
  d = data.frame(count(quarto, quarto$metodo))
  colnames(d) = c("metodo", "total")
  d = cbind(ordem = "4º.Lugar", d)
  
  e = data.frame(count(quinto, quinto$metodo))
  colnames(e) = c("metodo", "total")
  e = cbind(ordem = "5º.Lugar", e)
  
  f = data.frame(count(sexto, sexto$metodo))
  colnames(f) = c("metodo", "total")
  f = cbind(ordem = "6º.Lugar", f)
  
  all = rbind(a,b,c,d,e,f)
  
  write.csv(all, paste(diretorios$folderPartitions, "/", 
                          dataset_name, "-best-frequency-methods.csv", sep=""), 
            row.names = FALSE)
  
  final = cbind(primeiro, segundo, terceiro, quarto, quinto, sexto)
  final = final[,c(-4,-7,-10,-13, -16)]
  names(final) = c("fold", 
                      "1º.Lugar", "1.coef", 
                      "2º.Lugar", "2.coef",
                      "3º.Lugar", "3.coef",
                      "4º.Lugar", "4.coef",
                      "5º.Lugar", "5.coef",
                      "6º.Lugar", "6.coef")
  
  write.csv(final, paste(diretorios$folderPartitions, "/", 
                       dataset_name, "-todos.csv", sep=""), 
            row.names = FALSE)
  
  
  gc()
  cat("\n##################################################################################################")
  cat("\n# END OF THE FUNCTION HIGIEST COEFFICIENT                                                        #")
  cat("\n##################################################################################################")
  cat("\n\n\n\n")
}



analisaParticoes <- function(ds,
                            dataset_name,
                            number_dataset, 
                            number_cores, 
                            number_folds, 
                            folderResults,
                            resLS,
                            namesLabels){
  
  diretorios = directories(dataset_name, folderResults)
  metodos = c("average", "complete", "mcquitty", "single", "ward.D", "ward.D2")    
  FolderPart = paste(diretorios$folderReports, "/Partitions", sep="")
  
  todos = data.frame()
  final_labels = data.frame()
  final_groups = data.frame()
  final_metodos = data.frame()
  
  rotulos2 = c()
  
  f = 1
  while(f<=number_folds){
    cat("\n\nFOLD:", f)
    FolderSplit = paste(FolderPart, "/Split-", f, sep="")
    num.part = ds$Labels - 1
    
    fold_labels = data.frame()
    fold_groups = data.frame()
    fold_metodos = data.frame()
    
    p = 2
    while (p <= num.part) {
      cat("\n\tPARTITION:", p)
      
      FolderP = paste(FolderSplit, "/Partition-", p, sep="")
      if(dir.exists(FolderP)==FALSE){dir.create(FolderP)}
      
      apagar = c(0)
      part_metodos = data.frame(apagar)
      part_labels = data.frame(apagar)
      part_groups = data.frame(apagar)
      
      rotulos = c()
      grupos =c()
      
      m = 1
      while (m <= length(metodos)) {
        cat("\n\t\tMETODO:", metodos[m])
        FolderMetodo = paste(FolderSplit, "/", metodos[m], sep = "")
        FolderCluster = paste(FolderMetodo, "/Clusters", sep = "")
        
        setwd(FolderCluster)
        res = data.frame(read.csv(paste("partition-", p, ".csv", sep = "")))
        
        groups = res$group
        labels = res$label
        
        rotulos2 = labels
        rotulos = labels
        grupos = groups
        
        part_groups = cbind(part_groups, groups)
        part_labels = cbind(part_labels, labels)
        
        #names(part_groups)[2] = toString(metodos[m])
        #names(part_labels)[2] = toString(metodos[m])
        
        res2 = cbind(fold = f, partition = p,
                     method = metodos[m], res)
        todos = rbind(todos, res2)
        
        colnames(res)[1] = paste(metodos[m], ".groups", sep="")
        colnames(res)[2] = paste(metodos[m], ".labels", sep="")
        
        part_metodos = cbind(part_metodos, res)
        
        rm(res)
        m = m + 1
        gc()
      } # fim do método
      
      part_metodos = part_metodos[,-1]
      part_labels = part_labels [,-1]
      part_groups = part_groups[,-1]
      
      colnames(part_labels) = metodos
      colnames(part_groups) = metodos
      
      part_labels = cbind(fold = f, partition = p, 
                          groups = grupos, part_labels)
      
      part_groups = cbind(fold = f, partition = p, 
                          labels = rotulos, part_groups)
      
      fold_groups = rbind(fold_groups, part_groups)
      fold_labels = rbind(fold_labels, part_labels)
      fold_metodos = rbind(fold_metodos, part_metodos)
      
      write.csv(fold_metodos, paste(FolderP, "/fold-", f, 
                                    "-part-", p, "-groups-labels.csv", sep=""), 
                row.names = FALSE)
      
      write.csv(fold_groups, paste(FolderP, "/fold-", f, 
                          "-part-", p, "-groups.csv", sep=""), 
                row.names = FALSE)
      
      write.csv(fold_labels, paste(FolderP, "/fold-", f, 
                          "-part-", p, "-labels.csv", sep=""), 
                row.names = FALSE)
      
      p = p + 1
      gc()
    } # fim das partições
    
    final_labels = rbind(final_labels, fold_labels)
    final_groups = rbind(final_groups, fold_groups)
    final_metodos = rbind(final_metodos , fold_metodos)
    
    f = f + 1
    gc()
  } # fim do fold
  
  write.csv(final_labels, paste(diretorios$folderPartitions, 
                                "/all-labels.csv", sep=""), 
            row.names = FALSE)
  
  write.csv(final_groups, paste(diretorios$folderPartitions, 
                                "/all-groups.csv", sep=""), 
            row.names = FALSE)
  
  write.csv(final_metodos, paste(diretorios$folderPartitions, 
                                "/all-methods.csv", sep=""), 
            row.names = FALSE)

  
  num.part = ds$Labels-1
  cat("\n")
  a = 2
  while(a<=num.part){
    cat("\nPARTITION: ",a)
    res1 = filter(final_labels, final_labels$partition == a)
    res2 = filter(final_groups, final_groups$partition == a)
    setwd(diretorios$folderPartitions)
    write.csv(res1, paste("all-partitions-", a, "-labels.csv", sep=""), 
              row.names = FALSE)
    write.csv(res2, paste("all-partitions-", a, "-groups.csv", sep=""), 
              row.names = FALSE)
    a = a + 1
    gc()
  }
  
  
  b = 1
  while(b<=ds$Labels){
    #cat("\nB", b)
    
    rotulo = rotulos2[b]
    
    cat("\nLABEL:", rotulo)
    
    setwd(diretorios$folderPartitions)
    grupos = data.frame(read.csv("all-groups.csv"))
    
    gr = filter(grupos, grupos$labels == rotulo)
    gr = gr[order(gr$partition, decreasing = TRUE),]
    
    
    setwd(diretorios$folderPartitions)
    write.csv(gr, paste(rotulo, ".csv", sep=""), 
              row.names = FALSE)
   
    todos = data.frame()
    num.part = ds$Labels-1
    
    g = 2
    while(g<=num.part){
      
      #cat("\nGROUP:", g)
      
      gr2 = filter(gr, gr$partition == g)
      
      av = data.frame(count(gr2, gr2$average))
      colnames(av) = c("grupo", "total")
      av = cbind(metodo = "average", av)
      
      cm = count(gr2, gr2$complete)
      colnames(cm) = c("grupo", "total")
      cm = cbind(metodo = "complete", cm)
      
      mc = count(gr2, gr2$mcquitty)
      colnames(mc) = c("grupo", "total")
      mc = cbind(metodo = "mcquitty", mc)
      
      sin = count(gr2, gr2$single)
      colnames(sin) = c("grupo", "total")
      sin = cbind(metodo = "single", sin)
      
      wd = count(gr2, gr2$ward.D)
      colnames(wd) = c("grupo", "total")
      wd = cbind(metodo = "wardD", wd)
      
      wd2 = count(gr2, gr2$ward.D2)
      colnames(wd2) = c("grupo", "total")
      wd2 = cbind(metodo = "wardD2", wd2)
      
      all = rbind(av, cm, mc, sin, wd, wd2)
      all = cbind(partition = g, all)
      
      todos = rbind(todos, all)
      
      g = g + 1
      gc()
    }
    
    setwd(diretorios$folderPartitions)
    write.csv(todos, paste(rotulo, "-frequency.csv", sep=""), 
              row.names = FALSE)
    
    b = b + 1
    gc()
  }
  
  
  cor = "#ed2131"
  Folder = paste(diretorios$folderPartitions, "/Plots", sep="")
  if(dir.exists(Folder)==FALSE){dir.create(Folder)}
  setwd(Folder)
  pdf("label-distribution.pdf", width = 10, height = 7)
  op <- par(mar = c(2.5, 2.5, 2.5, 2.5), oma=c(1.5,1.5,1.5,1.5))
  par(mfrow = c(ds$Labels,ds$Labels))
    b = 1
    while(b<=ds$Labels){
      rotulo = rotulos2[b]
      setwd(diretorios$folderPartitions)
      grupos = data.frame(read.csv("all-groups.csv"))
      gr = filter(grupos, grupos$labels == rotulo)
      gr = gr[order(gr$partition, decreasing = TRUE),]
      num.part = ds$Labels-1
      g = 2
      while(g<=num.part){
        gr2 = filter(gr, gr$partition == g)
        plot(gr2$average, type="b" , lwd=2 , col=cor, 
             ylab="Groups" , xlab="Folds" , 
             main=paste(rotulo, " partition-", g, sep=""),
             bty="l" , pch=20 , cex=1, cex.main=0.7,
             cex.axis = 0.6, mgp = c(1.5, 0.5, 0), cex.lab = 0.7,
             xlim = c(1, 10), ylim = c(1, g), )
        #abline(h=seq(0,100,10) , col=cor, lwd=0.1)  
        g = g + 1
        gc()
      }
      b = b + 1
      gc()
    }
  par(op)
  print(par(op))
  dev.off()
  cat("\n")
  
  
  cat("\n\n#########################################################")
  cat("\n# FIM ANALISA PARTICOES                                   #")
  cat("\n###########################################################\n\n")
  
  gc()
}


separaBests <- function(ds,
                            dataset_name,
                            number_dataset, 
                            number_cores, 
                            number_folds, 
                            folderResults,
                            resLS,
                            namesLabels){
  
  diretorios = directories(dataset_name, folderResults)
  
  setwd(diretorios$folderPartitions)
  best = data.frame(read.csv(paste(dataset_name,"-1-Coefs.csv", sep="")))
  
  FolderDestino = paste(diretorios$folderOutput, 
                        "/Best", sep="")
  if(dir.exists(FolderDestino)==FALSE){dir.create(FolderDestino)}
  
  f = 1
  while(f<=number_folds){
    best_fold = best[f,]
    cat("\n", f, " - ", best_fold$metodo)
    FolderOrigem = paste(diretorios$folderOutputDataset, 
                         "/", best_fold$metodo, 
                         "/Split-", f, "/*", sep="")
    FolderD = paste(FolderDestino, "/Split-", f, sep="")
    if(dir.exists(FolderD)==FALSE){dir.create(FolderD)}
    
    print(system(paste("cp -r ", FolderOrigem, " ", FolderD, sep="")))
    
    f = f + 1
    gc()
  }
  
  
}




#############################################################################
# Please, any errors, contact us: elainececiliagatto@gmail.com              #
# Thank you very much!                                                      #
#############################################################################
