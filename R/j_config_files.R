rm(list = ls())


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


###############################################################################
# LOAD LIBRARY/PACKAGE                                                        #
###############################################################################
library(stringr)


###############################################################################
# READING DATASET INFORMATION FROM DATASETS-ORIGINAL.CSV                      #
###############################################################################
setwd(FolderRoot)
datasets = data.frame(read.csv("datasets-original.csv"))
n = nrow(datasets)


###############################################################################
# CREATING FOLDER TO SAVE CONFIG FILES                                        #
###############################################################################
FolderCF = paste(FolderRoot, "/j-config-files", sep="")
if(dir.exists(FolderCF)==FALSE){dir.create(FolderCF)}


###############################################################################
# CREATING CONFIG FILES FOR EACH DATASET                                      #
###############################################################################
i = 1
while(i<=n){

  # specific dataset
  ds = datasets[i,]

  # print the dataset name
  cat("\ndataset = ", ds$Name)

  # Confi File Name
  file_name = paste(FolderCF, "/j-", ds$Name, ".csv", sep="")

  # Starts building the configuration file
  output.file <- file(file_name, "wb")

  # Config file table header
  write("Config, Value",
        file = output.file, append = TRUE)

  # Absolute path to the folder where the dataset's "tar.gz" is stored

  # write("Dataset_Path, \"/home/u704616/Datasets\"",
  #       file = output.file, append = TRUE)

  write("Dataset_Path, /home/elaine/Datasets",
        file = output.file, append = TRUE)

  # job name
  job_name = paste("gpj-", ds$Name, sep = "")

  # directory name

  # folder_name = paste("\"/scratch/", job_name, "\"", sep = "")
  # folder_name = paste("/scratch/", job_name, sep = "")
  # folder_name = paste("~/Exhaustive-MiF1-ECC/", job_name, sep = "")
  # folder_name = paste("~/tmp/", job_name, sep = "")
  folder_name = paste("/dev/shm/", job_name, sep = "")

  # Absolute path to the folder where temporary processing will be done.
  # You should use "scratch", "tmp" or "/dev/shm", it will depend on the
  # cluster model where your experiment will be run.
  str1 = paste("Temporary_Path, ", folder_name, sep="")
  write(str1,file = output.file, append = TRUE)

  # dataset name
  str3 = paste("dataset_name, ", ds$Name, sep="")
  write(str3, file = output.file, append = TRUE)

  # Dataset number according to "datasets-original.csv" file
  str2 = paste("number_dataset, ", ds$Id, sep="")
  write(str2, file = output.file, append = TRUE)

  # Number used for X-Fold Cross-Validation
  write("number_folds, 10", file = output.file, append = TRUE)

  # Number of cores to use for parallel processing
  write("number_cores, 10", file = output.file, append = TRUE)

  # finish writing to the configuration file
  close(output.file)

  # increment
  i = i + 1

  # clean
  gc()
}

###############################################################################
# Please, any errors, contact us: elainececiliagatto@gmail.com                #
# Thank you very much!                                                        #                                #
###############################################################################
