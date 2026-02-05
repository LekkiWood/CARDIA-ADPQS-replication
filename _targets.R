library('targets')
library('tarchetypes')
library(crew)
tar_option_set(controller = crew_controller_local(workers = 10)) 
tar_option_set(seed = 11042012)
Sys.setenv(VROOM_CONNECTION_SIZE = as.character(10 * 1024 * 1024)) #For any large datasets

####This runs a workflow for cleaning the MESA TOPMed multi-omics project proteomics data (O-link), 
#### formatting it into the following format: wide for metabolites & long for exams, 
#### and providing basic quality control (QC) metrics



tar_option_set(packages = c("dplyr", "tidyr", "tibble", "readr", "data.table", "bit64", 
                            "foreign", "quarto", "rlang", "purrr", "rcompanion", "knitr", "gtsummary", "labelled",
                            "kableExtra", "gt", "cli", "quarto", "sandwich", "lme4", "lmerTest", "ggplot2", "glmnet", "fastDummies", "doParallel"))

tar_source("/media/Analyses/CARDIA-ADPQS-replication/R")

list(
  #---------------------------------------------------------------------------------------#
  #--------------------------------1. Build proteomics table------------------------------#
  #---------------------------------------------------------------------------------------#
  
  #--------------------------------1A. Initial table
  
  #Protein intensities
  tar_target(path_proteins, "/media/RawData/MESA/MESA-Multiomics/MESA-Multiomics_Proteomics/Olink/SMP_IntensityNormalized_20251005.csv", format = "file"),
  #Mapping info
  tar_target(path_protein_info, "/media/RawData/MESA/MESA-Multiomics/MESA-Multiomics_Proteomics/Olink/Mapping_SMP_Plate_20251005.csv", format = "file"),
  #Protein keys
  tar_target(path_protein_keys, "/media/RawData/MESA/MESA-Multiomics/MESA-Multiomics_Proteomics/Olink/MESAOlink3k_proteinKeys_03292023.csv", format = "file"),
  
  #Bridging file
  tar_target(path_bridge,"/media/RawData/MESA/MESA-Phenotypes/MESA-Website-Phenos/MESA-SHARE_IDList_Labeled.csv", format = "file"),
  
  #Run function
  tar_target(build_proteins_out,   
             build_protein_table_function(
               path_proteins = path_proteins,
               path_protein_info = path_protein_info,
               path_protein_keys = path_protein_keys,
               path_bridge = path_bridge)
  ),
  
  #Save outputs
  #QC info
  tar_target(build_proteins_QC, build_proteins_out$QC_info_out),
  #Proteins
  tar_target(Proteins_long, build_proteins_out$Formatted_proteins_out),
  #N in formatted file by exam
  tar_target(Proteins_long_N, build_proteins_out$N_by_exam),
 
   #--------------------------------1B. Mapping file
  
  #QC info
  tar_target(protein_info, QC_proteins_function(path_proteins = path_proteins,
                                                path_protein_info = path_protein_info,
                                                path_protein_keys = path_protein_keys,
                                                Proteins_long = Proteins_long)),
  tar_target(Protein_mapping_file, protein_info$final_proteins_mapping),
  
  #--------------------------------1C. Final clean protein table
  
  #Build proteins
  tar_target(proteins_clean, final_proteins_function(Proteins_long = Proteins_long, 
                                                     QC_file = Protein_mapping_file)
  ),
  
  
  #Save proteins
  tar_target(Proteins_long_clean, proteins_clean$Final_proteins),
  #N in cleaned & formatted file by exam
  tar_target(Proteins_clean_N, proteins_clean$N_by_exam),
  
  
  #---------------------------------------------------------------------------------------#
  #--------------------------------2. Build traits table----------------------------------#
  #---------------------------------------------------------------------------------------#
  
  #------------Files-------
  #E1_covs
  tar_target(path_E1_covs, "/media/RawData/MESA/MESA-Phenotypes/MESA-Website-Phenos/MESAe1FinalLabel02092016.dta", format = "file"),
  
  #E1_nutrients
  tar_target(path_E1_nutr, "/media/RawData/MESA/MESA-Phenotypes/MESA-SHARe-Phenos/E1_nutrients_new.csv", format = "file"),
  
  #CVD events
  tar_target(path_cvd, "/media/RawData/MESA/MESA-Phenotypes/MESA-Website-Phenos/MESAEvThru2020AllCohort_20241120.dta", format = "file"),
  
  #------------Function-------
  tar_target(build_traits, build_traits_function(path_E1_covs = path_E1_covs,
                                                 path_E1_nutr = path_E1_nutr,
                                                 path_bridge = path_bridge,
                                                 path_cvd = path_cvd,
                                                 cleaned_proteins = Proteins_long)),

  #------------Outputs-------
  
  #Filenames and missingness
  tar_target(traits_QC_info, build_traits$QC_info_out),
  #Trait data
  tar_target(traits_db, build_traits$Traits_table),
  
  #---------------------------------------------------------------------------------------#
  #--------------------------------Quarto output------------------------------------------#
  #---------------------------------------------------------------------------------------#
  
  tarchetypes::tar_quarto(
  build_quarto,
  path = "/media/Analyses/CARDIA-ADPQS-replication/CARDIA-ADPQS-replication.qmd",
  quiet = FALSE
  )
  
  
  
)
  

  
 