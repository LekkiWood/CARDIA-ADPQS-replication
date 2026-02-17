#----------------------Build function -----------------------------------------#

#path_bridge = "/media/RawData/MESA/MESA-Phenotypes/MESA-Website-Phenos/MESA-SHARE_IDList_Labeled.csv"
#path_E1_covs = "/media/RawData/MESA/MESA-Phenotypes/MESA-Website-Phenos/MESAe1FinalLabel02092016.dta"
#path_E1_nutr =  "/media/RawData/MESA/MESA-Phenotypes/MESA-SHARe-Phenos/E1_nutrients_new.csv"
#path_cvd = "/media/RawData/MESA/MESA-Phenotypes/MESA-Website-Phenos/MESAEvThru2020AllCohort_20241120.dta"
#cleaned_proteins <- tar_read(Proteins_long)

#-----------------------------------------------------------------------------#
#-----------------------------------------------------------------------------#




build_traits_function <- function(path_E1_covs, path_E1_nutr, path_bridge, path_cvd, cleaned_proteins)
  
{
  
  #Bridging file
  bridge <- data.table::fread(path_bridge) |>
    dplyr::select(`SHARE ID Number`, `MESA Participant ID`) |>
    dplyr::rename(sidno = `SHARE ID Number`, idno = `MESA Participant ID`)
  

  
  #---------------------------Build outcomes file-----------------------#
  

  ###############
  # Diet
  ###############

  #Used for selection
  
  diet <- read.csv(path_E1_nutr) |>
    dplyr::select(sidno, enrgyn1c) |>
    dplyr::rename(energy = enrgyn1c)
  
  ###############
  # CVD
  ###############
  
  cvd <- foreign::read.dta(path_cvd) |>
    dplyr::mutate(cvdh = as.factor(dplyr::case_when(cvdh=="No" ~ 0,
                                                    cvdh=="Yes" ~ 1,
                                                   TRUE ~ NA_real_)
                                   )
                  ) |>
    dplyr::select(idno, cvdh, cvdhtt, dth, dthtype, dthtt)
  
  ###############
  # E1
  ###############
  
  E1 <- foreign::read.dta(path_E1_covs) |>
    dplyr::select(idno, age1c, gender1, race1c, bmi1c, dm031c, pkyrs1c, pamvcm1c, egfr1c, sbp1c, hdl1, chol1, htnmed1c, lipid1c) |>
    dplyr::rename(age = age1c,
                  BMI = bmi1c,
                  smoking = pkyrs1c,
                  PA = pamvcm1c,
                  eGFR = egfr1c, 
                  SBP = sbp1c,
                  HDL = hdl1,
                  cholesterol = chol1
                  ) |>
    dplyr::mutate(sex = as.factor(dplyr::case_when(gender1=="0: FEMALE" ~ 1,
                                                   gender1=="1: MALE" ~ 2,
                                                   TRUE ~ NA_real_)),
                  race = as.factor(dplyr::case_when(race1c == "1: white, CAUCASIAN" ~ 1,
                                                    race1c == "2: CHINESE-AMERICAN" ~ 2,
                                                    race1c == "3: black, AFRICAN-AMERICAN" ~ 3,
                                                    race1c == "4: HISPANIC" ~ 4,
                                                    TRUE ~ NA_real_)),
                  diabetes = dplyr::case_when(dm031c=="NORMAL" ~ 0,
                                              dm031c=="IFG" ~ 0,
                                              dm031c=="Untreated DIABETES" ~ 1,
                                              dm031c=="Treated DIABETES" ~ 1,
                                              TRUE ~ NA_real_),
                  htnmeds = as.factor(dplyr::case_when(htnmed1c=="No" ~0,
                                                      htnmed1c== "Yes" ~ 1,
                                                      TRUE ~ NA_real_)),
                  lipidmeds = as.factor(dplyr::case_when(lipid1c=="No" ~0,
                                                        lipid1c== "Yes" ~ 1,
                                                        TRUE ~ NA_real_))
                  )|>
                  dplyr::select(-race1c, -gender1, -dm031c, -htnmed1c, -lipid1c)

  
  #------------------------------------------------------------#
  #---------------Merge to full file --------------------------#
  #------------------------------------------------------------#
  
      Traits <- E1 |>
      dplyr::full_join(bridge, dplyr::join_by(idno)) |>
      dplyr::full_join(diet, dplyr::join_by(sidno)) |>
      dplyr::full_join(cvd, dplyr::join_by(idno)) 
  
  
    #------------------------------------------------------------#
    #---------------Select Pps.        --------------------------#
    #------------------------------------------------------------#    
    
    
    
#---------------------------N with protein data--------------------------------------------#
    
    protein_ids <- cleaned_proteins |>
      dplyr::filter(Exam==1 & !is.na(idno)) |>
      dplyr::select(idno, sidno)
    
    Traits <- Traits |>
      dplyr::filter(idno %in% protein_ids$idno)
    
#---------------------------N with protein & cvd data-------------------------------------#
    
    
    nocvd_ids <- Traits |>
      dplyr::filter(is.na(cvdhtt)) |>
      dplyr::select(idno, sidno)
    
    Traits <- Traits |>
      dplyr::filter(!idno %in% nocvd_ids$idno)
    
    cvd_ids <- Traits |>
      dplyr::select(idno, sidno)
    
#---------------------------N with protein & cvd & acceptable energy -----------------------#
    
    nodiet_ids <- Traits |>
      dplyr::filter(is.na(energy)) |>
      dplyr::select(idno, sidno)
    
    Traits <- Traits |>
      dplyr::filter(!idno %in% nodiet_ids$idno)
    
    highenergy_ids <- Traits |>
      dplyr::filter((sex==1 & (energy > 6000 | energy < 600)) | (sex==2 & (energy > 8000 | energy < 800))) |>
      dplyr::select(idno, sidno)
    
    Traits <- Traits |>
      dplyr::filter(!idno %in% highenergy_ids$idno)
    
    diet_ids <- Traits |>
      dplyr::select(idno, sidno)
    
  #------------------------------------------------------------#
  #---------------Info for Quarto -----------------------------#
  #------------------------------------------------------------#
  
  QC_info <- list(
    filenames = list(
      E1_covs_filename = path_E1_covs,
      E1_nutr_filename = path_E1_nutr,
      cvd_filename = path_cvd
                     ),
    IDs = list(
      protein_ids = protein_ids,
      nocvd_ids = nocvd_ids,
      cvd_ids = cvd_ids,
      nodiet_ids = nodiet_ids,
      highenergy_ids =  highenergy_ids,
      diet_ids = diet_ids
               )
    )
  
  #----------------Outputs -----------------------------#
  
  list(
    QC_info_out = QC_info,
    Traits_table = Traits
   )
  
  
}