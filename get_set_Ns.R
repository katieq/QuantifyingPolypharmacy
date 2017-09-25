# Quinn & Shah, A dataset quantifying polypharmacy in the United States (2017)
# Function used to get N-combos where each of the N N-1 combos have >= min_exposure_count
# KQuinn 2016-2017

get_set_Ns <- function(N, exposure_count_Nminus1s, min_exposure_count = 100){
  
  if (N == 2){
    drug_ids = (exposure_count_Nminus1s %>% dplyr::filter(pw_count > min_exposure_count))$drug_id
    set_Ns = as.data.frame(t(combn(sort(drug_ids), 2))) %>%
      dplyr::select(drug_id_A = V1, drug_id_B = V2)
    
  } else if (N == 3) {
    exposure_count_Nminus1s = exposure_count_Nminus1s %>%
      dplyr::select(drug_id_A, drug_id_B, exposure_count) %>% # remove unhelpful columns for now.
      dplyr::filter(exposure_count > min_exposure_count)
    
    drug_1s = sort(unique(c(exposure_count_Nminus1s$drug_id_B, exposure_count_Nminus1s$drug_id_A)))
    
    set_Ns = exposure_count_Nminus1s %>%
      dplyr::filter(!(drug_id_B == max(drug_1s))) %>%
      dplyr::select(drug_id_A, drug_id_B) %>%
      dplyr::rowwise() %>%
      dplyr::do(data.frame(drug_id_A = .$drug_id_A,
                           drug_id_B = .$drug_id_B,
                           drug_id_C = drug_1s[which(drug_1s > .$drug_id_B)])) %>%
      dplyr::semi_join(exposure_count_Nminus1s, 
                       by = c("drug_id_A"="drug_id_A", "drug_id_C"="drug_id_B")) %>%
      dplyr::semi_join(exposure_count_Nminus1s,
                       by = c("drug_id_B"="drug_id_A", "drug_id_C"="drug_id_B"))
    
  } else if (N == 4){
    exposure_count_Nminus1s = exposure_count_Nminus1s %>%
      dplyr::select(drug_id_A, drug_id_B, drug_id_C, exposure_count) %>% # remove unhelpful columns for now.
      dplyr::filter(exposure_count > min_exposure_count) %>%
      dplyr::arrange(drug_id_A, drug_id_B, drug_id_C)
    
    drug_1s = sort(unique(c(exposure_count_Nminus1s$drug_id_A, 
                             exposure_count_Nminus1s$drug_id_B,
                             exposure_count_Nminus1s$drug_id_C)))

    set_Ns = exposure_count_Nminus1s %>%
      dplyr::filter(!(drug_id_C == max(drug_1s))) %>%
      dplyr::select(drug_id_A, drug_id_B, drug_id_C) %>%
      dplyr::rowwise() %>%
      dplyr::do(data.frame(drug_id_A = .$drug_id_A,
                           drug_id_B = .$drug_id_B,
                           drug_id_C = .$drug_id_C,
                           drug_id_D = drug_1s[which(drug_1s > .$drug_id_C)])) %>%
      dplyr::semi_join(exposure_count_Nminus1s, 
                       by = c("drug_id_A"="drug_id_A", "drug_id_B"="drug_id_B", "drug_id_D"="drug_id_C")) %>%
      dplyr::semi_join(exposure_count_Nminus1s,
                       by = c("drug_id_A"="drug_id_A", "drug_id_C"="drug_id_B", "drug_id_D"="drug_id_C")) %>%
      dplyr::semi_join(exposure_count_Nminus1s,
                       by = c("drug_id_B"="drug_id_A", "drug_id_C"="drug_id_B", "drug_id_D"="drug_id_C"))
    
  } else if (N == 5){
    exposure_count_Nminus1s = exposure_count_Nminus1s %>%
      dplyr::filter(exposure_count > min_exposure_count) %>%
      dplyr::select(drug_id_A, drug_id_B, drug_id_C, drug_id_D) # remove unhelpful columns for now.
    
    drug_1s = sort(unique(c(exposure_count_Nminus1s$drug_id_A, 
                             exposure_count_Nminus1s$drug_id_B,
                             exposure_count_Nminus1s$drug_id_C,
                             exposure_count_Nminus1s$drug_id_D)))
    
    set_Ns = exposure_count_Nminus1s %>%
      dplyr::filter(!(drug_id_D == max(drug_1s))) %>%
      dplyr::select(drug_id_A, drug_id_B, drug_id_C, drug_id_D) %>%
      dplyr::rowwise() %>%
      dplyr::do(data.frame(drug_id_A = .$drug_id_A,
                           drug_id_B = .$drug_id_B,
                           drug_id_C = .$drug_id_C,
                           drug_id_D = .$drug_id_D,
                           drug_id_E = drug_1s[which(drug_1s > .$drug_id_D)])) %>%
      dplyr::semi_join(exposure_count_Nminus1s, 
                       by = c("drug_id_A"="drug_id_A", "drug_id_B"="drug_id_B", "drug_id_C"="drug_id_C", "drug_id_E"="drug_id_D")) %>%
      dplyr::semi_join(exposure_count_Nminus1s,
                       by = c("drug_id_A"="drug_id_A", "drug_id_B"="drug_id_B", "drug_id_D"="drug_id_C", "drug_id_E"="drug_id_D")) %>%
      dplyr::semi_join(exposure_count_Nminus1s,
                       by = c("drug_id_A"="drug_id_A", "drug_id_C"="drug_id_B", "drug_id_D"="drug_id_C", "drug_id_E"="drug_id_D")) %>%
      dplyr::semi_join(exposure_count_Nminus1s,
                       by = c("drug_id_B"="drug_id_A", "drug_id_C"="drug_id_B", "drug_id_D"="drug_id_C", "drug_id_E"="drug_id_D"))
    
  }
}


