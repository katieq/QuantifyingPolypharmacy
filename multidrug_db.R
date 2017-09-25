# Quinn & Shah, A dataset quantifying polypharmacy in the United States (2017)
# Code to generate Data Record 1: Drug ingredient combinations
# KQuinn 2016-2017

library(dplyr)
library(tidyr)
library(data.table)

# To protect the privacy of the format of Truven Health Marketscan Databases, this available 
# code begins with prescription data in the format:
# data_drug : patient_id | drug_id | age_in_days | days_supply
# Where patient_id is a patient identifier
#       drug_id is a drug identifier (for a local mapping)
#       age_in_days is the age of the prescription, in days
#       days_supply is the days of supply in the prescription.

# Load data and define time window variable -----------------------------------------------------
# Load the data frame data_drug, which was extracted from the Truven Health Marketscan Database,
# with columns: patient_id | age_in_days | days_supply | drug_id

# To-do here: load your dataset

# Load the drug mapping (Data Record 3), saved as the data frame Q_drug_id
# with columns: drug_id | drug_name | rxcui | atc | atc_class_name | estimate_drug_cost_per_day

# To-do here: load drug mapping as Q_drug_id and atc_classes

# Set variable:
t_window = 30 # Select 30-day non-overlapping windows to define concommitant drug exposure.


# Scan prescriptions into exact-exposure windows ----------------------------------------------------------
# Convert prescription timing from age_in_days to windows:
data_drug = data_drug %>%
  dplyr::mutate(window_id = (age_in_days %/% t_window) * t_window) %>%
  dplyr::select(patient_id, window_id, drug_id, days_supply)

# Create additional exposures for cases where the prescription days_supply extends beyond one window:

scan_data <- function(data_drug, t_window = fix_twindow){
  # As an approximation: very few days_supply are longer than 90 days, so limit to adding up to 3 windows

  data_drug$days_supply[which(data_drug$days_supply > 120)] = 120
  data_drug$days_supply[which(data_drug$days_supply < 1)] = 1
  
  n_windows_create = 3

  data_drug_scan = data_drug
  for (i in 1:n_windows_create){
    data_drug_scan = data_drug_scan %>%
      dplyr::bind_rows(data_drug %>%
                         dplyr::filter(days_supply > i*t_window) %>%
                         dplyr::mutate(window_id = window_id + i*t_window))
  }
  drug_exposures = data_drug_scan %>%
    dplyr::select(patient_id, window_id, drug_id) %>%
    dplyr::distinct()
}

drug_exposures = scan_data(data_drug, t_window = t_window)
remove(data_drug)

# Group exposures by window, representing as a string of drug_id, separated by "_"
get_exact_drug_code <- function(drug_id_vec){
  drug_id_vec <- sort(drug_id_vec)
  exact_drug_code <- paste0("_",paste(drug_id_vec, collapse = "_"),"_")
}

exact_drugs <- drug_exposures %>%
  dplyr::group_by(patient_id, window_id) %>%
  dplyr::summarise(exact_drug_code = get_exact_drug_code(drug_id),
                   n_pw_drugs = n())

summary_exact_drugs <- exact_drugs %>% 
  dplyr::count(exact_drug_code) %>%
  dplyr::rename(exact_exposure_count = n) %>%
  dplyr::arrange(-exact_exposure_count) %>%
  dplyr::mutate(n_drugs = stringr::str_count(exact_drug_code, "_") - 1) 


# Note: If data_drugs is very large, it can be processed in batches, then combined after this step.





# Count at-least-exposure for N-drug combinations ----------------------------------------------

# shift this to separate file run from here? probably not.
# or shift subsections to a different file? e.g. get_atleast_1s, get_atleast_2s, ... 

# Create index, storing the rows containing each drug_id, to enable fast counting of atleast_exposure
index_row_drug_ndrug_list = list()

tmp <- summary_exact_drugs %>%
  dplyr::mutate(rn = row_number())

for (i in 1:max(summary_exact_drugs$n_drugs)){
  
  print(paste0("Index for n_drugs = ",i," ..."))
  
  summary_exact_drugs_i = tmp %>% dplyr::filter(n_drugs == i)
  
  index_row_drug_ndrug_list[[i]] = summary_exact_drugs_i %>%
    tidyr::separate(exact_drug_code, as.character(c(0:(i+1))), sep="_", extra="warn",fill="warn",remove=TRUE, convert=TRUE) %>%
    tidyr::gather(row_pos, drug_id, 1+c(1:i), na.rm=TRUE) %>%
    dplyr::select(rn, drug_id, exact_exposure_count)
}
index_row_drug = as.data.frame(rbindlist(index_row_drug_ndrug_list))
remove(index_row_drug_ndrug_list, i, tmp)

# Count single drug atleast exposures from the index:
atleast_1s = index_row_drug %>%
  dplyr::group_by(drug_id) %>%
  dplyr::summarise(atleast_exposure_count = sum(exact_exposure_count)) %>%
  dplyr::arrange(-atleast_exposure_count)

# To count atleast exposure to N-drug combinations, 
# Create an index of location of drug_id 
index_row_drug_list <- list()
for (i in 1:max(Q_drug_id$drug_id)){
  # Progress:
  if (i %% 100 == 0){
    print(paste(round(i/max(Q_drug_id$drug_id),2)*100,"%"))
  }
  index_row_drug_list[[i]] <- (index_row_drug %>% dplyr::semi_join(data.frame(drug_id = i), by = "drug_id"))$rn
}

# Load function that makes a set of N-drug combinations,
# where each of the N N-1 drug subsets have at least min_exposure_count exposures.
source("get_set_Ns.R") 

get_2s_count <- function(drug_id_A_in, drug_id_B_in){
  sum(summary_exact_drugs[intersect(index_row_drug_list[[drug_id_A_in]],
                                    index_row_drug_list[[drug_id_B_in]]),"exact_exposure_count"])
}

print(" COUNTING 2s ... ")
set_2s <- get_set_Ns(2, atleast_1s, min_exposure_count = 1e3) # all 2-drug combos made from single drugs with >1e3 exposures
print(paste0("Generated set_2s, with ", dim(set_2s)[1], " rows"))

atleast_2s <- set_2s %>% # columns need to be drug_id_A and drug_id_B
  dplyr::rowwise() %>%
  dplyr::mutate(atleast_exposure_count = get_2s_count(drug_id_A, drug_id_B))
save(atleast_2s, file = "truven_atleast_2s.RData")

# _3s
get_3s_count <- function(drug_id_A_in, drug_id_B_in, drug_id_C_in){
  sum(summary_exact_drugs[Reduce(intersect, list(index_row_drug_list[[drug_id_A_in]],
                                                 index_row_drug_list[[drug_id_B_in]],
                                                 index_row_drug_list[[drug_id_C_in]])),"exact_exposure_count"])
}

set_3s <- get_set_Ns(3, atleast_2s, min_exposure_count = 1e4) 
print(paste0("Generated set_3s, with ", dim(set_3s)[1], " rows, from 2s " ))

print(" COUNTING 3s ... ")
atleast_3s <- set_3s %>%
  dplyr::rowwise() %>%
  dplyr::mutate(atleast_exposure_count = get_3s_count(drug_id_A, drug_id_B, drug_id_C))
save(atleast_3s, file = "truven_atleast_3s.RData")

# _4s
get_4s_count <- function(drug_id_A_in, drug_id_B_in, drug_id_C_in, drug_id_D_in){
  sum(summary_exact_drugs[Reduce(intersect, list(index_row_drug_list[[drug_id_A_in]],
                                                 index_row_drug_list[[drug_id_B_in]],
                                                 index_row_drug_list[[drug_id_C_in]],
                                                 index_row_drug_list[[drug_id_D_in]])),"exact_exposure_count"])
}

set_4s <- get_set_Ns(4, atleast_3s, min_exposure_count = 1e4)
print(paste0("Generated set_4s, with ", dim(set_4s)[1], " rows, from 3s " ))

print(" COUNTING 4s ... ")
atleast_4s <- set_4s %>%
  dplyr::rowwise() %>%
  dplyr::mutate(atleast_exposure_count = get_4s_count(drug_id_A, drug_id_B, drug_id_C, drug_id_D))
save(atleast_4s, file = "truven_atleast_4s.RData")

# _5s
get_5s_count <- function(drug_id_A_in, drug_id_B_in, drug_id_C_in, drug_id_D_in, drug_id_E_in){
  sum(summary_exact_drugs[Reduce(intersect, list(index_row_drug_list[[drug_id_A_in]],
                                                 index_row_drug_list[[drug_id_B_in]],
                                                 index_row_drug_list[[drug_id_C_in]],
                                                 index_row_drug_list[[drug_id_D_in]],
                                                 index_row_drug_list[[drug_id_E_in]])),"exact_exposure_count"])
}

set_5s <- get_set_Ns(5, atleast_4s, min_exposure_count = 1e4) #<- test that
print(paste0("Generated set_5s, with ", dim(set_5s)[1], " rows, from 4s "))

print(" COUNTING 5s ... ")
atleast_5s <- set_5s %>%
  dplyr::rowwise() %>%
  dplyr::mutate(atleast_exposure_count = get_5s_count(drug_id_A, drug_id_B, drug_id_C, drug_id_D, drug_id_E))
save(atleast_5s, file = "truven_atleast_5s.RData")






# Assemble: Join atleast and exact and label with drug names: --------------------------------------------------------------------
atleast_1s <- atleast_1s %>% dplyr::arrange(-atleast_exposure_count)
atleast_2s <- atleast_2s %>% dplyr::arrange(-atleast_exposure_count)
atleast_3s <- atleast_3s %>% dplyr::arrange(-atleast_exposure_count)
atleast_4s <- atleast_4s %>% dplyr::arrange(-atleast_exposure_count)
atleast_5s <- atleast_5s %>% dplyr::arrange(-atleast_exposure_count)

exact_1s <- summary_exact_drugs %>% filter(n_drugs ==1)
exact_2s <- summary_exact_drugs %>% filter(n_drugs ==2)
exact_3s <- summary_exact_drugs %>% filter(n_drugs ==3)
exact_4s <- summary_exact_drugs %>% filter(n_drugs ==4)
exact_5s <- summary_exact_drugs %>% filter(n_drugs ==5)

# Separate the drug_id into columns
exact_1s <- exact_1s %>% 
  tidyr::separate(col=exact_drug_code, into=c("tmp","drug_id","tmp2"),"_", convert=TRUE) %>%
  dplyr::select(-tmp, -tmp2)
exact_2s <- exact_2s %>% 
  tidyr::separate(col=exact_drug_code, into=c("tmp","drug_id_A","drug_id_B","tmp2"),"_", convert=TRUE) %>%
  dplyr::select(-tmp, -tmp2)
exact_3s <- exact_3s %>% 
  tidyr::separate(col=exact_drug_code, into=c("tmp","drug_id_A","drug_id_B","drug_id_C","tmp2"),"_", convert=TRUE) %>%
  dplyr::select(-tmp, -tmp2)
exact_4s <- exact_4s %>% 
  tidyr::separate(col=exact_drug_code, into=c("tmp","drug_id_A","drug_id_B","drug_id_C","drug_id_D","tmp2"),"_", convert=TRUE) %>%
  dplyr::select(-tmp, -tmp2)
exact_5s <- exact_5s %>% 
  tidyr::separate(col=exact_drug_code, into=c("tmp","drug_id_A","drug_id_B","drug_id_C","drug_id_D","drug_id_E","tmp2"),"_", convert=TRUE) %>%
  dplyr::select(-tmp, -tmp2)

# Join the atleast and exact exposure counts into 5 tables:
# For 1s, include all 
drugs_1s <- atleast_1s %>%
  dplyr::left_join(exact_1s %>% dplyr::select(-n_drugs), by = "drug_id")
# For 2s, include all
drugs_2s <- atleast_2s %>%
  dplyr::left_join(exact_2s %>% dplyr::select(-n_drugs), by = c("drug_id_A","drug_id_B"))

drugs_3s <- atleast_3s %>%
  dplyr::filter(exposure_count > 1e4) %>%
  dplyr::left_join(exact_3s %>% dplyr::select(-n_drugs), by = c("drug_id_A","drug_id_B","drug_id_C"))

drugs_4s <- atleast_4s %>%
  dplyr::filter(exposure_count > 1e4) %>%
  dplyr::left_join(exact_4s %>% dplyr::select(-n_drugs), by = c("drug_id_A","drug_id_B","drug_id_C","drug_id_D"))

drugs_5s <- atleast_5s %>%
  dplyr::filter(exposure_count > 1e4) %>%
  dplyr::left_join(exact_5s %>% dplyr::select(-n_drugs), by = c("drug_id_A","drug_id_B","drug_id_C","drug_id_D","drug_id_E"))

# Join to drug names, instead of drug_id
drugs_1s <- drugs_1s %>%
  dplyr::left_join(Q_drug_id, by = "drug_id") %>%
  dplyr::select(drug_name_A = drug_concept_name, atleast_exposure_count, exact_exposure_count, estimate_drug_cost_per_day = drug_cost_per_day)

drugs_2s <- drugs_2s %>%
  dplyr::left_join(Q_drug_id, by = c("drug_id_A"="drug_id")) %>%
  dplyr::left_join(Q_drug_id, by = c("drug_id_B"="drug_id")) %>%
  dplyr::mutate(estimate_drug_cost_per_day = drug_cost_per_day.x + drug_cost_per_day.y) %>%
  dplyr::select(drug_name_A = drug_concept_name.x, drug_name_B = drug_concept_name.y, 
                atleast_exposure_count, exact_exposure_count, estimate_drug_cost_per_day)

drugs_3s <- drugs_3s %>%
  dplyr::left_join(Q_drug_id, by = c("drug_id_A"="drug_id")) %>%
  dplyr::left_join(Q_drug_id, by = c("drug_id_B"="drug_id")) %>%
  dplyr::left_join(Q_drug_id, by = c("drug_id_C"="drug_id")) %>%
  dplyr::mutate(estimate_drug_cost_per_day = drug_cost_per_day.x + drug_cost_per_day.y + drug_cost_per_day) %>%
  dplyr::select(drug_name_A = drug_concept_name.x, drug_name_B = drug_concept_name.y, drug_name_C = drug_concept_name,
                atleast_exposure_count, exact_exposure_count, estimate_drug_cost_per_day)

drugs_4s <- drugs_4s %>%
  dplyr::left_join(Q_drug_id, by = c("drug_id_A"="drug_id")) %>%
  dplyr::left_join(Q_drug_id, by = c("drug_id_B"="drug_id")) %>%
  dplyr::left_join(Q_drug_id, by = c("drug_id_C"="drug_id")) %>%
  dplyr::left_join(Q_drug_id, by = c("drug_id_D"="drug_id")) %>%
  dplyr::mutate(estimate_drug_cost_per_day = drug_cost_per_day.x + drug_cost_per_day.y + drug_cost_per_day.x.x + drug_cost_per_day.y.y) %>%
  dplyr::select(drug_name_A = drug_concept_name.x, drug_name_B = drug_concept_name.y, drug_name_C = drug_concept_name.x.x,
                drug_name_D = drug_concept_name.y.y,
                atleast_exposure_count, exact_exposure_count, estimate_drug_cost_per_day)

drugs_5s <- drugs_5s %>%
  dplyr::left_join(Q_drug_id, by = c("drug_id_A"="drug_id")) %>%
  dplyr::left_join(Q_drug_id, by = c("drug_id_B"="drug_id")) %>%
  dplyr::left_join(Q_drug_id, by = c("drug_id_C"="drug_id")) %>%
  dplyr::left_join(Q_drug_id, by = c("drug_id_D"="drug_id")) %>%
  dplyr::left_join(Q_drug_id, by = c("drug_id_E"="drug_id")) %>%
  dplyr::mutate(estimate_drug_cost_per_day = drug_cost_per_day.x + drug_cost_per_day.y + drug_cost_per_day.x.x + drug_cost_per_day.y.y
                + drug_cost_per_day) %>%
  dplyr::select(drug_name_A = drug_concept_name.x, drug_name_B = drug_concept_name.y, drug_name_C = drug_concept_name.x.x,
                drug_name_D = drug_concept_name.y.y, drug_name_E = drug_concept_name,
                atleast_exposure_count, exact_exposure_count, estimate_drug_cost_per_day)

# Assemble: Calculate fractions: --------------------------------------------------------------------------------
n_windows = sum(summary_exact_drugs$exact_exposure_count) # = 1.7e9
drugs_1s = drugs_1s %>%
  dplyr::mutate(fraction_exact = round(exact_exposure_count / atleast_exposure_count, 4),
                fraction_allwindows = round(atleast_exposure_count / n_windows, 4))
drugs_2s = drugs_2s %>%
  dplyr::mutate(fraction_exact = round(exact_exposure_count / atleast_exposure_count, 4),
                fraction_allwindows = round(atleast_exposure_count / n_windows, 4))
drugs_3s = drugs_3s %>%
  dplyr::mutate(fraction_exact = round(exact_exposure_count / atleast_exposure_count, 4),
                fraction_allwindows = round(atleast_exposure_count / n_windows, 4))
drugs_4s = drugs_4s %>%
  dplyr::mutate(fraction_exact = round(exact_exposure_count / atleast_exposure_count, 4),
                fraction_allwindows = round(atleast_exposure_count / n_windows, 4))
drugs_5s = drugs_5s %>%
  dplyr::mutate(fraction_exact = round(exact_exposure_count / atleast_exposure_count, 4),
                fraction_allwindows = round(atleast_exposure_count / n_windows, 4))

# Assemble: Overrepresentation: --------------------------------------------------------------------------------
atleast_2s = atleast_2s %>%
  dplyr::left_join(atleast_1s %>% dplyr::select(drug_id_A = drug_id, count_A = atleast_exposure_count), by = "drug_id_A") %>%
  dplyr::left_join(atleast_1s %>% dplyr::select(drug_id_B = drug_id, count_B = atleast_exposure_count), by = "drug_id_B") %>%
  dplyr::mutate(overrep = round(atleast_exposure_count / (count_A / n_windows * count_B), 3))

atleast_3s = atleast_3s %>%
  dplyr::left_join(atleast_1s %>% dplyr::select(drug_id_A = drug_id, count_A = atleast_exposure_count), by = "drug_id_A") %>%
  dplyr::left_join(atleast_1s %>% dplyr::select(drug_id_B = drug_id, count_B = atleast_exposure_count), by = "drug_id_B") %>%
  dplyr::left_join(atleast_1s %>% dplyr::select(drug_id_C = drug_id, count_C = atleast_exposure_count), by = "drug_id_C") %>%
  dplyr::left_join(atleast_2s %>% dplyr::select(drug_id_A = drug_id_A, drug_id_B = drug_id_B, count_AB = atleast_exposure_count), 
                   by = c("drug_id_A", "drug_id_B")) %>%
  dplyr::left_join(atleast_2s %>% dplyr::select(drug_id_A = drug_id_A, drug_id_C = drug_id_B, count_AC = atleast_exposure_count), 
                   by = c("drug_id_A", "drug_id_C")) %>%
  dplyr::left_join(atleast_2s %>% dplyr::select(drug_id_B = drug_id_A, drug_id_C = drug_id_B, count_BC = atleast_exposure_count), 
                   by = c("drug_id_B", "drug_id_C")) %>%
  dplyr::mutate(overrep_from_1 = round(atleast_exposure_count / (count_A / n_windows * count_B / n_windows * count_C), 3)) %>%
  dplyr::mutate(overrep_from_AB_C = round(atleast_exposure_count / (count_AB / n_windows * count_C), 3),
                overrep_from_AC_B = round(atleast_exposure_count / (count_AC / n_windows * count_B), 3),
                overrep_from_BC_A = round(atleast_exposure_count / (count_BC / n_windows * count_A), 3),
                overrep_min = pmin(overrep_from_AB_C, overrep_from_AC_B, overrep_from_BC_A))

atleast_4s = atleast_4s %>%
  dplyr::left_join(atleast_1s %>% dplyr::select(drug_id_A = drug_id, count_A = atleast_exposure_count), by = "drug_id_A") %>%
  dplyr::left_join(atleast_1s %>% dplyr::select(drug_id_B = drug_id, count_B = atleast_exposure_count), by = "drug_id_B") %>%
  dplyr::left_join(atleast_1s %>% dplyr::select(drug_id_C = drug_id, count_C = atleast_exposure_count), by = "drug_id_C") %>%
  dplyr::left_join(atleast_1s %>% dplyr::select(drug_id_D = drug_id, count_D = atleast_exposure_count), by = "drug_id_D") %>%
  dplyr::left_join(atleast_3s %>% 
                     dplyr::select(drug_id_A = drug_id_A, drug_id_B = drug_id_B, drug_id_C = drug_id_C, count_ABC = atleast_exposure_count), 
                   by = c("drug_id_A", "drug_id_B", "drug_id_C")) %>%
  dplyr::left_join(atleast_3s %>% 
                     dplyr::select(drug_id_A = drug_id_A, drug_id_B = drug_id_B, drug_id_D = drug_id_C, count_ABD = atleast_exposure_count), 
                   by = c("drug_id_A", "drug_id_B", "drug_id_D")) %>%
  dplyr::left_join(atleast_3s %>% 
                     dplyr::select(drug_id_A = drug_id_A, drug_id_C = drug_id_B, drug_id_D = drug_id_C, count_ACD = atleast_exposure_count), 
                   by = c("drug_id_A", "drug_id_C", "drug_id_D")) %>%
  dplyr::left_join(atleast_3s %>% 
                     dplyr::select(drug_id_B = drug_id_A, drug_id_C = drug_id_B, drug_id_D = drug_id_C, count_BCD = atleast_exposure_count), 
                   by = c("drug_id_B", "drug_id_C", "drug_id_D")) %>%
  dplyr::mutate(overrep_from_1 = round(atleast_exposure_count / (count_A / n_windows * count_B / n_windows * count_C / n_windows * count_D), 3),
                overrep_from_ABC_D = round(atleast_exposure_count / (count_ABC / n_windows * count_D), 3),
                overrep_from_ABD_C = round(atleast_exposure_count / (count_ABD / n_windows * count_C), 3),
                overrep_from_ACD_B = round(atleast_exposure_count / (count_ACD / n_windows * count_B), 3),
                overrep_from_BCD_A = round(atleast_exposure_count / (count_BCD / n_windows * count_A), 3),
                overrep_min = pmin(overrep_from_ABC_D, overrep_from_ABD_C, overrep_from_ACD_B, overrep_from_BCD_A))

atleast_5s = atleast_5s %>%
  dplyr::left_join(atleast_1s %>% dplyr::select(drug_id_A = drug_id, count_A = atleast_exposure_count), by = "drug_id_A") %>%
  dplyr::left_join(atleast_1s %>% dplyr::select(drug_id_B = drug_id, count_B = atleast_exposure_count), by = "drug_id_B") %>%
  dplyr::left_join(atleast_1s %>% dplyr::select(drug_id_C = drug_id, count_C = atleast_exposure_count), by = "drug_id_C") %>%
  dplyr::left_join(atleast_1s %>% dplyr::select(drug_id_D = drug_id, count_D = atleast_exposure_count), by = "drug_id_D") %>%
  dplyr::left_join(atleast_1s %>% dplyr::select(drug_id_E = drug_id, count_E = atleast_exposure_count), by = "drug_id_E") %>%
  dplyr::left_join(atleast_4s %>% 
                     dplyr::select(drug_id_A = drug_id_A, drug_id_B = drug_id_B, drug_id_C = drug_id_C, drug_id_D = drug_id_D, count_ABCD = atleast_exposure_count), 
                   by = c("drug_id_A", "drug_id_B", "drug_id_C", "drug_id_D")) %>%
  dplyr::left_join(atleast_4s %>% 
                     dplyr::select(drug_id_A = drug_id_A, drug_id_B = drug_id_B, drug_id_C = drug_id_C, drug_id_E = drug_id_D, count_ABCE = atleast_exposure_count), 
                   by = c("drug_id_A", "drug_id_B", "drug_id_C", "drug_id_E")) %>%
  dplyr::left_join(atleast_4s %>% 
                     dplyr::select(drug_id_A = drug_id_A, drug_id_B = drug_id_B, drug_id_D = drug_id_C, drug_id_E = drug_id_D, count_ABDE = atleast_exposure_count), 
                   by = c("drug_id_A", "drug_id_B", "drug_id_D", "drug_id_E")) %>%
  dplyr::left_join(atleast_4s %>% 
                     dplyr::select(drug_id_A = drug_id_A, drug_id_C = drug_id_B, drug_id_D = drug_id_C, drug_id_E = drug_id_D, count_ACDE = atleast_exposure_count), 
                   by = c("drug_id_A", "drug_id_C", "drug_id_D", "drug_id_E")) %>%
  dplyr::left_join(atleast_4s %>% 
                     dplyr::select(drug_id_B = drug_id_A, drug_id_C = drug_id_B, drug_id_D = drug_id_C, drug_id_E = drug_id_D, count_BCDE = atleast_exposure_count), 
                   by = c("drug_id_B", "drug_id_C", "drug_id_D", "drug_id_E")) %>%
  dplyr::mutate(overrep_from_1 = round(atleast_exposure_count / (count_A / n_windows * count_B / n_windows * count_C / n_windows * count_D / n_windows * count_E), 3),
                overrep_from_ABCD_E = round(atleast_exposure_count / (count_ABCD / n_windows * count_E), 3),
                overrep_from_ABCE_D = round(atleast_exposure_count / (count_ABCE / n_windows * count_D), 3),
                overrep_from_ABDE_C = round(atleast_exposure_count / (count_ABDE / n_windows * count_C), 3),
                overrep_from_ACDE_B = round(atleast_exposure_count / (count_ACDE / n_windows * count_B), 3),
                overrep_from_BCDE_A = round(atleast_exposure_count / (count_BCDE / n_windows * count_A), 3),
                overrep_min = pmin(overrep_from_ABCD_E, overrep_from_ABCE_D, overrep_from_ABDE_C, overrep_from_ACDE_B, overrep_from_BCDE_A))

atleast_2s = atleast_2s %>% dplyr::select(drug_id_A, drug_id_B, atleast_exposure_count, overrep_1 = overrep, overrep_N = overrep)
atleast_3s = atleast_3s %>% dplyr::select(drug_id_A, drug_id_B, drug_id_C, atleast_exposure_count, overrep_1 = overrep_from_1, overrep_N = overrep_min)
atleast_4s = atleast_4s %>% dplyr::select(drug_id_A, drug_id_B, drug_id_C, drug_id_D, atleast_exposure_count, overrep_1 = overrep_from_1, overrep_N = overrep_min)
atleast_5s = atleast_5s %>% dplyr::select(drug_id_A, drug_id_B, drug_id_C, drug_id_D, drug_id_E, atleast_exposure_count, overrep_1 = overrep_from_1, overrep_N = overrep_min)

# Join to drugs_2s, _3s, _4s, _5s:
atleast_2s = atleast_2s %>% dplyr::rename(observe_per_expect_N1 = overrep_N)
atleast_3s = atleast_3s %>% dplyr::rename(observe_per_expect_1s = overrep_1, observe_per_expect_N1 = overrep_N)
atleast_4s = atleast_4s %>% dplyr::rename(observe_per_expect_1s = overrep_1, observe_per_expect_N1 = overrep_N)
atleast_5s = atleast_5s %>% dplyr::rename(observe_per_expect_1s = overrep_1, observe_per_expect_N1 = overrep_N)

Q_drug_name = Q_drug_id %>% dplyr::select(drug_id, drug_name = drug_concept_name)
atleast_2s = atleast_2s %>%
  dplyr::left_join(Q_drug_name %>% dplyr::select(drug_id_A = drug_id, drug_name_A = drug_name), by = "drug_id_A") %>%
  dplyr::left_join(Q_drug_name %>% dplyr::select(drug_id_B = drug_id, drug_name_B = drug_name), by = "drug_id_B") %>%
  dplyr::select(drug_name_A, drug_name_B, observe_per_expect_N1)
atleast_3s = atleast_3s %>%
  dplyr::left_join(Q_drug_name %>% dplyr::select(drug_id_A = drug_id, drug_name_A = drug_name), by = "drug_id_A") %>%
  dplyr::left_join(Q_drug_name %>% dplyr::select(drug_id_B = drug_id, drug_name_B = drug_name), by = "drug_id_B") %>%
  dplyr::left_join(Q_drug_name %>% dplyr::select(drug_id_C = drug_id, drug_name_C = drug_name), by = "drug_id_C") %>%
  dplyr::select(drug_name_A, drug_name_B, drug_name_C, observe_per_expect_1s, observe_per_expect_N1)
atleast_4s = atleast_4s %>%
  dplyr::left_join(Q_drug_name %>% dplyr::select(drug_id_A = drug_id, drug_name_A = drug_name), by = "drug_id_A") %>%
  dplyr::left_join(Q_drug_name %>% dplyr::select(drug_id_B = drug_id, drug_name_B = drug_name), by = "drug_id_B") %>%
  dplyr::left_join(Q_drug_name %>% dplyr::select(drug_id_C = drug_id, drug_name_C = drug_name), by = "drug_id_C") %>%
  dplyr::left_join(Q_drug_name %>% dplyr::select(drug_id_D = drug_id, drug_name_D = drug_name), by = "drug_id_D") %>%
  dplyr::select(drug_name_A, drug_name_B, drug_name_C, drug_name_D, observe_per_expect_1s, observe_per_expect_N1)
atleast_5s = atleast_5s %>%
  dplyr::left_join(Q_drug_name %>% dplyr::select(drug_id_A = drug_id, drug_name_A = drug_name), by = "drug_id_A") %>%
  dplyr::left_join(Q_drug_name %>% dplyr::select(drug_id_B = drug_id, drug_name_B = drug_name), by = "drug_id_B") %>%
  dplyr::left_join(Q_drug_name %>% dplyr::select(drug_id_C = drug_id, drug_name_C = drug_name), by = "drug_id_C") %>%
  dplyr::left_join(Q_drug_name %>% dplyr::select(drug_id_D = drug_id, drug_name_D = drug_name), by = "drug_id_D") %>%
  dplyr::left_join(Q_drug_name %>% dplyr::select(drug_id_E = drug_id, drug_name_E = drug_name), by = "drug_id_E") %>%
  dplyr::select(drug_name_A, drug_name_B, drug_name_C, drug_name_D, drug_name_E, observe_per_expect_1s, observe_per_expect_N1)

drugs_2s = drugs_2s %>%
  dplyr::left_join(atleast_2s, by = c("drug_name_A", "drug_name_B"))
drugs_3s = drugs_3s %>%
  dplyr::left_join(atleast_3s, by = c("drug_name_A", "drug_name_B", "drug_name_C"))
drugs_4s = drugs_4s %>%
  dplyr::left_join(atleast_4s, by = c("drug_name_A", "drug_name_B", "drug_name_C", "drug_name_D"))
drugs_5s = drugs_5s %>%
  dplyr::left_join(atleast_5s, by = c("drug_name_A", "drug_name_B", "drug_name_C", "drug_name_D", "drug_name_E"))



# Assemble: Hide counts <100: --------------------------------------------------------------------------------
# if exact_exposure_count < 100 -> (i) "<100" and (ii) fraction_exact -> ">[value if exact=100]
drugs_2s = drugs_2s %>%
  dplyr::mutate(observe_per_expect_N1 = ifelse(atleast_exposure_count < 100, " ", observe_per_expect_N1),
                atleast_exposure_count = ifelse(atleast_exposure_count < 100, "<100", atleast_exposure_count),
                fraction_exact = ifelse(is.na(exact_exposure_count) | exact_exposure_count < 100, " ", fraction_exact),
                fraction_allwindows = ifelse(atleast_exposure_count < 100, " ", fraction_allwindows),
                exact_exposure_count = ifelse(is.na(exact_exposure_count), 0, exact_exposure_count),
                exact_exposure_count = ifelse(exact_exposure_count < 100, "<100", exact_exposure_count))

drugs_3s = drugs_3s %>%
  dplyr::mutate(observe_per_expect_N1 = ifelse(atleast_exposure_count < 100, " ", observe_per_expect_N1),
                observe_per_expect_1s = ifelse(atleast_exposure_count < 100, " ", observe_per_expect_1s),
                atleast_exposure_count = ifelse(atleast_exposure_count < 100, "<100", atleast_exposure_count),
                fraction_exact = ifelse(is.na(exact_exposure_count) | exact_exposure_count < 100, " ", fraction_exact),
                fraction_allwindows = ifelse(atleast_exposure_count < 100, " ", fraction_allwindows),
                exact_exposure_count = ifelse(is.na(exact_exposure_count), 0, exact_exposure_count),
                exact_exposure_count = ifelse(exact_exposure_count < 100, "<100", exact_exposure_count))

drugs_4s = drugs_4s %>%
  dplyr::mutate(observe_per_expect_N1 = ifelse(atleast_exposure_count < 100, " ", observe_per_expect_N1),
                observe_per_expect_1s = ifelse(atleast_exposure_count < 100, " ", observe_per_expect_1s),
                atleast_exposure_count = ifelse(atleast_exposure_count < 100, "<100", atleast_exposure_count),
                fraction_exact = ifelse(is.na(exact_exposure_count) | exact_exposure_count < 100, " ", fraction_exact),
                fraction_allwindows = ifelse(atleast_exposure_count < 100, " ", fraction_allwindows),
                exact_exposure_count = ifelse(is.na(exact_exposure_count), 0, exact_exposure_count),
                exact_exposure_count = ifelse(exact_exposure_count < 100, "<100", exact_exposure_count))

drugs_5s = drugs_5s %>%
  dplyr::mutate(observe_per_expect_N1 = ifelse(atleast_exposure_count < 100, " ", observe_per_expect_N1),
                observe_per_expect_1s = ifelse(atleast_exposure_count < 100, " ", observe_per_expect_1s),
                atleast_exposure_count = ifelse(atleast_exposure_count < 100, "<100", atleast_exposure_count),
                fraction_exact = ifelse(is.na(exact_exposure_count) | exact_exposure_count < 100, " ", fraction_exact),
                fraction_allwindows = ifelse(atleast_exposure_count < 100, " ", fraction_allwindows),
                exact_exposure_count = ifelse(is.na(exact_exposure_count), 0, exact_exposure_count),
                exact_exposure_count = ifelse(exact_exposure_count < 100, "<100", exact_exposure_count))


# Assemble: Save: --------------------------------------------------------------------------------
save(drugs_1s, drugs_2s, drugs_3s, drugs_4s, drugs_5s, file="db_drugs_12345.RData")


write.table(drugs_1s, file="db_final/db_drugs_1s.tsv", quote=FALSE, row.name = FALSE, sep = "\t")
write.table(drugs_2s, file="db_final/db_drugs_2s.tsv", quote=FALSE, row.name = FALSE, sep = "\t")
write.table(drugs_3s, file="db_final/db_drugs_3s.tsv", quote=FALSE, row.name = FALSE, sep = "\t")
write.table(drugs_4s, file="db_final/db_drugs_4s.tsv", quote=FALSE, row.name = FALSE, sep = "\t")
write.table(drugs_5s, file="db_final/db_drugs_5s.tsv", quote=FALSE, row.name = FALSE, sep = "\t")





