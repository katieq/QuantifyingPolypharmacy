# Quinn & Shah, A dataset quantifying polypharmacy in the United States (2017)
# Code to generate Data Record 2: Drug class combinations
# This file follows on, and uses variables generated in the workspace of the file for Data Record 1 (multidrug_db.R), 
# KQuinn 2016-2017

# Convert exact-exposures for drug ingredients to classes -------------------------

get_exact_class_code <- function(class_id_vec){
  class_id_vec <- sort(class_id_vec)
  exact_drug_code <- paste0("_",paste(class_id_vec, collapse = "_"),"_")
}

Q_drug_id$atc12 <- substr(Q_drug_id$atc,1,3)

tmp <- summary_exact_drugs %>%
  dplyr::mutate(rn = row_number())

summary_exact_class_list = list()
for (i in 1:max(summary_exact_drugs$n_drugs)){
  print(paste0("Class codes for n_drugs = ",i," ..."))
  summary_exact_class_list[[i]] <- tmp %>%
    dplyr::filter(n_drugs == i) %>%
    tidyr::separate(exact_drug_code, as.character(c(0:(i+1))), 
                    sep="_", extra="warn",fill="warn",remove=TRUE, convert=TRUE) %>%
    tidyr::gather(row_pos, drug_id, 1+c(1:i), na.rm=TRUE) %>%
    dplyr::select(rn, drug_id, exact_exposure_count) %>%
    dplyr::left_join(Q_drug_id %>% dplyr::select(drug_id, atc12), by = "drug_id") %>%
    dplyr::group_by(rn) %>%
    dplyr::summarise(exact_class_code = get_exact_class_code(atc12),
                     n_drugs = i,
                     exact_exposure_count = exact_exposure_count[1])
}
remove(tmp)
summary_exact_class <- as.data.frame(rbindlist(summary_exact_class_list))

summary_exact_class <- summary_exact_class %>%
  dplyr::group_by(exact_class_code) %>%
  dplyr::summarise(exact_exposure_count = sum(exact_exposure_count),
                   n_drugs = n_drugs[1]) %>%
  dplyr::arrange(-exact_exposure_count)



# Count atleast-exposure for N-class combinations ------------------------

# Create an index of the rows containing n-drug combinations
index_row_drug_ndrug_list = list()
print("Creating row numbers.. ")
tmp <- summ_exact_class %>%
  dplyr::mutate(rn = row_number())
for (i in 1:max(summary_exact_class$n_drugs)){
  print(paste0("Index for n_drugs = ",i," ..."))
  summ_exact_drugs_i <- tmp %>% dplyr::filter(n_drugs == i)
  index_row_drug_ndrug_list[[i]] <- summ_exact_drugs_i %>%
    tidyr::separate(exact_class_code, as.character(c(0:(i+1))), 
                    sep="_", extra="warn",fill="warn",remove=TRUE, convert=TRUE) %>%
    tidyr::gather(row_pos, atc_class, 1+c(1:i), na.rm=TRUE) %>%
    dplyr::select(rn, atc_class, exact_exposure_count)
}
remove(tmp)
index_row_drug <- as.data.frame(rbindlist(index_row_drug_ndrug_list))
remove(index_row_drug_ndrug_list, i)

# Count single-drug atleast exposures
atleast_1s <- index_row_drug %>%
  dplyr::group_by(atc_class) %>%
  dplyr::summarise(atleast_exposure_count = sum(exact_exposure_count)) %>%
  dplyr::arrange(-atleast_exposure_count)

save(atleast_1s, file = "class_truven_atleast_1s.RData")

# Create index of the location of drug_id, used to efficiently count N-drug atleast exposures
index_row_class_list <- list()
for (i in 1:dim(drug_class_12)[1]){
  index_row_class_list[[i]] <- (index_row_drug %>% dplyr::semi_join(data.frame(atc_class = drug_class_12$atc12[i]), by = "atc_class"))$rn
}
# save(index_row_class_list, file = "../../../../scratch/users/kjquinn/drug_combo_db/full_class_truven_index_row_drug_list.RData")

drug_class_12 <- drug_class_12 %>%
  dplyr::mutate(rn = row_number())

# This function creates set of N-drug-class combinations, 
# where each N N-1 drug class combination has at lest min_exposure_count exposures.
source("code/get_set_Ns.R")

get_2s_count <- function(drug_id_A_in, drug_id_B_in){
  sum(summary_exact_class[intersect(index_row_class_list[[drug_id_A_in]],
                                    index_row_class_list[[drug_id_B_in]]),"exact_exposure_count"])
}

print(" COUNTING 2s ... ")
set_2s <- get_set_Ns(2, atleast_1s %>% dplyr::rename(drug_id = atc_class), min_exposure_count = 1e3) %>%
  dplyr::rename(atc_A = drug_id_A, atc_B = drug_id_B)
print(paste0("Generated set_2s, with ", dim(set_2s)[1], " rows"))


atleast_2s <- set_2s %>% # columns need to be drug_id_A and drug_id_B
  dplyr::left_join(drug_class_12 %>% dplyr::select(atc12, drug_id_A = rn), by = c("atc_A"="atc12")) %>%
  dplyr::left_join(drug_class_12 %>% dplyr::select(atc12, drug_id_B = rn), by = c("atc_B"="atc12")) %>%
  dplyr::rowwise() %>%
  dplyr::mutate(atleast_exposure_count = get_2s_count(drug_id_A, drug_id_B)) %>%
  dplyr::select(atc_A, atc_B, atleast_exposure_count)
save(atleast_2s, file = "class_truven_atleast_2s.RData")

# _3s
get_3s_count <- function(drug_id_A_in, drug_id_B_in, drug_id_C_in){
  sum(summary_exact_class[Reduce(intersect, list(index_row_class_list[[drug_id_A_in]],
                                                 index_row_class_list[[drug_id_B_in]],
                                                 index_row_class_list[[drug_id_C_in]])),"exact_exposure_count"])
}

set_3s <- get_set_Ns(3, atleast_2s %>% dplyr::rename(drug_id_A = atc_A, drug_id_B = atc_B), min_exposure_count = 1e3) %>%
  dplyr::rename(atc_A = drug_id_A, atc_B = drug_id_B, atc_C = drug_id_C)
print(paste0("Generated set_3s, with ", dim(set_3s)[1], " rows, from 2s " ))

print(" COUNTING 3s ... ")
atleast_3s <- set_3s %>%
  dplyr::left_join(drug_class_12 %>% dplyr::select(atc12, drug_id_A = rn), by = c("atc_A"="atc12")) %>%
  dplyr::left_join(drug_class_12 %>% dplyr::select(atc12, drug_id_B = rn), by = c("atc_B"="atc12")) %>%
  dplyr::left_join(drug_class_12 %>% dplyr::select(atc12, drug_id_C = rn), by = c("atc_C"="atc12")) %>%
  dplyr::rowwise() %>%
  dplyr::mutate(atleast_exposure_count = get_3s_count(drug_id_A, drug_id_B, drug_id_C)) %>%
  dplyr::select(atc_A, atc_B, atc_C, atleast_exposure_count)
save(atleast_3s, file = "class_truven_atleast_3s.RData")

# _4s
get_4s_count <- function(drug_id_A_in, drug_id_B_in, drug_id_C_in, drug_id_D_in){
  sum(summary_exact_class[Reduce(intersect, list(index_row_class_list[[drug_id_A_in]],
                                                 index_row_class_list[[drug_id_B_in]],
                                                 index_row_class_list[[drug_id_C_in]],
                                                 index_row_class_list[[drug_id_D_in]])),"exact_exposure_count"])
}

set_4s <- get_set_Ns(4, atleast_3s %>% dplyr::rename(drug_id_A = atc_A, drug_id_B = atc_B, drug_id_C = atc_C), min_exposure_count = 5e4) %>% # prev 1e3
  dplyr::rename(atc_A = drug_id_A, atc_B = drug_id_B, atc_C = drug_id_C, atc_D = drug_id_D)
print(paste0("Generated set_4s, with ", dim(set_4s)[1], " rows, from 3s " ))

print(" COUNTING 4s ... ")
atleast_4s <- set_4s %>%
  dplyr::left_join(drug_class_12 %>% dplyr::select(atc12, drug_id_A = rn), by = c("atc_A"="atc12")) %>%
  dplyr::left_join(drug_class_12 %>% dplyr::select(atc12, drug_id_B = rn), by = c("atc_B"="atc12")) %>%
  dplyr::left_join(drug_class_12 %>% dplyr::select(atc12, drug_id_C = rn), by = c("atc_C"="atc12")) %>%
  dplyr::left_join(drug_class_12 %>% dplyr::select(atc12, drug_id_D = rn), by = c("atc_D"="atc12")) %>%
  dplyr::rowwise() %>%
  dplyr::mutate(atleast_exposure_count = get_4s_count(drug_id_A, drug_id_B, drug_id_C, drug_id_D)) %>%
  dplyr::select(atc_A, atc_B, atc_C, atc_D, atleast_exposure_count)
save(atleast_4s, file = "class_truven_atleast_4s.RData")

# _5s
get_5s_count <- function(drug_id_A_in, drug_id_B_in, drug_id_C_in, drug_id_D_in, drug_id_E_in){
  sum(summary_exact_class[Reduce(intersect, list(index_row_class_list[[drug_id_A_in]],
                                                 index_row_class_list[[drug_id_B_in]],
                                                 index_row_class_list[[drug_id_C_in]],
                                                 index_row_class_list[[drug_id_D_in]],
                                                 index_row_class_list[[drug_id_E_in]])),"exact_exposure_count"])
}

set_5s <- get_set_Ns(5, atleast_4s %>% 
                       dplyr::rename(drug_id_A = atc_A, drug_id_B = atc_B, drug_id_C = atc_C, drug_id_D = atc_D), 
                     min_exposure_count = 5e4) %>% # prev 1e3, then 5e4, then 5e5 for Cs
  dplyr::rename(atc_A = drug_id_A, atc_B = drug_id_B, atc_C = drug_id_C, atc_D = drug_id_D, atc_E = drug_id_E)
print(paste0("Generated set_5s, with ", dim(set_5s)[1], " rows, from 4s "))

print(" COUNTING 5s ... ")
atleast_5s <- set_5s %>%
  dplyr::left_join(drug_class_12 %>% dplyr::select(atc12, drug_id_A = rn), by = c("atc_A"="atc12")) %>%
  dplyr::left_join(drug_class_12 %>% dplyr::select(atc12, drug_id_B = rn), by = c("atc_B"="atc12")) %>%
  dplyr::left_join(drug_class_12 %>% dplyr::select(atc12, drug_id_C = rn), by = c("atc_C"="atc12")) %>%
  dplyr::left_join(drug_class_12 %>% dplyr::select(atc12, drug_id_D = rn), by = c("atc_D"="atc12")) %>%
  dplyr::left_join(drug_class_12 %>% dplyr::select(atc12, drug_id_E = rn), by = c("atc_E"="atc12")) %>%
  dplyr::rowwise() %>%
  dplyr::mutate(atleast_exposure_count = get_5s_count(drug_id_A, drug_id_B, drug_id_C, drug_id_D, drug_id_E)) %>%
  dplyr::select(atc_A, atc_B, atc_C, atc_D, atc_E, atleast_exposure_count)
save(atleast_5s, file = "class_truven_atleast_5s.RData")





# Assemble: Join atleast and exact and label with drug names: ---------------------------

atleast_1s <- atleast_1s %>% dplyr::arrange(-atleast_exposure_count)
atleast_2s <- atleast_2s %>% dplyr::arrange(-atleast_exposure_count)
atleast_3s <- atleast_3s %>% dplyr::arrange(-atleast_exposure_count)
atleast_4s <- atleast_4s %>% dplyr::arrange(-atleast_exposure_count)
atleast_5s <- atleast_5s %>% dplyr::arrange(-atleast_exposure_count)

exact_1s <- summ_exact_class %>% filter(n_drugs ==1)
exact_2s <- summ_exact_class %>% filter(n_drugs ==2)
exact_3s <- summ_exact_class %>% filter(n_drugs ==3)
exact_4s <- summ_exact_class %>% filter(n_drugs ==4)
exact_5s <- summ_exact_class %>% filter(n_drugs ==5)

exact_1s <- exact_1s %>% 
  tidyr::separate(col=exact_class_code, into=c("tmp","atc","tmp2"),"_", convert=TRUE) %>%
  dplyr::select(-tmp, -tmp2)
exact_2s <- exact_2s %>% 
  tidyr::separate(col=exact_class_code, into=c("tmp","atc_A","atc_B","tmp2"),"_", convert=TRUE) %>%
  dplyr::select(-tmp, -tmp2)
exact_3s <- exact_3s %>% 
  tidyr::separate(col=exact_class_code, into=c("tmp","atc_A","atc_B","atc_C","tmp2"),"_", convert=TRUE) %>%
  dplyr::select(-tmp, -tmp2)
exact_4s <- exact_4s %>% 
  tidyr::separate(col=exact_class_code, into=c("tmp","atc_A","atc_B","atc_C","atc_D","tmp2"),"_", convert=TRUE) %>%
  dplyr::select(-tmp, -tmp2)
exact_5s <- exact_5s %>% 
  tidyr::separate(col=exact_class_code, into=c("tmp","atc_A","atc_B","atc_C","atc_D","atc_E","tmp2"),"_", convert=TRUE) %>%
  dplyr::select(-tmp, -tmp2)

atc_classes_1s <- atleast_1s %>%
  dplyr::left_join(exact_1s %>% dplyr::select(-n_drugs), by = "atc")

atc_classes_2s <- atleast_2s %>%
  dplyr::left_join(exact_2s %>% dplyr::select(-n_drugs), by = c("atc_A","atc_B"))

atc_classes_3s <- atleast_3s %>%
  dplyr::filter(exposure_count > 1e4) %>%
  dplyr::left_join(exact_3s %>% dplyr::select(-n_drugs), by = c("atc_A","atc_B","atc_C"))

atc_classes_4s <- atleast_4s %>%
  dplyr::filter(exposure_count > 1e4) %>%
  dplyr::left_join(exact_4s %>% dplyr::select(-n_drugs), by = c("atc_A","atc_B","atc_C","atc_D"))

atc_classes_5s <- atleast_5s %>%
  dplyr::filter(exposure_count > 1e4) %>%
  dplyr::left_join(exact_5s %>% dplyr::select(-n_drugs), by = c("atc_A","atc_B","atc_C","atc_D","atc_E"))

drug_class_12 <- atc_classes %>% 
  dplyr::select(atc = atc_class, atc_class_name)

atc_classes_1s <- atc_classes_1s %>%
  dplyr::left_join(drug_class_12, by = "atc")
  dplyr::select(atc_A = atc, atc_class_name_A = atc_class_name, atleast_exposure_count, exact_exposure_count)

atc_classes_2s <- atc_classes_2s %>%
  dplyr::left_join(drug_class_12, by = c("atc_A"="atc")) %>%
  dplyr::left_join(drug_class_12, by = c("atc_B"="atc")) %>%
  dplyr::select(atc_A, atc_class_name_A = atc_class_name.x, atc_B, atc_class_name_B = atc_class_name.y, 
                atleast_exposure_count, exact_exposure_count)

atc_classes_3s <- atc_classes_3s %>%
  dplyr::left_join(drug_class_12, by = c("atc_A"="atc")) %>%
  dplyr::left_join(drug_class_12, by = c("atc_B"="atc")) %>%
  dplyr::left_join(drug_class_12, by = c("atc_C"="atc")) %>%
  dplyr::select(atc_A, atc_class_name_A = atc_class_name.x, 
                atc_B, atc_class_name_B = atc_class_name.y, 
                atc_C, atc_class_name_C = atc_class_name,
                atleast_exposure_count, exact_exposure_count)

atc_classes_4s <- atc_classes_4s %>%
  dplyr::left_join(drug_class_12, by = c("atc_A"="atc")) %>%
  dplyr::left_join(drug_class_12, by = c("atc_B"="atc")) %>%
  dplyr::left_join(drug_class_12, by = c("atc_C"="atc")) %>%
  dplyr::left_join(drug_class_12, by = c("atc_D"="atc")) %>%
  dplyr::select(atc_A, atc_class_name_A = atc_class_name.x, 
                atc_B, atc_class_name_B = atc_class_name.y, 
                atc_C, atc_class_name_C = atc_class_name.x.x,
                atc_D, atc_class_name_D = atc_class_name.y.y,
                atleast_exposure_count, exact_exposure_count)

atc_classes_5s <- atc_classes_5s %>%
  dplyr::left_join(drug_class_12, by = c("atc_A"="atc")) %>%
  dplyr::left_join(drug_class_12, by = c("atc_B"="atc")) %>%
  dplyr::left_join(drug_class_12, by = c("atc_C"="atc")) %>%
  dplyr::left_join(drug_class_12, by = c("atc_D"="atc")) %>%
  dplyr::left_join(drug_class_12, by = c("atc_E"="atc")) %>%
  dplyr::select(atc_A, atc_class_name_A = atc_class_name.x, 
                atc_B, atc_class_name_B = atc_class_name.y, 
                atc_C, atc_class_name_C = atc_class_name.x.x,
                atc_D, atc_class_name_D = atc_class_name.y.y,
                atc_E, atc_class_name_E = atc_class_name,
                atleast_exposure_count, exact_exposure_count)

# Assemble: Calculate fractions: -----------------------------------------------------
n_windows = sum(summary_exact_class$exact_exposure_count) # = 1.7e9
atc_classes_1s = atc_classes_1s %>%
  dplyr::mutate(fraction_exact = round(exact_exposure_count / atleast_exposure_count, 4),
                fraction_allwindows = round(atleast_exposure_count / n_windows, 4))
atc_classes_2s = atc_classes_2s %>%
  dplyr::mutate(fraction_exact = round(exact_exposure_count / atleast_exposure_count, 4),
                fraction_allwindows = round(atleast_exposure_count / n_windows, 4))
atc_classes_3s = atc_classes_3s %>%
  dplyr::mutate(fraction_exact = round(exact_exposure_count / atleast_exposure_count, 4),
                fraction_allwindows = round(atleast_exposure_count / n_windows, 4))
atc_classes_4s = atc_classes_4s %>%
  dplyr::mutate(fraction_exact = round(exact_exposure_count / atleast_exposure_count, 4),
                fraction_allwindows = round(atleast_exposure_count / n_windows, 4))
atc_classes_5s = atc_classes_5s %>%
  dplyr::mutate(fraction_exact = round(exact_exposure_count / atleast_exposure_count, 4),
                fraction_allwindows = round(atleast_exposure_count / n_windows, 4))

# Assemble: Calculate overrepresentation: -----------------------------------------------------

atleast_2s = atleast_2s %>%
  dplyr::left_join(atleast_1s %>% dplyr::select(atc_A = atc_class, count_A = atleast_exposure_count), by = "atc_A") %>%
  dplyr::left_join(atleast_1s %>% dplyr::select(atc_B = atc_class, count_B = atleast_exposure_count), by = "atc_B") %>%
  dplyr::mutate(overrep = round(atleast_exposure_count / (count_A / n_windows * count_B), 3))

atleast_3s = atleast_3s %>%
  dplyr::left_join(atleast_1s %>% dplyr::select(atc_A = atc_class, count_A = atleast_exposure_count), by = "atc_A") %>%
  dplyr::left_join(atleast_1s %>% dplyr::select(atc_B = atc_class, count_B = atleast_exposure_count), by = "atc_B") %>%
  dplyr::left_join(atleast_1s %>% dplyr::select(atc_C = atc_class, count_C = atleast_exposure_count), by = "atc_C") %>%
  dplyr::left_join(atleast_2s %>% dplyr::select(atc_A = atc_A, atc_B = atc_B, count_AB = atleast_exposure_count), 
                   by = c("atc_A", "atc_B")) %>%
  dplyr::left_join(atleast_2s %>% dplyr::select(atc_A = atc_A, atc_C = atc_B, count_AC = atleast_exposure_count), 
                   by = c("atc_A", "atc_C")) %>%
  dplyr::left_join(atleast_2s %>% dplyr::select(atc_B = atc_A, atc_C = atc_B, count_BC = atleast_exposure_count), 
                   by = c("atc_B", "atc_C")) %>%
  dplyr::mutate(overrep_from_1 = round(atleast_exposure_count / (count_A / n_windows * count_B / n_windows * count_C), 3)) %>%
  dplyr::mutate(overrep_from_AB_C = round(atleast_exposure_count / (count_AB / n_windows * count_C), 3),
                overrep_from_AC_B = round(atleast_exposure_count / (count_AC / n_windows * count_B), 3),
                overrep_from_BC_A = round(atleast_exposure_count / (count_BC / n_windows * count_A), 3),
                overrep_min = pmin(overrep_from_AB_C, overrep_from_AC_B, overrep_from_BC_A))

atleast_4s = atleast_4s %>%
  dplyr::left_join(atleast_1s %>% dplyr::select(atc_A = atc_class, count_A = atleast_exposure_count), by = "atc_A") %>%
  dplyr::left_join(atleast_1s %>% dplyr::select(atc_B = atc_class, count_B = atleast_exposure_count), by = "atc_B") %>%
  dplyr::left_join(atleast_1s %>% dplyr::select(atc_C = atc_class, count_C = atleast_exposure_count), by = "atc_C") %>%
  dplyr::left_join(atleast_1s %>% dplyr::select(atc_D = atc_class, count_D = atleast_exposure_count), by = "atc_D") %>%
  dplyr::left_join(atleast_3s %>% 
                     dplyr::select(atc_A = atc_A, atc_B = atc_B, atc_C = atc_C, count_ABC = atleast_exposure_count), 
                   by = c("atc_A", "atc_B", "atc_C")) %>%
  dplyr::left_join(atleast_3s %>% 
                     dplyr::select(atc_A = atc_A, atc_B = atc_B, atc_D = atc_C, count_ABD = atleast_exposure_count), 
                   by = c("atc_A", "atc_B", "atc_D")) %>%
  dplyr::left_join(atleast_3s %>% 
                     dplyr::select(atc_A = atc_A, atc_C = atc_B, atc_D = atc_C, count_ACD = atleast_exposure_count), 
                   by = c("atc_A", "atc_C", "atc_D")) %>%
  dplyr::left_join(atleast_3s %>% 
                     dplyr::select(atc_B = atc_A, atc_C = atc_B, atc_D = atc_C, count_BCD = atleast_exposure_count), 
                   by = c("atc_B", "atc_C", "atc_D")) %>%
  dplyr::mutate(overrep_from_1 = round(atleast_exposure_count / (count_A / n_windows * count_B / n_windows * count_C / n_windows * count_D), 3),
                overrep_from_ABC_D = round(atleast_exposure_count / (count_ABC / n_windows * count_D), 3),
                overrep_from_ABD_C = round(atleast_exposure_count / (count_ABD / n_windows * count_C), 3),
                overrep_from_ACD_B = round(atleast_exposure_count / (count_ACD / n_windows * count_B), 3),
                overrep_from_BCD_A = round(atleast_exposure_count / (count_BCD / n_windows * count_A), 3),
                overrep_min = pmin(overrep_from_ABC_D, overrep_from_ABD_C, overrep_from_ACD_B, overrep_from_BCD_A))

atleast_5s = atleast_5s %>%
  dplyr::left_join(atleast_1s %>% dplyr::select(atc_A = atc_class, count_A = atleast_exposure_count), by = "atc_A") %>%
  dplyr::left_join(atleast_1s %>% dplyr::select(atc_B = atc_class, count_B = atleast_exposure_count), by = "atc_B") %>%
  dplyr::left_join(atleast_1s %>% dplyr::select(atc_C = atc_class, count_C = atleast_exposure_count), by = "atc_C") %>%
  dplyr::left_join(atleast_1s %>% dplyr::select(atc_D = atc_class, count_D = atleast_exposure_count), by = "atc_D") %>%
  dplyr::left_join(atleast_1s %>% dplyr::select(atc_E = atc_class, count_E = atleast_exposure_count), by = "atc_E") %>%
  dplyr::left_join(atleast_4s %>% 
                     dplyr::select(atc_A = atc_A, atc_B = atc_B, atc_C = atc_C, atc_D = atc_D, count_ABCD = atleast_exposure_count), 
                   by = c("atc_A", "atc_B", "atc_C", "atc_D")) %>%
  dplyr::left_join(atleast_4s %>% 
                     dplyr::select(atc_A = atc_A, atc_B = atc_B, atc_C = atc_C, atc_E = atc_D, count_ABCE = atleast_exposure_count), 
                   by = c("atc_A", "atc_B", "atc_C", "atc_E")) %>%
  dplyr::left_join(atleast_4s %>% 
                     dplyr::select(atc_A = atc_A, atc_B = atc_B, atc_D = atc_C, atc_E = atc_D, count_ABDE = atleast_exposure_count), 
                   by = c("atc_A", "atc_B", "atc_D", "atc_E")) %>%
  dplyr::left_join(atleast_4s %>% 
                     dplyr::select(atc_A = atc_A, atc_C = atc_B, atc_D = atc_C, atc_E = atc_D, count_ACDE = atleast_exposure_count), 
                   by = c("atc_A", "atc_C", "atc_D", "atc_E")) %>%
  dplyr::left_join(atleast_4s %>% 
                     dplyr::select(atc_B = atc_A, atc_C = atc_B, atc_D = atc_C, atc_E = atc_D, count_BCDE = atleast_exposure_count), 
                   by = c("atc_B", "atc_C", "atc_D", "atc_E")) %>%
  dplyr::mutate(overrep_from_1 = round(atleast_exposure_count / (count_A / n_windows * count_B / n_windows * count_C / n_windows * count_D / n_windows * count_E), 3),
                overrep_from_ABCD_E = round(atleast_exposure_count / (count_ABCD / n_windows * count_E), 3),
                overrep_from_ABCE_D = round(atleast_exposure_count / (count_ABCE / n_windows * count_D), 3),
                overrep_from_ABDE_C = round(atleast_exposure_count / (count_ABDE / n_windows * count_C), 3),
                overrep_from_ACDE_B = round(atleast_exposure_count / (count_ACDE / n_windows * count_B), 3),
                overrep_from_BCDE_A = round(atleast_exposure_count / (count_BCDE / n_windows * count_A), 3),
                overrep_min = pmin(overrep_from_ABCD_E, overrep_from_ABCE_D, overrep_from_ABDE_C, overrep_from_ACDE_B, overrep_from_BCDE_A))

atleast_2s = atleast_2s %>% dplyr::select(atc_A, atc_B, atleast_exposure_count, overrep_1 = overrep, overrep_N = overrep)
atleast_3s = atleast_3s %>% dplyr::select(atc_A, atc_B, atc_C, atleast_exposure_count, overrep_1 = overrep_from_1, overrep_N = overrep_min)
atleast_4s = atleast_4s %>% dplyr::select(atc_A, atc_B, atc_C, atc_D, atleast_exposure_count, overrep_1 = overrep_from_1, overrep_N = overrep_min)
atleast_5s = atleast_5s %>% dplyr::select(atc_A, atc_B, atc_C, atc_D, atc_E, atleast_exposure_count, overrep_1 = overrep_from_1, overrep_N = overrep_min)

atleast_2s = atleast_2s %>% dplyr::rename(observe_per_expect_N1 = overrep_N)
atleast_3s = atleast_3s %>% dplyr::rename(observe_per_expect_1s = overrep_1, observe_per_expect_N1 = overrep_N)
atleast_4s = atleast_4s %>% dplyr::rename(observe_per_expect_1s = overrep_1, observe_per_expect_N1 = overrep_N)
atleast_5s = atleast_5s %>% dplyr::rename(observe_per_expect_1s = overrep_1, observe_per_expect_N1 = overrep_N)

atc_classes_2s = atc_classes_2s %>%
  dplyr::left_join(atleast_2s, by = c("atc_A", "atc_B"))
atc_classes_3s = atc_classes_3s %>%
  dplyr::left_join(atleast_3s, by = c("atc_A", "atc_B", "atc_C"))
atc_classes_4s = atc_classes_4s %>%
  dplyr::left_join(atleast_4s, by = c("atc_A", "atc_B", "atc_C", "atc_D"))
atc_classes_5s = atc_classes_5s %>%
  dplyr::left_join(atleast_5s, by = c("atc_A", "atc_B", "atc_C", "atc_D", "atc_E"))

atc_classes_2s = atc_classes_2s %>% dplyr::select(-exposure_count)
atc_classes_3s = atc_classes_3s %>% dplyr::select(-exposure_count)
atc_classes_4s = atc_classes_4s %>% dplyr::select(-exposure_count)
atc_classes_5s = atc_classes_5s %>% dplyr::select(-exposure_count)

# Assemble: Hide counts <100: -----------------------------------------------------
atc_classes_2s = atc_classes_2s %>%
  dplyr::mutate(observe_per_expect_N1 = ifelse(atleast_exposure_count < 100, " ", observe_per_expect_N1),
                atleast_exposure_count = ifelse(atleast_exposure_count < 100, "<100", atleast_exposure_count),
                fraction_exact = ifelse(is.na(exact_exposure_count) | exact_exposure_count < 100, " ", fraction_exact),
                fraction_allwindows = ifelse(atleast_exposure_count < 100, " ", fraction_allwindows),
                exact_exposure_count = ifelse(is.na(exact_exposure_count), 0, exact_exposure_count),
                exact_exposure_count = ifelse(exact_exposure_count < 100, "<100", exact_exposure_count))

atc_classes_3s = atc_classes_3s %>%
  dplyr::mutate(observe_per_expect_N1 = ifelse(atleast_exposure_count < 100, " ", observe_per_expect_N1),
                observe_per_expect_1s = ifelse(atleast_exposure_count < 100, " ", observe_per_expect_1s),
                atleast_exposure_count = ifelse(atleast_exposure_count < 100, "<100", atleast_exposure_count),
                fraction_exact = ifelse(is.na(exact_exposure_count) | exact_exposure_count < 100, " ", fraction_exact),
                fraction_allwindows = ifelse(atleast_exposure_count < 100, " ", fraction_allwindows),
                exact_exposure_count = ifelse(is.na(exact_exposure_count), 0, exact_exposure_count),
                exact_exposure_count = ifelse(exact_exposure_count < 100, "<100", exact_exposure_count))

atc_classes_4s = atc_classes_4s %>%
  dplyr::mutate(observe_per_expect_N1 = ifelse(atleast_exposure_count < 100, " ", observe_per_expect_N1),
                observe_per_expect_1s = ifelse(atleast_exposure_count < 100, " ", observe_per_expect_1s),
                atleast_exposure_count = ifelse(atleast_exposure_count < 100, "<100", atleast_exposure_count),
                fraction_exact = ifelse(is.na(exact_exposure_count) | exact_exposure_count < 100, " ", fraction_exact),
                fraction_allwindows = ifelse(atleast_exposure_count < 100, " ", fraction_allwindows),
                exact_exposure_count = ifelse(is.na(exact_exposure_count), 0, exact_exposure_count),
                exact_exposure_count = ifelse(exact_exposure_count < 100, "<100", exact_exposure_count))

atc_classes_5s = atc_classes_5s %>%
  dplyr::mutate(observe_per_expect_N1 = ifelse(atleast_exposure_count < 100, " ", observe_per_expect_N1),
                observe_per_expect_1s = ifelse(atleast_exposure_count < 100, " ", observe_per_expect_1s),
                atleast_exposure_count = ifelse(atleast_exposure_count < 100, "<100", atleast_exposure_count),
                fraction_exact = ifelse(is.na(exact_exposure_count) | exact_exposure_count < 100, " ", fraction_exact),
                fraction_allwindows = ifelse(atleast_exposure_count < 100, " ", fraction_allwindows),
                exact_exposure_count = ifelse(is.na(exact_exposure_count), 0, exact_exposure_count),
                exact_exposure_count = ifelse(exact_exposure_count < 100, "<100", exact_exposure_count))

# Assemble: Save: -----------------------------------------------------
save(atc_classes_1s, atc_classes_2s, atc_classes_3s, atc_classes_4s, atc_classes_5s, file="db_final/db_drug_classes_12345.RData")

write.table(atc_classes_1s, file="db_final/db_atc_classes_1s.tsv", quote=FALSE,  row.name = FALSE,sep="\t")
write.table(atc_classes_2s, file="db_final/db_atc_classes_2s.tsv", quote=FALSE,  row.name = FALSE,sep="\t")
write.table(atc_classes_3s, file="db_final/db_atc_classes_3s.tsv", quote=FALSE,  row.name = FALSE,sep="\t")
write.table(atc_classes_4s, file="db_final/db_atc_classes_4s.tsv", quote=FALSE,  row.name = FALSE,sep="\t")
write.table(atc_classes_5s, file="db_final/db_atc_classes_5s.tsv", quote=FALSE, row.name = FALSE, sep="\t")



