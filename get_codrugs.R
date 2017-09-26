# Quinn & Shah, A dataset quantifying polypharmacy in the United States (2017)
# Script to identify drugs used concomitantly with a given drug of interest
# KQuinn 2017

library(readr)

drugs_1s = readr::read_tsv("https://datadryad.org/bitstream/handle/10255/dryad.158543/db_drugs_1s.tsv")
drugs_2s = readr::read_tsv("https://datadryad.org/bitstream/handle/10255/dryad.158544/db_drugs_2s.tsv")

# Define drug_of_interest: (For example: metformin, oxycodone)
drug_of_interest = "metformin"

drugs_1s_of_interest = drugs_1s %>%
  dplyr::filter(drug_name_A == drug_of_interest)

total_exposure_count = drugs_1s_of_interest$atleast_exposure_count

drugs_2s_of_interest = drugs_2s %>%
  dplyr::filter(drug_name_A == drug_of_interest | drug_name_B == drug_of_interest) %>%
  dplyr::mutate(fraction_all_drug_of_interest = round(atleast_exposure_count / total_exposure_count, 4))

print(head(drugs_2s_of_interest, 10))
print(head(drugs_2s_of_interest %>% arrange(-observe_per_expect_N1), 10))
