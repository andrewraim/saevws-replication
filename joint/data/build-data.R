## Code to build a 2023 SAIPE dataset

library(dplyr)
library(tidyverse)
library(tidycensus)
library(data.table)
library(readxl)

# -------------------------------------------------------- #
# ----------------- Pull ACS 5-year 2023 ----------------- #
# -------------------------------------------------------- #
acs5 = tidycensus::load_variables(2023, "acs5")
# view(acs5)

acs5yr = get_acs(geography = "county", table="S1701", year = 2023, survey="acs5") %>%
	filter(variable=="S1701_C02_001") %>%
	mutate(pov_se=moe/1.645) %>% # Convert from 90% moe to se
	rename(pov_count = estimate)

# -------------------------------------------------------- #
# ----------------- Pull ACS 1-year 2023 ----------------- #
# -------------------------------------------------------- #
acs1 = tidycensus::load_variables(2023, "acs1")
# view(acs5)

## This data is severely reduced, probably many counties are privacy-suppressed
acs1yr = get_acs(geography = "county", table="S1701", year = 2023, survey="acs1") %>%
	filter(variable=="S1701_C02_001") %>%
	mutate(pov_se=moe/1.645) %>% # Convert from 90% moe to se
	rename(pov_count = estimate)

# -------------------------------------------------------- #
# ----------------- Load SNAP benefit data --------------- #
# -------------------------------------------------------- #

snap2022 = fread('cntysnap.csv',
		skip=2,
		colClasses = c('State' = 'character', 'County' = 'character')) %>%
	mutate(GEOID = paste0(State, County)) %>%
	rename(snap = 'July 2022') %>%
	dplyr::select(GEOID, snap)

dat_2023 = acs5yr %>%
	inner_join(snap2022, by=join_by(GEOID)) %>%
	dplyr::select(-c(variable, moe))

# -------------------------------------------------------- #
# ----------------- Load Housing Unit Sample ------------- #
# -------------------------------------------------------- #
hu2023 = get_acs(geography = "county", table="B98001", year = 2023, survey="acs5") %>%
	filter(variable=="B98001_002") %>%
	rename(hu_sampled = estimate) %>%
	select(-NAME)

# This gives population samples, not sample units, so will skip for now
# totalsamp2021 = get_acs(geography = "county", table="B98003", year = 2021, survey="acs5") %>%
# 	filter(variable=="B98003_001") %>%
# 	rename(total_sampled = estimate)

dat_2023 = dat_2023 %>%
	inner_join(hu2023, by=join_by(GEOID)) %>%
	dplyr::select(-variable)

# -------------------------------------------------------- #
# ---------------------- Load PEP data ------------------- #
# -------------------------------------------------------- #

pep2023 = fread('co-est2023-alldata.csv',
	colClasses = c('SUMLEV' = 'character', 'REGION' = 'character',
		'DIVISION' = "character", 'STATE' = 'character',
		'COUNTY' = "character")) %>%
	mutate(GEOID = paste0(STATE, COUNTY)) %>%
	rename(pep = 'POPESTIMATE2023') %>%
	dplyr::select(GEOID, pep)

dat_2023 = dat_2023 %>%
	inner_join(pep2023, by=join_by(GEOID))


# -------------------------------------------------------- #
# ---------------- Pull Response Rate Data --------------- #
# -------------------------------------------------------- #
acs5 = tidycensus::load_variables(2023, "acs5")
view(acs5)

acs5_rr = get_acs(geography = "county", table="B98021", year = 2023, survey="acs5") %>%
	filter(variable=="B98021_001") %>%
	dplyr::select(GEOID, estimate) %>%
	rename(rrate = estimate)

dat_2023 = dat_2023 %>%
	inner_join(acs5_rr, by=join_by(GEOID))

write_csv(dat_2023, file = "saipe.csv")

