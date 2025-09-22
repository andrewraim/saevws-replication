SAIPE dataset (2023) put together on 3/3/2025 using `build-data.R`.

Data and sources are:

- ACS 5-year Table S1701
  Tidycensus

- SNAP 2022 - following SAIPE guidance, using one year previous (`cntysnap.csv`)
  <https://www.census.gov/data/datasets/time-series/demo/saipe/model-tables.html>
  
- 2023 Population Estimates from PEP (`co-est2023-alldata.csv`)
  <https://www.census.gov/data/datasets/time-series/demo/popest/2020s-counties-total.html>

- Sampled Housing Units in ACS
  Tidycensus Table B98001
   
- ACS response rate - precent of households that provided a response
  Tidycensus Table B98021

