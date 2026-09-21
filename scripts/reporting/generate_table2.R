#!/usr/bin/env Rscript

library(kableExtra)
library(tableone)
library(dplyr)
library(tibble)

load("data/obsdat.RData")

DATA_COLUMNS <- c("AGE", "WHITE", "MARRIED", "HIGH_SES", "METROPOLITAN", "NODES_EXAMINED", "FIGO_STAGE_IA", "DEATH", "ELAPSED_TIME", "TOTAL_COST", "CLAIM_COUNT", "AVERAGE_GAPTIME", "HIGH_COMORBIDITY", "ANY_OTHER_TX", "TREATMENT_INITIATED", "AVERAGE_COST")
TREATMENT_OPTIONS <- c("1= None", "2= EBRT", "3= VBT")
END_REASONS <- c("Death", "last Medicare claim", "12/31/2019")

patient_data <- obsdat %>%
  filter(!is.na(GAPTIME)) %>%
  arrange(ID, INDEX) %>%
  summarise(
    AGE = first(AGE),
    WHITE = first(RACE_WBO == 1),
    MARRIED = first(as.logical(MARRIED)),
    HIGH_SES = first(as.numeric(MEDIAN_INCOME) > 2),
    METROPOLITAN = first(as.logical(METROPOLITAN)),
    NODES_EXAMINED = first(as.logical(NODES_EXAMINED)),
    FIGO_STAGE_IA = first(FIGO_STAGE_IB == 0),
    DEATH = sum(DEATH, na.rm = TRUE) > 0,
    ELAPSED_TIME = last(TIME),
    TOTAL_COST = sum(UPCOMING_COST),
    CLAIM_COUNT = n(),
    AVERAGE_GAPTIME = mean(GAPTIME),
    HIGH_COMORBIDITY = any(as.logical(HIGH_COMORBIDITY)),
    ANY_OTHER_TX = any(as.logical(ANY_OTHER_TX)),
    TREATMENT_INITIATED = any(as.logical(TREATMENT_READY)),
    ADJUVANT_TREATMENT = factor(last(ADJUVANT_TREATMENT, order_by = INDEX), levels = TREATMENT_OPTIONS),
    AVERAGE_COST = mean(UPCOMING_COST),
    .by = ID
  )

tab1 <- CreateTableOne(
  vars = DATA_COLUMNS,
  strata = "ADJUVANT_TREATMENT",
  data = patient_data,
  test = FALSE,
  smd = TRUE,
  addOverall = TRUE
)

tab1Mat <- print(tab1, smd = TRUE, quote = FALSE, noSpaces = TRUE, printToggle = FALSE)

tab1Mat %>%
  as.data.frame() %>%
  rownames_to_column() %>%
  kbl(format = "latex", booktabs = TRUE, escape = FALSE, align = "lccccc", linesep = "") %>%
  save_kable("results/table2.tex")
