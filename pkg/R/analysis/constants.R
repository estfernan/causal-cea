#
# Constant functions and variables used within project functions.
#

# constant functions: do not modify
MST <- function(T_, Y, ...) mean(T_, na.rm = TRUE)
MC <- function(T_, Y, ...) mean(Y, na.rm = TRUE)
SM <- function(T_, Y, ...) c(MST = MST(T_, Y, ...), MC = MC(T_, Y, ...))
MV <- function(T_, Y, kappa = 1, ...) kappa * mean(T_, na.rm = TRUE) - mean(Y, na.rm = TRUE)

# constant placeholders: do not modify
EMPTY_STRING <- ""
FALSE_BOOLEAN_INTEGER <- 0
NA_VAR <- NA
NULL_VAR <- NULL
PLACEHOLDER_CLAIM_ID <- "XXXXXXXXXXXXXXX"
PLACEHOLDER_CLAIM_SOURCE <- "CREATED"
TRUE_BOOLEAN_INTEGER <- 1

# constant variables: do not modify
CENSOR_EVENT_NAME <- "Censoring"
CENSOR_EVENT_TYPE <- -1
DEATH_EVENT_NAME <- "Death"
DEATH_EVENT_TYPE <- 0
DEFAULT_START_INDEX <- 1
DEFAULT_START_TIME <- 0
NO_TREATMENT <- "None"
P_BASELINE <- 3
P_TIMEVARYING <- 2
REGULAR_EVENT_NAME <- "Regular"
REGULAR_EVENT_TYPE <- 1
STAN_EXTENSION <- ".stan"
STAN_PATH <- "pkg/stan"

# column names: do not modify
BASELINE_COVARNAMES <- c("X1", "X2", "X3")
DEFAULT_COLUMNS <- c("TIME", "GAPTIME", "STATUS", "COST", "BALANCE", "TREATMENT")
KEYS <- c("ID", "INDEX")
LAG_VARS <- c("COST", "TREATMENT", "L1", "L2")
OUTCOME_COLUMNS <- c("TIME", "BALANCE")
TIMEVARYING_COVARNAMES <- c("L1", "L2")

# factor levels: do not modify
CENSORING_REASONS <- c("12/31/2019", "last Medicare claim")
COMORBIDITY_LEVELS <- c("None", "Low", "High")
EVENT_TYPES <- c("Regular", "Death", "Censoring")
FIGO_STAGES <- c("IA", "IB")
HISTOL_GRADES <- c("1= Grade 1", "2= Grade 2", "3= Grade 3")
INCOME_LEVELS <- c("< $40K", "$40K to < $55K", "$55K to < $75K", "$75K +")
RACE_CATEGORIES <- c("1= White", "2= Black", "3= Other")
TREATMENT_OPTIONS <- c("1= None", "2= EBRT", "3= VBT", "4= Both")

# Stan model names: do not modify
BERNOULLI_GLM <- "bernoulli_glmreg"
BERNOULLI_SIMPLE <- "bernoulli"
CATEGORICAL_GLM <- "categorical_glmreg"
CATEGORICAL_SIMPLE <- "categorical"
CONTINUATION_RATIO <- "continuation_ratio_dthreg"
EXPONENTIAL_PHREG <- "exponential_phreg"
GAMMA_GOMPERTZ_PMREG <- "gamma_gompertz_pmreg"
GAMMA_WEIBULL_PMREG <- "gamma_weibull_pmreg"
GAMMA_PIECEWISE_PMREG <- "gamma_piecewise_pmreg"
GAUSSIAN_GLM <- "gaussian_glmreg"
GAUSSIAN_SIMPLE <- "gaussian"
PIECEWISE_EXPONENTIAL_PHREG <- "piecewise_exponential_phreg"
PIECEWISE_EXPONENTIAL_GAR1_PHREG <- "piecewise_exponential_gAR1_phreg"
WEIBULL_PHREG <- "weibull_econ_phreg"
