#
# Dynamic treatment rules for applied analysis.
#

g_EBRT <- function(data)
{
  V <- data[["TIME"]]

  A <- data[["LAG_ADJUVANT_TREATMENT"]]
  S <- data[["TREATMENT_READY"]]

  idx <- data[["INDEX"]]

  a_dtr <- ifelse(
    idx == DEFAULT_START_INDEX,
    "1= None",
    ifelse(S == 1 & V < 0.5, "2= EBRT", "1= None")
  )

  if (!is.null(A))
  {
    a_dtr <- ifelse(
      A == "1= None",
      a_dtr,
      as.character(A)
    )
  }

  return(factor(a_dtr, levels = TREATMENT_OPTIONS))
}

g_VBT <- function(data)
{
  V <- data[["TIME"]]

  A <- data[["LAG_ADJUVANT_TREATMENT"]]
  S <- data[["TREATMENT_READY"]]

  idx <- data[["INDEX"]]

  a_dtr <- ifelse(
    idx == DEFAULT_START_INDEX,
    "1= None",
    ifelse(S == 1 & V < 0.5, "3= VBT", "1= None")
  )

  if (!is.null(A))
  {
    a_dtr <- ifelse(
      A == "1= None",
      a_dtr,
      as.character(A)
    )
  }

  return(factor(a_dtr, levels = TREATMENT_OPTIONS))
}

g_Both <- function(data)
{
  V <- data[["TIME"]]

  A <- data[["LAG_ADJUVANT_TREATMENT"]]
  S <- data[["TREATMENT_READY"]]

  idx <- data[["INDEX"]]

  a_dtr <- ifelse(
    idx == DEFAULT_START_INDEX,
    "1= None",
    ifelse(S == 1 & V < 0.5, "4= Both", "1= None")
  )

  if (!is.null(A))
  {
    a_dtr <- ifelse(
      A == "1= None",
      a_dtr,
      as.character(A)
    )
  }

  return(factor(a_dtr, levels = TREATMENT_OPTIONS))
}

g_None <- function(data)
{
  return(factor("1= None", levels = TREATMENT_OPTIONS))
}
