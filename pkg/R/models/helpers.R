# """
# File name: helpers.R
#
# Description:
#   Functions used to facilitate the model fitting procedures.
# """

split_surv <- function(data, time_cuts, event_code, formula_rhs) {
  # """
  # Helper to split survival data.
  # """
  fm <- sprintf("Surv(gaptime, event_type == %i) ~ %s", event_code, formula_rhs)

  outdat <- survival::survSplit(
    data = data,
    formula = as.formula(fm),
    cut = time_cuts,
    id = "patient_id",
    episode = "interval",
    end = "tstop",
    event = "event_type"
  )

  return(outdat)
}
