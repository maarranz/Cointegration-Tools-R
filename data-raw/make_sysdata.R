# Assembles R/sysdata.rda from the fitted surfaces. Run from data-raw/ after
# eg_response_surfaces.R.
eg_rs <- readRDS("eg_rs.rds")
save(eg_rs, file = "../R/sysdata.rda", compress = "xz")
cat("R/sysdata.rda written (eg_rs for nvars =", paste(names(eg_rs), collapse = ", "), ")\n")
