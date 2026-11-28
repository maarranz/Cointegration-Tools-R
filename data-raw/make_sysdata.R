# Assembles R/sysdata.rda from the fitted surfaces. Run from data-raw/ after
# eg_response_surfaces.R and, when available, ecm_response_surfaces.R and
# po_response_surfaces.R.
eg_rs <- readRDS("eg_rs.rds")
objs <- "eg_rs"
if (file.exists("ecm_rs.rds")) {
  ecm_rs <- readRDS("ecm_rs.rds")
  objs <- c(objs, "ecm_rs")
}
if (file.exists("po_rs.rds")) {
  po_rs <- readRDS("po_rs.rds")
  objs <- c(objs, "po_rs")
}
save(list = objs, file = "../R/sysdata.rda", compress = "xz")
cat("R/sysdata.rda written:", paste(objs, collapse = ", "), "\n")
