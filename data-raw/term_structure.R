# The term-structure data of the cointegration chapter's applied example:
# monthly 3-month Treasury bill rate (TB3MS) and 10-year Treasury constant
# maturity rate (GS10) from FRED. Run from data-raw/ on a machine with a FRED
# API key in ~/.Renviron (FRED_API_KEY=...):  Rscript term_structure.R
# Writes ../data/term_structure.rda; the download date is stored in the
# "retrieved" attribute and should be copied into R/data.R.
library(fredr)
fredr_set_key(Sys.getenv("FRED_API_KEY"))
get <- function(id) {
  d <- fredr(series_id = id)
  d <- d[!is.na(d$value), c("date", "value")]
  names(d)[2] <- tolower(id)
  d
}
term_structure <- merge(get("TB3MS"), get("GS10"), by = "date")
term_structure$date <- as.Date(term_structure$date)
term_structure <- term_structure[order(term_structure$date), ]
rownames(term_structure) <- NULL
attr(term_structure, "retrieved") <- Sys.Date()
dir.create("../data", showWarnings = FALSE)
save(term_structure, file = "../data/term_structure.rda", compress = "xz")
cat(sprintf("term_structure: %d months, %s to %s (retrieved %s)\n", nrow(term_structure),
            format(min(term_structure$date)), format(max(term_structure$date)), Sys.Date()))
