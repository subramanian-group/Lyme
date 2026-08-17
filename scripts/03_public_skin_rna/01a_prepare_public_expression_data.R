#Title: Clean Differential Expression Data
#Author: Clifford Rostomily

# Options -----------------------------------------------------------------

options(stringsAsFactors = FALSE)

# Paths and cache policy --------------------------------------------------

if (!requireNamespace("here", quietly = TRUE)) {
  stop("Package 'here' is required to resolve paths from the repository root.")
}

skin_raw_dir <- here::here("data", "raw", "03_public_skin_rna")
skin_processed_dir <- here::here("data", "processed", "03_public_skin_rna")
public_data_input <- file.path(skin_raw_dir, "Public Data.RData")
public_data_output <- file.path(skin_processed_dir, "publicDataClean.RData")
force_rebuild <- tolower(Sys.getenv("LYME_REBUILD_PREPROCESSING", "false")) %in%
  c("1", "true", "yes")

dir.create(skin_processed_dir, recursive = TRUE, showWarnings = FALSE)

if (file.exists(public_data_output) && !force_rebuild) {
  load(public_data_output)
  message("Using cached processed object: ", public_data_output)
} else {
  if (!file.exists(public_data_input)) {
    stop(
      "Cannot rebuild publicDataClean.RData because the source object is missing: ",
      public_data_input,
      "\nUse the deposited processed cache or place Public Data.RData at this path."
    )
  }

  required_packages <- c("GEOquery", "edgeR", "limma")
  missing_packages <- required_packages[
    !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
  ]
  if (length(missing_packages)) {
    stop("Missing required package(s): ", paste(missing_packages, collapse = ", "))
  }

  invisible(lapply(required_packages, library, character.only = TRUE))

# Functions ---------------------------------------------------------------

# No project helper functions are required by this script.

# Load Data ---------------------------------------------------------------

load(public_data_input)

data2 = data


# GSE77929 ----------------------------------------------------------------

dat = data2[["GSE77929"]]$GSE77929

factors = as.data.frame(pData(dat)[,c("geo_accession", "individual:ch1", "time:ch1")])
factors$type = rep("patient", nrow(factors))
colnames(factors) = c("gsm", "id", "time", "type")
factors$time[grep("V1", factors$time)] = 0
factors$time[grep("V2", factors$time)] = 21
factors$time[grep("V5", factors$time)] = 182
factors$id = gsub(" |-", "_", factors$id)
factors$id = gsub("patient_", "", factors$id)
factors = factors[, c("gsm", "id", "type", "time")]
factors$time = as.numeric(factors$time)

dat2 = exprs(dat)
meta = factors
features = rownames(dat2)

dat = list(data = dat2, meta = meta, features = features)
data2[["GSE77929"]] = dat

# GSE63085 ----------------------------------------------------------------

dat = data[["GSE63085"]]$GSE63085

factors = as.data.frame(pData(dat)[,c("geo_accession", "individual:ch1", "time:ch1",grep("characteristics",colnames(pData(dat)),value = T))])
factors$type = gsub("(.*) .*", "\\1", factors$`individual:ch1`)
factors = factors[,c(1,2,3,8)]
colnames(factors) = c("gsm", "id", "time", "type")
factors$time[grep("V1", factors$time)] = 0
factors$time[grep("V2", factors$time)] = 21
factors$time[grep("V5", factors$time)] = 182
factors$id = gsub(".* ", "", factors$id)
factors$id = gsub(" |-", "_", factors$id)
factors = factors[, c("gsm", "id", "type", "time")]
factors$time = gsub("control", 0, factors$time)
factors$time = as.numeric(factors$time)
factors$type_time = paste0(factors$type,"_",factors$time)

dat2 = dat@assayData$exprs
meta = factors
features = fData(dat)


## Normalize data
d0 <- DGEList(dat2)
d0 <- calcNormFactors(d0)
cutoff <- 2
drop <- which(apply(cpm(d0), 1, max) < cutoff)
d <- d0[-drop,]
mm <- model.matrix(~0 + meta$type_time)
y <- voom(d, mm, plot = T)
# tmp <- voom(d0, mm, plot = T)

dat = list(data = y, meta = meta, features = features)
data2[["GSE63085"]] = dat

# GSE68741 ----------------------------------------------------------------
dat = data[["GSE68741"]]$GSE68741

groups = c("none", "Borrelia_afzelii", "Borrelia_garinii", "Borrelia_burgdorferi_sensu_stricto")
factors = as.data.frame(pData(dat)[,c("geo_accession", "infection:ch1")])
factors$id = paste0("GSE68741", 1:nrow(factors))
factors$time = 1
colnames(factors) = c("gsm", "type", "id", "time")
factors$type = gsub(" ", "_", factors$type)
factors$type = gsub("Borrelia_afzeliiburgdorferi", "Borrelia_burgdorferi", factors$type)
factors$type = gsub("none", "control", factors$type)
factors = factors[, c("gsm", "id", "type", "time")]
factors$time = as.numeric(factors$time)

meta = factors
dat2 = dat@assayData$exprs
features = fData(dat)
rownames(dat2) = features$`Gene Symbol`

dat = list(data = dat2, meta = meta, features = features)
data2[["GSE68741"]] = dat

# GSE84479 ----------------------------------------------------------------

dat = data[["GSE84479"]]$GSE84479

factors = as.data.frame(pData(dat)[,c("geo_accession", 
                                      "subject id:ch1", 
                                      "visit number (1 = baseline):ch1", 
                                      "# of symptoms:ch1", 
                                      "multiple em (y = multiple erythema migrans lesions, n = single erythema migrans lesion):ch1", 
                                      "illness duration (days of illness prior to treatment beginning):ch1", 
                                      "ptlds status [ptlds status  0 = returned to health, 1 = symptoms only, 2 = ptlds]:ch1",
                                      "Sex:ch1",
                                      "mcgill pain total score (range 0-45):ch1",
                                      "fatigue severity [range from 9-63]  scores ≥36 are considered \"\"high\"\"]:ch1",
                                      "elevated liver enzyme [ast, alk, or alk-phos at pre-treatment visit (v1). 0 = no, 1 = yes]:ch1")])
colnames(factors) = c("gsm", "id", "time", "n_symptons", "multiple_em", "illness_duration", 
                      "status", "sex", "mcgill pain total score (range 0-45)",
                      "fatigue severity [range from 9-63]  scores ≥36 are considered \"high\"]:ch1",
                      "elevated liver enzyme [ast, alk, or alk-phos at pre-treatment visit (v1). 0 = no, 1 = yes]")

for(i in 1:6){
  factors$time = gsub(as.character(i), c("0", "21", "28", "84", "182", "365")[i], factors$time)
}
factors$time[factors$time == "7"] = 999

factors$type = rep("patient", nrow(factors))
controls = factors$id[factors$time == "0"][factors$n_symptons[factors$time == "0"] == "N/A"]
factors$type[factors$id %in% controls] = "control"
factors$time = as.numeric(factors$time)
factors$type_time = paste0(factors$type,"_",factors$time)

meta = factors
dat2 = dat@assayData$exprs
features = fData(dat)
rownames(dat2) = fData(dat)$SPOT_ID

dat = list(data = dat2, meta = meta, features = features)
data2[["GSE84479"]] = dat

# GSE55815 ----------------------------------------------------------------

dat = data[["GSE55815"]]$GSE55815

factors = as.data.frame(pData(dat)[,c("geo_accession", "subject id:ch1", "group:ch1", "visit number (1 = baseline):ch1", 
                                      "acute 2 tier status:ch1", "cd3+, cxcr3+, cd4+:ch1","cd3+, cxcr3+, cd8+:ch1", 
                                      "elevated liver enzyme:ch1", "illness duration (days):ch1", 
                                      "multiple em:ch1", "Sex:ch1", "symptom number:ch1")])
colnames(factors) = c("gsm", "id", "type", "time", 
                      "2 tier status", "cd3+cxcr3+cd4+", "cd3+cxcr3+cd8+",
                      "elevated liver enzyme", "illness duration",
                      "multiple_em", "sex", "n_symptom")
factors$type = gsub("Case", "patient", factors$type)
factors$type = gsub("Control", "control", factors$type)
factors$time = rep(0, nrow(factors))


meta = factors
dat2 = dat@assayData$exprs
features = fData(dat)
rownames(dat2) = features$SPOT_ID

dat = list(data = dat2, meta = meta, features = features)
data2[["GSE55815"]] = dat


# GSE154916 ----------------------------------------------------------------

dat1 = data[["GSE154916"]]$`GSE154916-GPL570_series_matrix.txt.gz`
dat2 = data[["GSE154916"]]$`GSE154916-GPL571_series_matrix.txt.gz`
dat = cbind(dat1@assayData$exprs,dat2@assayData$exprs)
meta = rbind(pData(dat1),pData(dat2))
features = rbind(fData(dat1),fData(dat2))
rownames(dat) = features$`Gene Symbol`[match(rownames(dat),features$ID)]

colnames(meta) = gsub(".ch1","",colnames(meta))
colnames(meta)[6] = "sample_type"

dat = list(data = dat, meta = meta, features = features)
data2[["GSE154916"]] = dat

# Output ------------------------------------------------------------------

publicData = data2[!names(data2) %in% c("GSE6092", "GSE31740")]

save(publicData, file = public_data_output)
message("Wrote processed object: ", public_data_output)
}

print(sessionInfo())
