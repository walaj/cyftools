load_redcap <- function() {
  
  # Read the data
  redcap <- fread("/Users/jeremiahwala/Sorger/projects/orion/clinical/74OrionCasesDeidentifiedMasterList2023_scOPupdate_better.csv")
  
  # Clean column names
  setnames(redcap, gsub("\\s+", "", colnames(redcap)))
  
  # Convert columns to lowercase
  redcap[, c("TIL", "MMRIHC") := lapply(.SD, tolower), .SDcols = c("TIL", "MMRIHC")]
  
  # Convert columns to numeric or integer
  convert_to_type <- function(col, type) {
    redcap[, (col) := sapply(get(col), function(x) {
      if(is.na(x) || x == "NA") {
        return(NA)
      } else {
        return(type(x))
      }
    })]
  }
  
  # convert some columns to numeric
  for (col in c("Renato_TMB", "OP_TMB") ) convert_to_type(col, as.numeric)
  for (col in c("Hypermutant","Recurrence", "PFSDays", "PFSCensor","OSDays","OSCensor",
             "BRAF","p53","NRAS","PIK3CA","APC","FBXW7","PTEN","TCF7L2","SOX9","CTNNB1")) 
    convert_to_type(col, as.integer)
 
  # Additional processing
  redcap[, tmb := pmax(Renato_TMB, OP_TMB, na.rm = TRUE)]
  redcap[, c("Renato_TMB", "OP_TMB") := NULL]
  redcap[, stage_num := roman_to_arabic(Stage_At_Diagnosis)]
  redcap[, pMMR := as.integer(MMRIHC %in% c("intact","Intact"))]
  redcap[, MMRIHC := fifelse(tolower(MMRIHC) == "intact", "proficient", 
                             fifelse(grepl("absent", tolower(MMRIHC)), "deficient", MMRIHC))]
  
  assign("redcap", redcap, envir = .GlobalEnv)
  invisible(NULL)
}
