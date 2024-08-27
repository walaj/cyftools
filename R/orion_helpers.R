generate_random_string <- function(string_length = 8) {
  characters <- c(letters, LETTERS, 0:9)  # alphanumeric characters
  random_string <- paste0(sample(characters, string_length, replace = TRUE), collapse = "")
  return(random_string)
}

# logistic regression plot
lr_plot__ <- function(lrr, threshold = 1e3) {
  
  # plot the results, first need to format the results
  exclude.models <- lrr[abs(Value) > threshold | is.na(Value), i]
  lrr[!i %in% exclude.models, se := sd(Value), by = .(Marker, Region)]
  lrr[!i %in% exclude.models, mean := mean(Value), by=.(Marker, Region)]
  lrrp <- lrr[Marker != "(Intercept)" & !i %in% exclude.models, .(Value = mean(Value), se = mean(se)), by = .(Marker, Region)]

  # Define the plot
  g <- ggplot(lrrp, aes(x = Marker, y = Value, fill = factor(Region))) +
    geom_bar(stat = "identity", position = "dodge") +
    geom_errorbar(aes(ymin = Value - se, ymax = Value + se), width = 0.2,
                  position = position_dodge(0.9)) + 
    scale_fill_manual(name="Region", values=c("#66c2a5","#fc8d62","#8da0cb","#e78ac3"), labels=c("Stroma","Tumor","Tumor border","Whole slide")) + 
    ylab("Log-odds Hypermutant") + xlab("") + coord_flip() + theme_bw()
  g
}
run_cox_model <- function(rd, numeric_columns, cat_columns) {
  
  redcap.norm <- data.table::copy(rd) # Assuming 'data' is your data.table and 'cols_to_standardize' is the vector of column names
  
  # standardize numeric columns
  redcap.norm[, (numeric_columns) := lapply(.SD, function(x) scale(x)), .SDcols = numeric_columns]

  ## pre-flipped
  if (redcap.norm[Slide_ID=="LSP10353", PFSCensor] == 1) { 
    rd[, PFSC := PFSCensor]
    ## not-yet-flipped
  } else { 
    redcap.norm[, PFSC := ifelse(PFSCensor==0,1,0)]
  }
  
  # Assuming 'response_var' is the name of the response variable, and 'cols' is the vector of predictor names
  predictors <- paste(c(numeric_columns,cat_columns), collapse = " + ")
  formula_str <- paste("Surv(PFSDays, PFSC)", "~", predictors)
  lm_formula.path <- as.formula(formula_str)
  
  # run the model
  model.cox <- coxph(formula = lm_formula.path, data = redcap.norm)
  
  return (model.cox)
}

# load all the data
load_cell_types__ <- function() {
  set.seed(1337)
  fff <- lapply(seq(40) , function(x) {
    str <- formatC(x, width = 2, flag = "0")
    #if (str=="33")
    #  str = "33_01"
    cat(paste("Working on",str),"\n")
    f <- fread(cmd=paste0("gunzip -c ~/Sorger/orion/jerry_40_orion/dataC",str,".csv.gz"))
    f[, c("kmean16","CD3pR1","CD3pR2","CD8pR1","CD8pR2","R0","R1",
          "R2","CD31nb","CD31intPan_CK","CD31intCD45","patch_id",
          "prediction_orion","prediction_adjacent", "Y","X","Hoechst",
          "Eccentricity","Solidity","Extent","Orientation","Area",
          "MajorAxisLength","MinorAxisLength") := NULL]
    f[,c("CellID","AF1","CD31","CD45","Argo550","CD4","FOXP3","CD8a","CD45RO","CD20",
         "PD_L1","CD3e","CD163","E_cadherin","PD_1","Ki67","Pan_CK","SMA") := NULL]
    #f[, names(f) := lapply(.SD, as.integer)]
    #f[, (cols_to_factor) := lapply(.SD, as.factor), .SDcols = cols_to_factor]
    f[, sample := x]
    return(f)
    #return(f[sample(.N,max(2.5e5,.N),replace=FALSE)])
  })
  cells <- rbindlist(fff)
  setnames(cells, c("PD-1","PD-L1","Pan-CK","E-cadherin"), c("PD1","PDL1","PanCK","Ecadherin"))
  saveRDS(cells, "~/Sorger/cache/pd40_no_fluorescence.rds", compress=FALSE)
  #return(cells)
}

roman_to_arabic <- function(y) {
  sapply(y, function(x) {
    if (grepl("^IV",x))
      return(4)
    if (grepl("^III",x))
      return(3)
    if (grepl("^II", x))
      return(2)
    if (grepl("^I",x))
      return(1)
    return(NA)
  })
}

# load the excel clinical table and cleanup
load_redcap <- function() {
  
  #redcap <- as.data.table(readxl::read_excel("~/Sorger/orion/74 Orion Cases DE-IDENTIFIED MASTER LIST 2023.xlsx"))
  redcap <- fread("/Users/jeremiahwala/Dropbox/Sorger/projects/orion/clinical/74OrionCasesDeidentifiedMasterList2023_scOPupdate_better.csv")
  
  # remove whitespaces from names
  setnames(redcap, gsub("\\s+", "", colnames(redcap)))
  
  # set to lower
  redcap[, TIL := tolower(TIL)]
  redcap[, MMRIHC := tolower(MMRIHC)]
  
  # convert some columns to numeric
  cols_to_convert <- c("Renato_TMB", "OP_TMB")
  for (col in cols_to_convert) {
    redcap[, (col) := sapply(get(col), function(x) {
      if(is.na(x) || x == "NA") {
        return(NA)
      } else {
        return(as.numeric(x))
      }
    })]
  }
  
  # convert some columns to integer
  cols_to_convert <- c("Hypermutant","Recurrence", "PFSDays",
                       "PFSCensor","OSDays","OSCensor",
                       "BRAF","p53","NRAS","PIK3CA","APC",
                       "FBXW7","PTEN","TCF7L2","SOX9",
                       "CTNNB1")
  for (col in cols_to_convert) {
    redcap[, (col) := sapply(get(col), function(x) {
      if(is.na(x) || x == "NA") {
        return(NA)
      } else {
        return(as.integer(x))
      }
    })]
  }
  
  redcap[, tmb := pmax(Renato_TMB, OP_TMB, na.rm = TRUE)]
  redcap[, c("Renato_TMB", "OP_TMB") := NULL]
 
  # give stage a number
  redcap[, stage_num := roman_to_arabic(Stage_At_Diagnosis)]
  
  # p/dMMR
  redcap[, pMMR := as.integer(MMRIHC %in% c("intact","Intact"))]
  redcap[, MMRIHC := ifelse(tolower(MMRIHC)=="intact", "proficient",MMRIHC)]
  redcap[, MMRIHC := ifelse(grepl("absent",tolower(MMRIHC)), "deficient",MMRIHC)]
  
  assign("redcap", redcap, envir = .GlobalEnv)
  invisible(NULL)
}

contour_plot__ <- function(dt) {
  ggplot(data=dt, aes(x=Xt, y=Yt, color=variable)) + 
    theme_bw() + geom_density_2d(linewidth=1) +
    coord_fixed(ratio = 1)
}

#' Load or run code block and cache result
#'
#' If a cached result is available, the function reads the cached result from
#' the cache path. If a cached result is not available, the function executes
#' the code block and saves the result in the cache path.
#'
#' @param cache_path A character string specifying the file path to the cache.
#' @param code_block A function that contains the code block to execute if a
#'   cached result is not available.
#' @param output_name A character string specifying the name of the variable
#'   where the output of the code block should be assigned in the global
#'   environment.
#'
#' @return The output of the code block is assigned to a variable in the global
#'   environment with the name specified by \code{output_name}.
#'
#' @examples
#' # Create a code block to execute
#' my_code_block <- function() {
#'   # some code here
#' }
#'
#' # Load or run the code block
#' load_or_run(cache_path = "my_cache_path.rds",
#'             code_block = my_code_block,
#'             output_name = "my_output")
#'
#' @export
load_or_run <- function(cache_path, code_block, output_name) {
  if (file.exists(cache_path) && use.cache) {
    message("Loading cached result from ", cache_path)
    out <- readRDS(cache_path)
  } else {
    message(paste("Running code block to create", output_name))
    out <- code_block() #eval(parse(text = code_block))
    saveRDS(out, cache_path)
  }
  assign(output_name, out, envir = .GlobalEnv)
  invisible(NULL)
}

#' KM Plot
#'
#' Create a Kaplan-Meier plot to visualize survival data.
#'
#' @param rd A data.table containing the clinical data, with each row as a patient
#' @param parm The parameter (column name in rd) to plot survival curves stratified by.
#' @param parm.name The name of the parameter used in the plot legend.
#' @param parm.labels Labels for the parameter categories used in the plot legend.
#' @param pfs [TRUE] Logical to say whether to plot the PFS data (vs OS)
#'
#' @return A ggplot2 object displaying the Kaplan-Meier plot.
#'
#' @examples
#' km_plot(data, "cd3_hi", "CD3 Expression", c("Low", "High"))
#'
#' @export
km_plot <- function(rd, parm, parm.name, parm.labels,
                    flip=FALSE, flip.labels = FALSE, return.model=FALSE, 
                    return.data = FALSE, 
                    stage_adjust = FALSE) {
  
  library(survival)
  
  # some data checking
  if (!("PFSDays" %in% colnames(rd) || "OSDays" %in% colnames(rd))) {
    stop("Error: Data table must contain either PFSDays or OSDays column.")
  }
  
  if (!(parm %in% colnames(rd))) {
    stop("Error: Data table must contain the requested parameter")
  }
  
  #parm = "cd3_hi"
  #parm.name = "cd3_hi"
  #parm.labels = c("Low", "Hi")
  parm.vals = c("#af8dc3","#7fbf7b")
  if (flip.labels)
    parm.vals <- rev(parm.vals)
  
  # get the survival data
  #s <- rep(pfs, nrow(rd))
  
  if (!"PFSC" %in% colnames(rd)) {
    cat("\n***Need PFSC column***\n")
    return(0)
  } 
  
  surv_model <- Surv(rd$PFSDays, rd$PFSC)
  # }
  # 
  # ## pre-flipped
  # if (rd[Slide_ID=="LSP10353", PFSCensor] == 1) { 
  #   #cat("preflipped\n")
  #   rd[, PFSC := ifelse(rd$PFSCensor==0,0,1)]
  #   surv_model <- Surv(ifelse(s, rd$PFSDays, rd$OSDays), 
  #                      ifelse(s, rd$PFSCensor,ifelse(rd$OSCensor ==0,0,1)))
  # ## not-yet-flipped
  # } else { 
  #   #cat("needs flipped\n")
  #   rd[, PFSC := ifelse(rd$PFSCensor==0,1,0)]
  #   surv_model <- Surv(ifelse(s, rd$PFSDays, rd$OSDays), 
  #                      ifelse(s, ifelse(rd$PFSCensor==0,1,0),ifelse(rd$OSCensor ==0,1,0)))
  # }
  
  # store the model data in a data.table for ggplot ease
  # fit the survival model
  if (stage_adjust)
    model_fit = survfit(surv_model ~ get(parm) + stage_num, data = rd)
  else
    model_fit = survfit(surv_model ~ get(parm) , data = rd)
  
  
  # Compute the optimal cutpoint
  if (rd[, is.numeric(get(parm))]) {
    res.cut <- survminer::surv_cutpoint(rd, time = "PFSDays", event = "PFSC", variables = parm)
    rd[, cut := as.integer(get(parm) > as.numeric(res.cut$cutpoint)[1])]
    parm="cut"
    parm.labels=c("0"="Low","1"="High")
    cat("Cutpoint: ", as.numeric(res.cut$cutpoint)[1], "\n")
    if (stage_adjust)
      model_fit = survfit(surv_model ~ get(parm) + stage_num, data = rd)
    else
      model_fit = survfit(surv_model ~ get(parm), data = rd)
  }
  
  dt.plot <- data.table(time=model_fit$time, 
                        surv=model_fit$surv,
                        strata=as.vector(unlist(sapply(names(model_fit$strata), function(x) rep(x, model_fit$strata[x])))),
                        lower=model_fit$lower, upper=model_fit$upper, censor=model_fit$n.censor)
  dt.plot[, strata := gsub("get\\(parm\\)=", "", strata)]
  dt.1 <- dt.plot[!duplicated(strata)]; dt.1$time=1; dt.1$surv=1; dt.1$lower=1; dt.1$upper=1; dt.1$censor=0
  dt.plot <- rbind(dt.plot, dt.1)
  
  # # Fit a Cox proportional hazards model to the data
  if (stage_adjust)
    cox_model <- coxph(surv_model ~ get(parm) + stage_num, data = rd)
  else
    cox_model <- coxph(surv_model ~ get(parm) , data = rd)
  
  
  # get the HRF
  coefs <- coef(cox_model)
  
  # Invert the HR if "flip" is TRUE
  if(flip) {
    coefs <- -coefs
  }
  
  se <- sqrt(diag(vcov(cox_model)))
  hr <- list(hr=exp(coefs), lower=exp(coefs - 1.96 * se), 
             upper=exp(coefs + 1.96 * se))
  hr.lab <- paste0("HR: ", round(hr[[1]], 3), " [", round(hr[[2]], 2),"-", round(hr[[3]], 2), "]")
  hr.x <- max(dt.plot$time)*0.9
  hr.y <- max(dt.plot$surv)*0.9
 
  cat(hr.lab, "\n")
  
  if (return.model)
    return(hr)
  
  if (return.data)
    return(dt.plot)
  
  # plot the K-M curves
  g <- ggplot(dt.plot[time < 2400], aes(x = time/365, y = surv, color = strata)) +
    geom_ribbon(aes(ymin = lower,ymax = upper, fill=strata),alpha = 0.15, linetype=0, show.legend=FALSE) +
    geom_step(aes(color=strata), linewidth=1) +
    geom_point(data = dt.plot[censor == 1][time < 2400], aes(x=time/365,y=surv), shape="+", size=3, color="black") + 
    labs(x = "\nYears", y = ifelse(TRUE,"Progression-Free Survival\n (probability)","Overall Survival (probability)\n")) + theme_bw() + 
    ggtitle(hr.lab) + 
    scale_color_manual(name=parm.name,values=parm.vals,labels=parm.labels) + 
    scale_x_continuous(breaks=s<-seq(0,6), labels=s, expand=c(0,0)) + 
    scale_fill_manual(values=parm.vals, guide="none") + 
    #annotate("text", x = hr.x, y = hr.y,label=hr.lab,hjust = 1, vjust = 1) + 
    scale_y_continuous(limits=c(0,1), breaks= s<-seq(0, 1, by = 0.2), labels = s, expand=c(0,0)) +
    theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())

  return(g)
  
}

#' Frequency Histogram
#'
#' This function generates a frequency histogram given a data table and two parameter names.
#'
#' @param dt data table with samples and columns of numerical data to plot
#' @param parm.x string name of numerical parameter column to plot
#' @param parm.fill string name of categorical parameter column to use for color fill
#' @return A ggplot object displaying the frequency histogram
#' @keywords histogram
#' @examples
#' \dontrun{
#' freq_histogram(my_dt, "CD3ep", "Hypermutant")
#' }
freq_histogram <- function(dt, parm.x, parm.fill, y.label) {
  
  if (!all(c(parm.x, parm.fill, "sample") %in% colnames(dt))) 
    stop("Required column names are not in the data.table. ", 
         "Should be formated with columns 'sample' and then ", 
         "whatever values you want to plot (numerical - ", parm.x, 
         ") and color by (categorical - ", parm.fill, ")")
  
  s <- unique(dt[order(dt[, get(parm.x)])]$sample)
  cols <- c("steelblue", "orange", "blue","red")[seq_along(dt[, unique(get(parm.fill))])]
  g <- ggplot(dt, aes(x = factor(sample, levels=s), y = get(parm.x), fill = factor(get(parm.fill))))+
    geom_col(position = position_dodge(width = 0.8)) +
    scale_fill_manual(values = cols, name = parm.fill) +
    xlab("Sample") +
    ylab(y.label) +
    theme_minimal() + coord_flip()
  g
}

load_mutations__ <- function() {
  
  # load and format the data
  mpath <- "/Users/jeremiahwala/Sorger/orion/orion_crc_dropbox/data_mutations_extended.txt"
  dt <- fread(mpath)
  setnames(dt, "Tumor_Sample_Barcode","ID")
  dt[, ID := gsub("^(\\w{1})\\w{2}(\\d+)$", "\\1\\2", ID)]
  dt[, sample := as.numeric(gsub("C", "", ID))]
  
  # get the mutation count per patient
  dt[, mut_count := .N, by=sample]
  
  # some patients have multiple mutations in one gene, ignore for now
  dt<- unique(dt, by = c("sample", "Hugo_Symbol"))
  
  # add a "p" at end just to note that comes from mutations file
  dt[, Hugo_Symbol := paste0(Hugo_Symbol,"p")]

  setkey(dt, sample, Hugo_Symbol)
  # limit to only the 15 most prevalent genes
  tab <- sort(table(dt$Hugo_Symbol), decreasing=TRUE)
  tab15 <- tab[seq(15)]
  dt <- dt[Hugo_Symbol %in% names(tab15)]
  
  # recast so samples are rows
  dt_wide <- dcast(dt[, .(sample, Hugo_Symbol, mut_count)], mut_count + sample ~ Hugo_Symbol, fun.aggregate = length, fill = 0)
  
  return(dt_wide)
}

# get the confusion matrix between two columns in redcap
gene_sample_conf_mat__ <- function(a,x) {
  small <- redcap[,c(a, x), with=FALSE]
  small <- small %>% mutate(!!x := as.numeric(get(x) > 0)) ## !! unquotes x
  small <- small %>% mutate(!!a := as.numeric(get(a) > 0)) ## !! unquotes a
  
  xx = small[!is.na(get(a)) & !is.na(get(x)), get(a)]
  yy = small[!is.na(get(a)) & !is.na(get(x)), get(x)]
  conf_mat <- table(xx,yy)
  return(conf_mat)
}

# convienence function to load the data and save it
subset_pd__ <- function() {
  pd <- readRDS("~/Sorger/orion/jerry_40_orion/jerry_trim.rds")
  
  # create a simplified and melted data.table "rrm"
  rr <- pd[, .(Xt, Yt, E_cadherinp, CD3ep, Pan_CKp, PD_1p, PD_L1p, sample, CellID)]
  rr[, PanCK_PDL1 := as.integer(Pan_CKp > 0 & PD_L1p > 0)]
  rr[, PanCK_CD3ep := as.integer(CD3ep > 0 & PD_L1p > 0)]
  rrm <- melt(rr, id.vars=c("Xt","Yt","CellID","sample"))
  rrm <- rrm[value > 0, ]
  return(rrm)
}

redcap_to_json <- function() {
  
  # have a data.table to change names on
  dt <- data.table::copy(redcap)
  
  ## load radials
  #rad <- load_radials()
  rad[, c("crcid", "file_name") := NULL]
  dt <- merge(dt, rad, by = "ID", all.x = TRUE)
  names(dt) <- gsub("\\+", "plus", names(dt))
  
  dt[, Gender := ifelse(Gender=="M", "male", Gender)]
  dt[, Gender := ifelse(Gender=="F", "female", Gender)]
  
  dt <- dt[,.(LVI,PNI,TIL,Age,Gender, Histology,Grade,stage_num, Death,KRAS,p53,NRAS,BRAF,PIK3CA,
              APC,FBXW7,PTEN,TCF7L2,SOX9,CTNNB1,ID,MMRIHC,tmb,Recurrence,
              CD3_4_8_20r,PD1_20r,PanCK_Ecad_20r,PDL1_20r,CD68_20r,
              CD163_20r,CD68plusPDL1_20r,CD3plusPD1_20r,PanCKplusPDL1_20r,CD3_4_8_200r,
              PD1_200r,PanCK_Ecad_200r,PDL1_200r,CD68_200r,
              CD163_200r,CD68plusPDL1_200r,CD3plusPD1_200r,PanCKplusPDL1_200r)]
  
  dt[, MMRIHC := tolower(MMRIHC)]
  
  setnames(dt, c("Recurrence","LVI","PNI","TIL","Age","Gender","Histology","Death","KRAS",
             "p53","NRAS","BRAF","PIK3CA","APC","FBXW7","PTEN",
             "TCF7L2","SOX9","CTNNB1","stage_num","Grade","MMRIHC","tmb"),
           c("CLI_recurrence","CLI_lvi","CLI_pni","CLI_til","CLI_age","CLI_gender","CLI_histology","CLI_vitalstatus","SMG_mutsig_KRAS",
             "SMG_mutsig_TP53","SMG_mutsig_NRAS","SMG_mutsig_BRAF","SMG_mutsig_PIK3CA",
             "SMG_mutsig_APC","SMG_mutsig_FBXW7",
             "SMG_mutsig_PTEN","SMG_mutsig_TCF7L2","SMG_mutsig_SOX9","SMG_mutsig_CTNNB1",
             "CLI_stagenum",
             "CLI_grade","CLI_mmr","rate_non"))  
  
  imgf <- c("CD3_4_8_20r","PD1_20r","PanCK_Ecad_20r","PDL1_20r","CD68_20r",
  "CD163_20r","CD68plusPDL1_20r","CD3plusPD1_20r","PanCKplusPDL1_20r","CD3_4_8_200r",
  "PD1_200r","PanCK_Ecad_200r","PDL1_200r","CD68_200r",
  "CD163_200r","CD68plusPDL1_200r","CD3plusPD1_200r","PanCKplusPDL1_200r")
  
  setnames(dt, imgf, paste0("CLI_",imgf))
  
  # tmp
  setnames(dt, "CLI_CD3_4_8_20r","cd3_non")
  setnames(dt, "CLI_CD3_4_8_20r","cd3_non")
  
  
  # scale for the plot
  dt[, rate_non := rate_non * 1e-6]

  # Create a named list for each row
  dt_list <- lapply(1:nrow(dt), function(i) {
    setNames(as.list(dt[i, -"ID", with = FALSE]), names(dt)[-which(names(dt)=="ID")])
  })
  
  # format the SMGs
  dt_list <- lapply(dt_list, function(sub_list) {
    
    # Iterate over each sublist
    sub_list <- sapply(names(sub_list), function(name) {
      
      # If name starts with "SMG"
      if (grepl("^SMG", name)) {
        
        # If value is "NA", remove element
        if (sub_list[[name]] == "NA") {
          return(NULL)
        }
        
        # Else, convert to numeric
        else {
          return(as.numeric(sub_list[[name]]))
        }
      }
      
      # If name does not start with "SMG", return element as is
      else {
        return(sub_list[[name]])
      }
    }, simplify = FALSE)  # Do not simplify results to a vector
    
    # Remove NULL elements (those that were "NA")
    sub_list[sapply(sub_list, is.null)] <- NULL
    
    return(sub_list)
  })
  

  # Name each list with the corresponding `ID`
  names(dt_list) <- c(dt$`ID`)
  
  # SMG list
  smg_list <- c("SMG_mutsig_KRAS", "SMG_mutsig_TP53", "SMG_mutsig_NRAS", 
                "SMG_mutsig_BRAF", "SMG_mutsig_PIK3CA", "SMG_mutsig_APC", 
                "SMG_mutsig_FBXW7", "SMG_mutsig_PTEN", "SMG_mutsig_TCF7L2", 
                "SMG_mutsig_SOX9", "SMG_mutsig_CTNNB1")
  # Initialize an empty list
  dt_list[["all_q"]] <- list()
  
  # Loop through the names
  for (name in smg_list) {
    # Perform calculation and add to list
    suppressWarnings(dt_list[["all_q"]][[name]] <- sum(as.numeric(dt[[name]]), na.rm=TRUE) / sum(!is.na(dt[[name]])))
  }
  
  # Create a list of samples
  samples <- as.vector(dt$`ID`)
  
  # Create a list of sets (column names)
  sets <- names(dt)[-which(names(dt)=="ID")]
  
  # Convert the list to JSON
  json_data <- list("data" = dt_list, "samples" = samples, "sets" = sets)
  
  # Convert to JSON and save to file
  json_file_name <- "~/git/xyz/data/fireBrowseJson/JW-TP.coMut_table.json"
  writeLines(toJSON(json_data, pretty = TRUE, auto_unbox = TRUE), json_file_name)
  
}

load_radials <- function() {
  
  mydir <- "~/Sorger/orion/processed_40_orion/header_pheno_spatial/"
  files <- list.files(path = mydir, pattern = "*radial.csv", full.names = TRUE)
  
  # Define the function to calculate the means and read into a data.table
  calculate_means <- function(file) {
    cat("working on", file, "\n")
    
    # 147456 is PanCK + Ecad
    # 2048 is PDL1
    # 4416 is CD3 or CD8 or CD4
    # Create the command string
    cmd <- paste0("~/git/cysift/src/cysift select -o 4416 - <", file, " | awk -F, '{ for (i=NF-17; i<=NF; i++) sum[i] += $i; ++n} END { for (i=NF-17; i<=NF; i++) printf \"%f%s\", sum[i]/n, (i==NF ? \"\\n\" : \",\") }'")
    #cmd <- paste0("grep -v ^@ ", file, " | awk -F, '{ for (i=NF-17; i<=NF; i++) sum[i] += $i; ++n} END { for (i=NF-17; i<=NF; i++) printf \"%f%s\", sum[i]/n, (i==NF ? \"\\n\" : \",\") }'")
    
    # Execute the command and read the output into a data.table
    tmp <- fread(cmd = cmd)
    
    # Add the file name as a new column
    tmp[, "file_name" := file]
    
    return(tmp)
  }
  
  dt_list <- mclapply(files, calculate_means, mc.cores = 4)
  dt <- rbindlist(dt_list, fill = TRUE)
  
  setnames(dt, paste0("V",seq(18)), c("CD3_4_8_20r","PD1_20r","PanCK_Ecad_20r",
                                      "PDL1_20r","CD68_20r","CD163_20r",
                                      "CD68+PDL1_20r","CD3+PD1_20r","PanCK+PDL1_20r",
                                      "CD3_4_8_200r","PD1_200r","PanCK_Ecad_200r",
                                      "PDL1_200r","CD68_200r","CD163_200r","CD68+PDL1_200r",
                                      "CD3+PD1_200r","PanCK+PDL1_200r"))
  
  # Assuming dt is your data.table and file is the column with the file paths
  dt[, crcid := sub('.*(CRC[0-9]+).*', '\\1', file_name)]
  
  dt2 = data.table(crcid = paste0("CR",redcap$SpecimenID), ID = redcap$ID)
  dt3 <- dt[dt2, on = .(crcid), nomatch = 0]
  
  return(dt3)
}

header_read <- function(infile) {
    
    f <- fread(cmd=paste("cysift view -H",infile), sep="\t", header=FALSE, fill=TRUE)
    setnames(f, c("V1","V2","V3","V4"), c("tag_type","id","type","description"))
    f[, description := sub("^...", "", description)]
    f[, tag_type := sub("^.", "", tag_type)]
    f[, type := sub("^...", "", type)]
    f[, id := sub("^...", "", id)]
    return(f)
}

ttest_pairs <- function(dt, group_col, data_col) {
  # Extract unique group levels
  groups <- unique(dt[[group_col]])
  
  # Create a list to store the t.test results for each pair
  results <- list()
  
  results <- data.table()
  # Loop through each unique pair of groups and perform t.test
  for (i in 1:(length(groups) - 1)) {
    for (j in (i + 1):length(groups)) {
      group1 <- dt[get(group_col) == groups[i], ..data_col, with = FALSE]
      group2 <- dt[get(group_col) == groups[j], ..data_col, with = FALSE]
      test_result <- t.test(group1[[1]], group2[[1]])
      w_result <- wilcox.test(group1[[1]],group2[[1]])
      
      results <- rbind(results, data.table(
        a=as.character(groups[i]),
        b=as.character(groups[j]),
        p_t=test_result$p.value,
        p_w=w_result$p.value,
        mean_a = mean(group1[[1]], na.rm=TRUE),
        mean_b = mean(group2[[1]], na.rm=TRUE)
      ))
    }
  }

  return(results)
}

library(boot)
# Function for bootstrapping
# Function for bootstrapping
bootstrapR2 <- function(data, indices, covariates) {
  covariates <- unique(c(covariates, "PFSDays", "PFSC"))
  
  # Subset the data table correctly for bootstrapping
  data_subset <- data[indices, ..covariates]
  cox_model <- suppressWarnings(coxph(Surv(PFSDays, PFSC) ~ ., data = data_subset, control = coxph.control(iter.max = 1000)))
  return(coxsnellR2(cox_model))
}

coxsnellR2 <- function(cox_model) {
  
  # Ensure that the input is a Cox model
  if (!inherits(cox_model, "coxph")) {
    stop("Input must be a coxph model")
  }

  loglik_full <- cox_model$loglik[2]  # log-likelihood of the full model
  loglik_null <- cox_model$loglik[1]  # log-likelihood of the null model
  
  # Number of observations
  n <- length(cox_model$residuals)
  
  # Calculate Likelihood Ratio Test Statistic
  lrt_statistic <- -2 * (loglik_null - loglik_full)
  
  # Calculate Cox & Snell R-squared
  cox_snell_r2 <- 1 - exp(-lrt_statistic / n)
  
  return(cox_snell_r2)
  
  # Calculate Nagelkerke R-squared
  max_cox_snell <- 1 - exp(-2 / n * loglik_null)
  nagelkerke_r2 <- cox_snell_r2 / max_cox_snell
  
  # Return Nagelkerke R-squared
  return(nagelkerke_r2)
}


partition_vector <- function(vec) {
  total_sum <- sum(vec)
  n <- length(vec)
  half_sum <- total_sum %/% 2
  
  # Create and initialize DP matrix
  dp <- matrix(FALSE, n + 1, half_sum + 1)
  dp[,1] <- TRUE
  
  # Fill the DP matrix
  for (i in 1:n) {
    for (j in 1:half_sum) {
      if (j >= vec[i]) {
        dp[i + 1, j + 1] <- dp[i, j + 1] || dp[i, j - vec[i] + 1]
      } else {
        dp[i + 1, j + 1] <- dp[i, j + 1]
      }
    }
  }
  
  # Find the closest sum to half_sum
  for (j in half_sum:1) {
    if (dp[n + 1, j + 1]) {
      closest_sum <- j
      break
    }
  }
  
  # Backtrack to find the elements of the first group
  group1 <- c()
  temp_sum <- closest_sum
  for (i in n:1) {
    if (!dp[i, temp_sum + 1]) {
      group1 <- c(group1, vec[i])
      temp_sum <- temp_sum - vec[i]
    }
  }
  
  # The other group is formed by the remaining elements
  group2 <- setdiff(vec, group1)
  
  list(Group1 = group1, Group2 = group2)
}

library(data.table)

compare_groups <- function(rrt, var) {
  # Ensure rrt is a data.table
  setDT(rrt)
  
  # Define the pairs for comparison
  pairs <- list(c("tipMMR", "pMMR"), c("dMMR", "pMMR"), c("dMMR", "tipMMR"))
  
  # Initialize an empty data.table for results
  results <- data.table(comparison = character(), 
                        wilcox_p_value = numeric(), 
                        t_test_p_value = numeric(),
                        mean_group1 = numeric(), 
                        median_group1 = numeric(),
                        mean_group2 = numeric(), 
                        median_group2 = numeric())
  
  # Perform Wilcoxon and t-tests for each pair and calculate mean and median
  for (pair in pairs) {
    group1 <- rrt[get("tipMMR") == pair[1], get(var)]
    group2 <- rrt[get("tipMMR") == pair[2], get(var)]
    
    # Wilcoxon test
    wilcox_test <- wilcox.test(group1, group2)
    
    # t-test
    t_test <- t.test(group1, group2)
    
    # Calculate mean and median
    mean1 <- mean(group1, na.rm = TRUE)
    median1 <- median(group1, na.rm = TRUE)
    mean2 <- mean(group2, na.rm = TRUE)
    median2 <- median(group2, na.rm = TRUE)
    
    # Append results
    results <- rbindlist(list(results, data.table(
      comparison = paste(pair[1], "vs", pair[2]),
      wilcox_p_value = wilcox_test$p.value,
      t_test_p_value = t_test$p.value,
      mean_group1 = mean1,
      median_group1 = median1,
      mean_group2 = mean2,
      median_group2 = median2
    )))
  }
  
  return(results)
}

# Function to create a beeswarm plot with significance comparisons
create_beeswarm_plot <- function(data, x_var, y_var, withbox=FALSE) {
  dx_combos <- list(c("tipMMR","tdpMMR"),c("tipMMR","dMMR"),c("dMMR","tdpMMR"))
 # dx_combos <- list(c("focal","no"))
  #dx_combos <- list(c("Pre","Post"))
  #dx_combos  <- list(c("PD","PR"),c("PD","SD"),c("SD","PR"))
  #dx_combos <- list(c("High","Low"))
  #dx_combos <- list(c("Mut","WT"))
  
  g <- ggplot(data, aes_string(x = x_var, y = y_var)) 
  
  g <- g + ggbeeswarm::geom_quasirandom(
      shape = 21, color = "white",
      alpha = 1, size = 1.5,
      aes(fill = .data[[x_var]])
    ) +
    #scale_fill_manual(values = jhu_colors) + 
    #scale_fill_manual(values = pMMR_colors) +
    #scale_fill_manual(values = gleason_colors) +
    theme_bw() +
    theme(
      legend.position = "none",
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.title = element_text(size = 10, face = "bold"),
      axis.text = element_text(size = 10)
    ) +
    stat_compare_means(comparisons = dx_combos, label = "p.signif",
                       method="wilcox.test") # Add significance levels

  if (withbox)
    g <- g + geom_boxplot(width = 0.2,                                                                                                                                                  
                          outlier.shape = NA,                                                                                                                                           
                          alpha = 0.5,
                          fill=NA,
                          color="black")
  g
}

sig_display <- function(p) {
  
  if (p > 0.05)
    return ("NS")
  if (p > 0.01)
    return ("*")
  if (p > 0.001)
    return("**")
  return("***")
}