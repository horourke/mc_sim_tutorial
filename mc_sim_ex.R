# Example script for O'Rourke (2023-2024) Monte Carlo simulation tutorial for path models

install.packages("MplusAutomation")
install.packages("dplyr")
library(MplusAutomation)
library(dplyr)

## Create master directory for simulation study ##
# Set up directories to switch easily between home, work, etc. - choose 1
#homepath <- "C:/users/home/mc_sim_ex/"
#filepath <- homepath
workpath <- "C:/users/work/mc_sim_ex/"
filepath <- workpath

# Create master directory on your machine if needed
if (!dir.exists(filepath)) {
  dir.create(filepath, recursive = TRUE)
}

# Conditions to iterate over
# factors = outcome distribution type
# n_values = sample sizes
factors <- c("ZINB", "ZIP")
n_values <- c(100, 250, 500, 750, 1500)
reps <- 500
dist <- "pi"
yvarspec <- "y"

# Nested loop over factors and values
# The sets of if statements account for the differing numbers of parameters 
# in the ZIP and ZINB output files.
# ZINB has 1 more parameter, the dispersion parameter (variance of Y)
for (f in factors) {
  # This loop allows the .inp scripts to have different text for the ZINB and ZIP models
  # ZINB has "y(nbi)" and ZIP has "y(pi)"
  # The variance of y ("y") is commented out of code (!) for ZIP models
  if (f == "ZINB") {
    dist <- "nbi"
    yvarspec <- "y"
  } else if (f == "ZIP") {
    dist <- "pi"
    yvarspec <- "!y"
  }
  for (n in n_values) {
    
    # Create the subdirectory path for each condition
    sub_path <- paste0(filepath, f, "/N = ", n)
    
    #Create actual subdirectory if not yet created
    if (!dir.exists(sub_path)) {
      dir.create(sub_path, recursive = TRUE)
    }
    
    ###############################    
    # Write data generation scripts
    ###############################  
    
    # Using writeLines to write text lines to a .inp file
    # sprintf() returns character objects containing a formatted combination of input values
    # sprintf() allows us to refer to our objects and factors within the loops
    # Within the sprintf function, %s refers to string value (text), %d refers to digit (numeric value)
    # the final line in sprintf() gives the order of the objects that are being referred to
    
    # Create the file to be written
    dgscr <- paste0(sub_path, "/", f, "_n", n, ".inp")
    writeLines(sprintf('title: MC sim example for %s n=%d;
montecarlo:			
	names = x m y;
	seed = 53487;
	nobs = %d;
	nreps = %d;
	generate = y(%s);
    COUNT = y(%s);
    REPSAVE = ALL;
	save = %s%d_*.dat;
ANALYSIS: 
estimator=ml;
integration=montecarlo;

model population:
        [m@0];
        m@1;
        [x@0];
        x@0.25;
  		y ON x*.01 m*0.14;
   		[y@3];
        %s@1;
        y#1 ON x*-.01 m*-0.14;
  		[y#1@0];
        m ON x*0.59;
MODEL:
        [m] (m);
        m*1;
   		[y*3] (int_c);
        %s*1;
  		[y#1*0] (int_z);
  		y ON x*.01 (cpc);
        y ON m*0.14 (bc);
        y#1 ON x*-.01 (cpz);
        y#1 ON m*-0.14 (bz);
        m ON x*0.59 (a);
OUTPUT: TECH9;', 
                       f, n, n, reps, dist, dist, f, n, yvarspec, yvarspec),
               con = dgscr)
    
    ## Simulate data ##
    runModels(sub_path)
    
    ## Create file suffix naming conventions for moving the files around ##
    # file suffix for the data gen .inp/.out files
    file_suffix_n <- paste0(f, "_n", n)
    # file suffix for the "list.dat" files that are created during runModels()
    file_suffix_list <- paste0(tolower(f), n, "_list")
    
    # Move data gen .inp/.out & list files to "datagen" folder before creating analysis scripts
    datagen <- paste0(filepath, "datagen")
    if (!dir.exists(datagen)) {
      dir.create(datagen, recursive = TRUE)
    }
    file.copy(paste0(sub_path, "/", file_suffix_n, ".inp"), datagen)
    file.remove(paste0(sub_path, "/", file_suffix_n, ".inp"))
    file.copy(paste0(sub_path, "/", file_suffix_n, ".out"), datagen)
    file.remove(paste0(sub_path, "/", file_suffix_n, ".out"))
    file.copy(paste0(sub_path, "/", file_suffix_list, ".dat"), datagen)
    file.remove(paste0(sub_path, "/", file_suffix_list, ".dat"))
    
    ########################################################################
    # Write & save scripts that create analysis scripts for each replication
    ########################################################################
    
    # Create the file to be written; ensure suffix matches to code below
    ascr <- paste0(sub_path, "/", f, n, "_r.inp")
    
    # Write the Mplus script to the file
    writeLines(sprintf('[[init]]
iterators = sample;
sample = 1:%d;
filename = "%s_n%d_[[sample]].inp";
outputDirectory = "%s";
[[/init]]

TITLE: MC sim example for %s n=%d rep [[sample]];
DATA:
    File is %s%d_[[sample]].dat;

  VARIABLE:
      NAMES = y m x;
      USEVARIABLES = y m x;
      count = y(%s);
ANALYSIS:
  estimator=ml;
  integration=montecarlo;
  bootstrap=500;
MODEL:
          [m] (m);
          m*1;
            [y*3] (int_c);
          %s*1;
            [y#1*0] (int_z);
            y ON x*.01 (cpc);
          y ON m*0.14 (bc);
          y#1 ON x*-.01 (cpz);
          y#1 ON m*-0.14 (bz);
          m ON x*0.59 (a);

  MODEL CONSTRAINT:
  NEW(ab_ll_x0 ab_lg_x0 ab_ll_x1 ab_lg_x1);

  ab_ll_x0 = a*bz*exp(int_z + bz*m + cpz*0);
  ab_lg_x0 = a*bz*(exp(int_z + bz*m + cpz*0)/((1+exp(int_z + bz*m + cpz*0))^2));

  ab_ll_x1 = a*bz*exp(int_z + bz*m + cpz*1);
  ab_lg_x1 = a*bz*(exp(int_z + bz*m + cpz*1)/((1+exp(int_z + bz*m + cpz*1))^2));
  OUTPUT: CINTERVAL(bcbootstrap);
  OUTPUT: TECH1;', 
                       reps, f, n, sub_path, f, n, f, n, dist, yvarspec),
               con = ascr)
    
    # file suffix for the creation scripts 
    file_suffix_r <- paste0(tolower(f), n, "_r") 
    
    ## Create .inp files for all replications in a subdirectory ##
    createModels(paste0(sub_path, "/", file_suffix_r, ".inp"))
    
    #Move createmodels .inp files out of subdirectory before running model scripts
    createm <- paste0(filepath,"/modelcreation")
    if (!dir.exists(createm)) {
      dir.create(createm, recursive = TRUE)
    }
    file.copy(paste0(sub_path, "/", file_suffix_r, ".inp"), createm)
    file.remove(paste0(sub_path, "/", file_suffix_r, ".inp")) 
    
    ## Run model scripts for all replications in subdirectory ##
    runModels(sub_path)
  }
}
  
#############################################
# SAVING OUTPUT ESTIMATES FOR FUTURE ANALYSES
#############################################

# This code saves raw parameter estimates, standard errors, p values, and confidence intervals

# Re-initialize loops because we need the attached if statements
# "matcols" refers to the # of columns of estimates from the Mplus tech1 output
# ZINB models have 1 more column for the Y dispersion parameter
for (f in factors) {
  if (f == "ZINB") {
    matcols <- 14
  } else if (f == "ZIP") {
    matcols <- 13
  }
  for (n in n_values) {
    # Create references to condition subdirectories & file replications within these loops
    sub_path <- paste0(filepath, f, "/N = ", n)
    file_suffix_n <- paste0(f, "_n", n)
    
    #Use MplusAuto to read in results from all output files in the referenced directory
    # Outputs are read in as lists by condition, with list length = reps
    output <- paste0(f, "_", n, "_output")
    assign(output, readModels(
      target = sub_path,  
      recursive = TRUE,
      what = "parameters",
      quiet = TRUE)
    )
    
    #Create empty matrix with for all replications by condition
    # Create one per desired set of output
    mat_ests <- matrix(0, nrow = reps, ncol = matcols)
    mat_se <- matrix(0, nrow = reps, ncol = matcols)
    mat_p <- matrix(0, nrow = reps, ncol = matcols)
    mat_l <- matrix(0, nrow = reps, ncol = matcols)
    mat_u <- matrix(0, nrow = reps, ncol = matcols)
    
    # f1-f3 convert filepath to list naming conventions used by readModels()
    f1 <- gsub(":",".",filepath)
    f2 <- gsub("/",".",f1)
    f3 <- gsub(" ",".",f2)
    
    # Use list data to create data frames containing output from each replication by condition
    # get() pulls values from the list
    # the following line converts the list values to a numeric object, 
    #then to a vector, then transposes the vector such that rows are replications
    # the third line referring to the matrix adds each replication to a new row in the matrix
    # This process is repeated for each set of desired output
    for (i in 1:reps) {
      filename <- paste0(f3, f, ".N...", n, ".", tolower(f), "_n", n, "_", i, ".out")
      ## Parameter estimates ##
      ests_1 <- get(output)[[filename]][["parameters"]][["unstandardized"]][["est"]]
      ests_2 <- t(as.matrix(unlist(ests_1)))
      mat_ests[i, ] <- as.numeric(ests_2)
      ## Standard errors ##
      se_1 <- get(output)[[filename]][["parameters"]][["unstandardized"]][["se"]]
      se_2 <- t(as.matrix(unlist(se_1)))
      mat_se[i, ] <- as.numeric(se_2)
      ## p values ##
      pval_1 <- get(output)[[filename]][["parameters"]][["unstandardized"]][["pval"]]
      pval_2 <- t(as.matrix(unlist(pval_1)))
      mat_p[i, ] <- as.numeric(pval_2)
      ## Lower CIs ##
      low_1 <- get(output)[[filename]][["parameters"]][["ci.unstandardized"]][["low2.5"]]
      low_2 <- t(as.matrix(unlist(low_1)))
      mat_l[i, ] <- as.numeric(low_2)
      ## Upper CIs ##
      up_1 <- get(output)[[filename]][["parameters"]][["ci.unstandardized"]][["up2.5"]]
      up_2 <- t(as.matrix(unlist(up_1)))
      mat_u[i, ] <- as.numeric(up_2)
    }
    #Create data frames from matrices
    ests <- as.data.frame(mat_ests)
    se <- as.data.frame(mat_se)
    pval <- as.data.frame(mat_p)
    cilow <- as.data.frame(mat_l)
    ciup <- as.data.frame(mat_u)
    
    # Save each output as file in condition subdirectories
    write.table(ests, file=paste0(sub_path, "/", file_suffix_n, "_ests.dat"), row.names=FALSE, sep="\t", quote=FALSE)
    write.table(se, file=paste0(sub_path, "/", file_suffix_n, "_stderrs.dat"), row.names=FALSE, sep="\t", quote=FALSE)
    write.table(pval, file=paste0(sub_path, "/", file_suffix_n, "_pvals.dat"), row.names=FALSE, sep="\t", quote=FALSE)
    write.table(cilow, file=paste0(sub_path, "/", file_suffix_n, "_cis_lower.dat"), row.names=FALSE, sep="\t", quote=FALSE)
    write.table(ciup, file=paste0(sub_path, "/", file_suffix_n, "_cis_upper.dat"), row.names=FALSE, sep="\t", quote=FALSE)
  }
}

## Combine condition-level data into one dataset per outcome with all replications, for all conditions ##

# Create empty data frame for each output
simests <- data.frame()
ses <- data.frame()
pvs <- data.frame()
cisl <- data.frame()
cisu <- data.frame()

# Row counter
row1 <- 1

# Re-initialize the factor loops, because we no longer need the "if" statements in the loops above
for (f in factors) {
  for (n in n_values) {
    sub_path <- paste0(filepath, f, "/N = ", n)
    file_suffix_n <- paste0(f, "_n", n)
    ## Parameter estimates ##
    filename_est <- paste0(sub_path, "/", f, "_n", n, "_ests.dat")
    esttable <- read.table(filename_est, header = TRUE)
    simests <- bind_rows(simests, esttable)
    ## Standard errors ##
    filename_se <- paste0(sub_path, "/", f, "_n", n, "_stderrs.dat")
    setable <- read.table(filename_se, header = TRUE)
    ses <- bind_rows(ses, setable)
    ## p values ##
    filename_p <- paste0(sub_path, "/", f, "_n", n, "_pvals.dat")
    ptable <- read.table(filename_p, header = TRUE)
    pvs <- bind_rows(pvs, ptable)
    ## Lower CIs ##
    filename_cil <- paste0(sub_path, "/", f, "_n", n, "_cis_lower.dat")
    ciltable <- read.table(filename_cil, header = TRUE)
    cisl <- bind_rows(cisl, ciltable)
    # Upper CIs ##
    filename_ciu <- paste0(sub_path, "/", f, "_n", n, "_cis_upper.dat")
    ciutable <- read.table(filename_ciu, header = TRUE)
    cisu <- bind_rows(cisu, ciutable)
    
    # create factor & rep variables for each dataset
    for (rep in 1:reps) {
      simests[row1, "rep"] <- rep
      simests[row1, "factor"] <- f
      simests[row1, "n"] <- n
      ses[row1, "rep"] <- rep
      ses[row1, "factor"] <- f
      ses[row1, "n"] <- n
      pvs[row1, "rep"] <- rep
      pvs[row1, "factor"] <- f
      pvs[row1, "n"] <- n
      cisl[row1, "rep"] <- rep
      cisl[row1, "factor"] <- f
      cisl[row1, "n"] <- n
      cisu[row1, "rep"] <- rep
      cisu[row1, "factor"] <- f
      cisu[row1, "n"] <- n
      # Add a row to the dataframe
      row1 <- row1 + 1
    }
  }
}

# This splitting step is not always necessary - split the data to correctly arrange variable names
# this must be done due to where the dispersion parameter is located in the Mplus Tech1 output for the ZINB models
# if you always have the same number of parameter estimates across conditions, skip this step
# You WILL still want to use the names() function to assign names to your dataset, use Mplus tech1 to do so
split_data_ests <- split(simests, simests$factor)
split_data_ests$ZIP <- split_data_ests$ZIP[,c(1:9,14,10:13,15:17)]
names(split_data_ests$ZIP) <- c("cp_counts", "b_counts", "cp_zeroes", "b_zeroes", 
                                "a", "m_int", "y_z_int", "y_c_int", 
                                "m_resvar", "dispers", 
                                "ab_ll_x0", "ab_lg_x0", "ab_ll_x1", "ab_lg_x1", 
                                "rep", "factor", "n")
names(split_data_ests$ZINB) <- c("cp_counts", "b_counts", "cp_zeroes", "b_zeroes", 
                                 "a", "m_int", "y_z_int", "y_c_int", 
                                 "m_resvar", "dispers", 
                                 "ab_ll_x0", "ab_lg_x0", "ab_ll_x1", "ab_lg_x1", 
                                 "rep", "factor", "n")
simests <- bind_rows(split_data_ests)

split_data_se <- split(ses, ses$factor)
split_data_se$ZIP <- split_data_se$ZIP[,c(1:9,14,10:13,15:17)]
names(split_data_se$ZIP) <- c("cp_counts_se", "b_counts_se", "cp_zeroes_se", "b_zeroes_se", 
                              "a_se", "m_int_se", "y_z_int_se", "y_c_int_se", 
                              "m_resvar_se", "dispers_se", 
                              "ab_ll_x0_se", "ab_lg_x0_se", "ab_ll_x1_se", "ab_lg_x1_se", 
                              "rep", "factor", "n")
names(split_data_se$ZINB) <- c("cp_counts_se", "b_counts_se", "cp_zeroes_se", "b_zeroes_se", 
                               "a_se", "m_int_se", "y_z_int_se", "y_c_int_se", 
                               "m_resvar_se", "dispers_se", 
                               "ab_ll_x0_se", "ab_lg_x0_se", "ab_ll_x1_se", "ab_lg_x1_se", 
                               "rep", "factor", "n")
ses <- bind_rows(split_data_se)

split_data_p <- split(pvs, pvs$factor)
split_data_p$ZIP <- split_data_p$ZIP[,c(1:9,14,10:13,15:17)]
names(split_data_p$ZIP) <- c("cp_counts_p", "b_counts_p", "cp_zeroes_p", "b_zeroes_p", 
                             "a_p", "m_int_p", "y_z_int_p", "y_c_int_p", 
                             "m_resvar_p", "dispers_p", 
                             "ab_ll_x0_p", "ab_lg_x0_p", "ab_ll_x1_p", "ab_lg_x1_p", 
                             "rep", "factor", "n")
names(split_data_p$ZINB) <- c("cp_counts_p", "b_counts_p", "cp_zeroes_p", "b_zeroes_p", 
                              "a_p", "m_int_p", "y_z_int_p", "y_c_int_p", 
                              "m_resvar_p", "dispers_p", 
                              "ab_ll_x0_p", "ab_lg_x0_p", "ab_ll_x1_p", "ab_lg_x1_p", 
                              "rep", "factor", "n")
pvs <- bind_rows(split_data_p)

split_data_cil <- split(cisl, cisl$factor)
split_data_cil$ZIP <- split_data_cil$ZIP[,c(1:9,14,10:13,15:17)]
names(split_data_cil$ZIP) <- c("cp_counts_l", "b_counts_l", "cp_zeroes_l", "b_zeroes_l", 
                               "a_l", "m_int_l", "y_z_int_l", "y_c_int_l", 
                               "m_resvar_l", "dispers_l", 
                               "ab_ll_x0_l", "ab_lg_x0_l", "ab_ll_x1_l", "ab_lg_x1_l", 
                               "rep", "factor", "n")
names(split_data_cil$ZINB) <- c("cp_counts_l", "b_counts_l", "cp_zeroes_l", "b_zeroes_l", 
                                "a_l", "m_int_l", "y_z_int_l", "y_c_int_l", 
                                "m_resvar_l", "dispers_l", 
                                "ab_ll_x0_l", "ab_lg_x0_l", "ab_ll_x1_l", "ab_lg_x1_l", 
                                "rep", "factor", "n")
cilow <- bind_rows(split_data_cil)

split_data_ciu <- split(cisu, cisu$factor)
split_data_ciu$ZIP <- split_data_ciu$ZIP[,c(1:9,14,10:13,15:17)]
names(split_data_ciu$ZIP) <- c("cp_counts_u", "b_counts_u", "cp_zeroes_u", "b_zeroes_u", 
                               "a_u", "m_int_u", "y_z_int_u", "y_c_int_u", 
                               "m_resvar_u", "dispers_u", 
                               "ab_ll_x0_u", "ab_lg_x0_u", "ab_ll_x1_u", "ab_lg_x1_u", 
                               "rep", "factor", "n")
names(split_data_ciu$ZINB) <- c("cp_counts_u", "b_counts_u", "cp_zeroes_u", "b_zeroes_u", 
                                "a_u", "m_int_u", "y_z_int_u", "y_c_int_u", 
                                "m_resvar_u", "dispers_u", 
                                "ab_ll_x0_u", "ab_lg_x0_u", "ab_ll_x1_u", "ab_lg_x1_u", 
                                "rep", "factor", "n")
ciupper <- bind_rows(split_data_ciu)

#Save estimates data as one large .dat file in the filepath main directory
write.table(simests, file=paste0(filepath, "/all_ests.dat"), row.names=FALSE, sep="\t", quote=FALSE)
write.table(ses, file=paste0(filepath, "/stderrs.dat"), row.names=FALSE, sep="\t", quote=FALSE)
write.table(pvs, file=paste0(filepath, "/pvals.dat"), row.names=FALSE, sep="\t", quote=FALSE)
write.table(cilow, file=paste0(filepath, "/cilow.dat"), row.names=FALSE, sep="\t", quote=FALSE)
write.table(ciupper, file=paste0(filepath, "/ciupper.dat"), row.names=FALSE, sep="\t", quote=FALSE)

