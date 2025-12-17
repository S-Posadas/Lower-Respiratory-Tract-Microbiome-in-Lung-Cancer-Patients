#' Modified from Rhea package (Serial-Group-Comparisons):
#' Original Script: Statistical comparison of taxonomic composition and alpha-diversity parameters between groups

    
    ##################################################################################
    ######              Script for unpaired data with 2 groups                  ######
    ##################################################################################
    
    ###################       Load all required libraries     ########################
    
    # Check if required packages are already installed, and install if missing
    packages <-c("plotrix", "grid","ggrepel","gtable","Matrix", "PMCMRplus", "PerformanceAnalytics") 

    # packages <-c("plotrix","PerformanceAnalytics","reshape","ggplot2","gridExtra","grid","ggrepel","gtable","Matrix","cowplot", "tidyr","PMCMRplus") 
    
    # Function to check whether the package is installed
    # InsPack <- function(pack)
    # {
    #   if ((pack %in% installed.packages()) == FALSE) {
    #     install.packages(pack,repos ="http://cloud.r-project.org/")
    #   } 
    # }
    
    # Applying the installation on the list of packages
    # lapply(packages, InsPack)
    
    # Make the libraries 
    lib <- lapply(packages, require, character.only = TRUE)
    
    # Check if it was possible to install all required libraries
    flag <- all(as.logical(lib))
    
    # Resolve dplyr conflicts explicitly
    conflictRules("dplyr", exclude = "lag")  # Prefer stats::lag()
    
    #####################################################################################################################
    ####                                        Functions to be used  in main Script.                            ########
    #####################################################################################################################
    
    # Function to calculate relative abundance -> not necessary for this script
    # rel.abundance <- function(data)
    # {
    #   total = sum(data)
    #   rel.data <- 100 * data / total
    #   return(rel.data)
    # }
    
    # Replace abundance with zero if value is below given cutoff
    abundance.fix <- function(data)
    {
      data[data < abundance_cutoff] <- 0
      return(data)
    }
    
    # Replace zero value with NA 
    fill_zero.NA <- function(data,ReplaceZero)
    {
      if (ReplaceZero == "NO") {
        return(data)
      } else if (ReplaceZero == "YES") {
        data[data == 0] <- NA
        return(data)
      } else {
        return(data)
      }
    }
    
    # Return maxima and minima values of the given input value
    max.fun <- function(data)
    {
      data.max <- max(as.numeric(as.character(data)), na.rm = TRUE)
      return (data.max)
    }
    
    # Calculate the prevalence of the given input table
    pre.fun.na <- function(data)
    {
      prevalence <- nnzero(data, na.counted = FALSE)
      return(prevalence)
    }
    
    # Return the maximum median value for each group
    max.med <- function(data)
    {
      max.median <-max(aggregate (data ~ independent_variable,FUN = median,simplify = TRUE)[,2])
      return(max.median)
    }
    
    # Return the prevalence ratio
    max.pre <- function(data)
    {
      
      # Return the number of samples (excluding NA) for each group separately
      found <- aggregate (data ~ independent_variable,FUN = pre.fun.na,simplify = TRUE,na.action = na.pass)[,2]
      
      # Return the total number of samples (including missing values) for each group separately
      all <- aggregate (data ~ independent_variable,FUN = base::length,simplify = TRUE,na.action = na.pass)[,2]
      
      # Calculate the ratio for each group and return the maximum ratio out of the groups
      max.ratio <- max(found / all)
      return(max.ratio)
    }
    
    # Set the theme to change text for plotting (ggplot - Gtable)
    mytheme <- gridExtra::ttheme_default(

      # Adjust settings for the text inside table
      core = list(fg_params = list(cex = 0.8)),

      # Adjust the test for column and row header
      colhead = list(fg_params = list(cex = 0.9)),
      rowhead = list(fg_params = list(cex = 1.0))
    )

    
    #####################################################################################################################
    ####                                      Pre-processing of OTUs Table                                       ########
    #####################################################################################################################
    
    # Convert independent variable into factor to avoid errors
    original_table[,independant_variable_name] <- as.factor(original_table[,independant_variable_name]) 
    
    # Store independent variable columns from original table 
    independent_variable <- original_table[[independant_variable_name]]
    
    # Check if group_order includes all groups if not the default order will be used
    if(dim(as.data.frame(group_order))[1]!=nlevels(independent_variable)) { 
      group_order <- levels(independent_variable)
    }
    
    # Store columns with sample ID information from original table 
    if(isTRUE(paired)) id_variable = as.factor(original_table[[id_name]])
    
    # Store metadata variable columns from original table 
    my_meta_data <- original_table[1:taxonomic_variables_start - 1]
    
    # Store relative abundance values of all OTUs
    my_otu_data <- original_table[taxonomic_variables_start:dim(original_table)[2]]
    my_otu_data <- my_otu_data[colSums(my_otu_data) != 0]

    # Transform data by zeroing very low abundances (based on given abundance cutoff - abundance_cutoff)
    my_otu_mod =  as.data.frame(apply(my_otu_data,2,abundance.fix))
    
    # Transform data by replacing all zero values with missing values
    # Column consisting of "NA" or "0" only are removed (see below)
    my_otu_mod_noz = as.data.frame(apply(my_otu_mod,2,fill_zero.NA,ReplaceZero))
    
    # Remove column if entire OTU column contain zeros or missing values
    my_otu_mod_noz <- my_otu_mod_noz[,!apply(my_otu_mod_noz , 2 , function(x) all(is.na(x) | (x == 0))), drop = FALSE]
    
    # Transform data by removing any OTU with median relative abundance below cutoff
    t_otu_mod_noz = as.data.frame(t(my_otu_mod_noz))
    
    # Calculate median for each OTUs
    t_otu_mod_noz$max.median <- apply(t_otu_mod_noz,1,max.med)
    
    # Select OTUs above median cutoff (med.cutoff)
    if(perform_max_median_cutoff == T){
      selected_max <- t_otu_mod_noz[t_otu_mod_noz$max.median >= max_median_cutoff,]
    }else{
      selected_max <- t_otu_mod_noz
    }
    
    # Remove calculated median column "max.median"
    selected_max$max.median <- NULL
    
    # Make a separate object as data frame (columns are OTUs and rows are samples)
    otu_mod_noz_max <- as.data.frame(t(selected_max))
    
    # Transpose data (columns are samples and rows are OTUs)
    t_otu_mod_noz_max = as.data.frame(t(otu_mod_noz_max))
    
    # Transform data by removing all OTUs with prevalence below the given cutoff
    t_otu_mod_noz_max$pre <- apply(t_otu_mod_noz_max,1,max.pre)
    selected_pre <- t_otu_mod_noz_max[t_otu_mod_noz_max$pre > prevalence_cutoff,]
    selected_pre$pre <- NULL
    
    # Transform and filter OTU table
    otu_mod_noz_max_pre <- as.data.frame(t(selected_pre))
    
    # Merge the metadata and OTU data in one table 
    # This table will be used as input for the analysis
    input_table <- cbind(my_meta_data,otu_mod_noz_max_pre) 
    
    # Only if paired
    if(isTRUE(paired)){
      
      # Check if all time points are available
      bool_vec <- table(input_table[[id_name]])==nlevels(as.factor(input_table[,dependant_variable_name]))
      
      # List with IDs without all time points
      bool_names <- names(subset(bool_vec,bool_vec==FALSE))
      
      length <- dim(input_table)[1]
      
      for( i in 1:length){
        # Check if time point information is needed for the ID
        if(input_table[i,id_name] %in% bool_names) {
          # Delete ID from the list if time point information was already added
          bool_names <- bool_names[bool_names != input_table[i,id_name]]
          id_sub <- subset(input_table, input_table[,id_name] == input_table[i,id_name])
          # Which timepoints are already available and which time points are missing
          missing <- levels(independent_variable)[which(levels(factor(input_table[[dependant_variable_name]],levels=group_order)) %in% id_sub[,dependant_variable_name]==FALSE)]
          # Add a line with all missing time points for the i-th ID
          for (j in 1:nlevels(as.factor(missing))) {
            # Add a row to to input_table
            input_table <- rbind(input_table,rep(NA,dim(input_table)[2]))
            # Add missing timepoint information
            input_table[dim(input_table)[1],dependant_variable_name] <- missing[j]
            # Assign ID to the missing time point
            input_table[dim(input_table)[1],id_name] <- input_table[i,id_name]
          }
        }
      }
      
      
      # Store repeated time variable columns from original table 
      independent_variable = factor(input_table[[dependant_variable_name]],levels=group_order)
      
      # Store columns with sample ID information from original table 
      id_variable = as.factor(input_table[[id_name]])
      
      # Convert independent variable into factor to avoid errors
      input_table[,dependant_variable_name] <- factor(input_table[,dependant_variable_name],levels=group_order) 
      
    }
    
    #####################################################################################################################
    ####                              Differential Statistical Analysis                                        ########
    #####################################################################################################################
    
    # The number of groups to compare
    groups_to_compare <- length(unique(input_table[,independant_variable_name]))
    if(groups_to_compare == 2){KW = F}else if(groups_to_compare > 2){KW = T}
    
    # The number of total observations per category (independant variable) to be used in the analysis
    total <- summary(independent_variable)
    
    # Create vector with group information 
    prevalence_list <- as.numeric(independent_variable)
    
    # Create vector with prevalence values
    # Iterate through all groups 
    for ( i in 1:nlevels(independent_variable)) {
      for ( j in 1:length(prevalence_list)) {
        
        # Assign the prevalence value of a group to each sample 
        if ( as.character(independent_variable[j]) == names(total[i])){
          prevalence_list[j] <- as.numeric(total[i])
        }
      }
    }
    
    # Create an empty dataframe for Kruskal-Wallis Rank Sum Test / Friedman Test
    if(KW == T) df <- data.frame(name = character(0),pvalue = numeric(0),sign = character(0))
    
    # Create an empty dataframe for Fisher's Exact Test
    Fdf <- data.frame(name = character(0),pvalue = numeric(0),sign = character(0))
    
    # Create an empty dataframe for Wilcoxon Rank Sum and Signed Rank Test
    all_pair_pval_table <- data.frame(measure = character(0),pair = character(0),Group1 = character(0),Group2 = character(0),pvalue = numeric(0),corrected = numeric(0))
    
    # Create an empty dataframe to the results of the paired Fisher's Exact Test
    all_pair_fpval_table <- data.frame(measure = character(0),pair = character(0),Group1 = character(0),Group2 = character(0),pvalue = numeric(0),corrected = numeric(0))
    
    ###################              Create empty lists to store information about signficant differences              ###################
    
    # Making a list for pairwise p-value table of Wilcoxon test
    pvaltable <- list()
    
    # Making a list for pairwise p-value table of Fisher's test
    fpvaltable <- list() #
    
    # Making a list for overall p-value table for Fisher's test
    allfpvaltable <- list()
    
    # Making a list for overall p-value table for Kruskal test
    pvaltableAll <- list()
    
    #######################               Friedman Rank Sum Test for each OTU for all Timepoints                     #######################
    #######################                   Skillings Mack Test for Missing Data                                   #######################
    
    fail <- FALSE
    
    if(isTRUE(paired) & isTRUE(KW)){
      
    # Start calculation with the dependant variable (e.g. Richness)
    for (i in dependant_variables_start:dim(input_table)[2])
    {
      
      # Take the values for all samples of the dependant variable/OTU
      my_test_vector <- input_table[,i]
      
      # Save the name of the observed variable/OTU
      my_name <- colnames(input_table)[i]
      test_table <- xtabs(my_test_vector[!is.na(my_test_vector)] ~ independent_variable[!is.na(my_test_vector)])
      
      # Number of groups to be compared
      num_of_represented_groups <- sum(as.vector(test_table) > 0)
      
      # Test whether any group is missing
      num_of_missing_groups <- nlevels(independent_variable) - num_of_represented_groups
      
      # Function to return "TRUE" or "FALSE" corresponding to "not missing" and "missing"
      if (num_of_missing_groups > 0) {
        isgroupmissing <- TRUE
      } else {
        isgroupmissing <- FALSE
      }
      
      # Empty dataframe with information for the analysis of repeated measurements
      friedman_df <- data.frame(value=double(), time=factor(),block=factor())
      
      # Save ID related values and time points in a dataframe
      tmp <- as.data.frame(cbind(my_test_vector,independent_variable,id_variable))
      
      # Create ordered dataframe for statistical analysis
      for(i in 1:nlevels(independent_variable)) {
        friedman_tmp <- subset(tmp,tmp[,2]==i)[order(subset(tmp,tmp[,2]==i)$id_variable),]
        friedman_df <- rbind(friedman_df,friedman_tmp)
      }
      
      # Save as matrix
      mat <- matrix(friedman_df[,1],nrow=nlevels(id_variable),ncol=nlevels(independent_variable))
      
      # Remove lines where no information is available for all time points
      mat <- mat[rowSums(is.na(mat)) !=  nlevels(independent_variable), ]
      
      if (ReplaceMissingValues == "NO"){
        # Remove lines with missing values completely
        mat <- mat[rowSums(is.na(mat))==0,]
        # Perfroms Firedman Test if more than one samples remain after removing missing values
        fit <- tryCatch (friedman.test(mat),error = function(i) {fail <<- TRUE})
      } else{
        # Performs a Friedman Rank Sum Test with Skillingsmack Test for missing Data
        fit <-  tryCatch (skillingsMackTest(mat),error = function(i) {fail <<- TRUE})
        
      }
      # Function to assign corrected and not-corrected pvalues  
      if (fail) {
        my_pvalue <- NaN
        fail <- FALSE
      } else {
        
        # Round the p-value down to four decimals 
        my_pvalue <- round(fit$p.value,8)
      }
      
      # Add the p-values to the table
      newRow <- data.frame(name = my_name,pvalue = my_pvalue,missing = isgroupmissing)
      df <- rbind(df,newRow)
    }
    
    # Applying Benjamini-Hochberg (1995) correction 
    df$corrected <- round(p.adjust(df$pvalue, method = "BH"),8)
    }
    
    #######################               Kruskal Wallis Test for each OTU for all Groups                #######################
    
    if(isFALSE(paired) & isTRUE(KW)){

      # Start calculation with the dependant variable (e.g. Richness)
      for (i in dependant_variables_start:dim(input_table)[2])
      {
        
        # Take the values for all samples of the dependant variable/OTU
        my_test_vector <- input_table[,i]
        
        # Save the name of the observed variable/OTU
        my_name <- colnames(input_table)[i]
        test_table <- xtabs(my_test_vector ~ independent_variable)
        
        # Number of groups to be compared
        num_of_represented_groups <- sum(as.vector(test_table) > 0)
        
        # Test whether any group is missing
        num_of_missing_groups <- nlevels(independent_variable) - num_of_represented_groups
        
        # Function to return "TRUE" or "FALSE" corresponding to "not missing" and "missing"
        if (num_of_missing_groups > 0) {
          isgroupmissing <- TRUE
        } else {
          isgroupmissing <- FALSE
        }
        
        # Performs a Kruskal-Wallis rank sum test
        fit <- tryCatch (kruskal.test(my_test_vector ~ independent_variable),error = function(i) {fail <<- TRUE})
        
        # Function to assign corrected and not-corrected pvalues  
        if (fail) {
          my_pvalue <- NaN
          #my_corrected_pvalue <- 0
          fail <- FALSE
        } else {
          
          # Round the p-value down to four decimals 
          my_pvalue <- round(fit$p.value,8)
        }
        
        # Add the p-values to the table
        newRow <- data.frame(name = my_name,pvalue = my_pvalue,missing = isgroupmissing)
        df <- rbind(df,newRow)
      }
      
      # Applying Benjamini-Hochberg (1995) correction 
      df$corrected <- round(p.adjust(df$pvalue, method = "BH"),8)
    }
    
    #######################   Paired-wilcoxon Test, (paired) Fisher's Exact Test    ########
    
    # Wilcoxon signed-rank test for pairwise comparisons
    count <- 1 
    x <- 0 
    
    # Vector with all possible group combinations
    idx <- combn(nlevels(independent_variable), 2)
    
    for (i in dependant_variables_start:dim(input_table)[2])
    {
      flag=TRUE
      
      # The vector of a dependant variable/OTU
      my_test_vector <- input_table[,i]
      
      # The name in the header of the dependant variable/OTU
      my_name <- colnames(input_table)[i]
      
      # Save the p-value and the corrected p-value
      if(isTRUE(KW)) pvalue <- df[count,2] ; cpvalue <- df[count,4]
      
      # Save information about missing group
      if(isTRUE(KW)) {
        missing_group <- df[count,3]
        }else{
          test_table <- xtabs(my_test_vector[!is.na(my_test_vector)] ~ independent_variable[!is.na(my_test_vector)])
          
          # Number of groups to be compared
          num_of_represented_groups <- sum(as.vector(test_table) > 0)
          
          # Test whether any group is missing
          num_of_missing_groups <- nlevels(independent_variable) - num_of_represented_groups
          
          # Function to return "TRUE" or "FALSE" corresponding to "not missing" and "missing"
          if (num_of_missing_groups > 0) {
            missing_group <- TRUE
          } else {
            missing_group <- FALSE
          }
        }
      
      count <- count + 1
      signif_pairs = data.frame( measure = as.character(),name = as.character(),Group1 = as.character(),Group2 = as.character(),pvalue = as.numeric())
      
      # Get the names of all group combinations
      idx_name <- combn(levels(independent_variable), 2)
      
      if (isFALSE(KW) || (isTRUE(KW) && !is.na(pvalue) && pvalue <= sig.cutoff)) {
        
        # Compute p-values from Wilcoxon test for all comparisons
        ppval_res <- numeric(ncol(idx))
        
        # Create an empty dataframe to hold the results of pairwise comparison 
        pair_pval_table <-data.frame(measure = character(0),pair = character(0),Group1 = character(0),Group2 = character(0),pvalue = numeric(0),corrected = numeric(0))
        if (isTRUE(KW)) pair_pval_table <- rbind(pair_pval_table,data.frame( measure = my_name,name = "All",Group1 = " -",Group2 = "- ",pvalue = pvalue, corrected = cpvalue))
        
        # Compute p-values of Wilcoxon test for all comparisons
        for (i in 1:ncol(idx)){
          # Performs a Wilcoxon rank sum test
          if(isTRUE(paired)){
            fit <- tryCatch (wilcox.test(my_test_vector[as.numeric(independent_variable) == idx[1,i]],my_test_vector[as.numeric(independent_variable) == idx[2,i]],paired=TRUE),error = function(i) {fail <<- TRUE})
          }else if(isFALSE(paired)){
            fit <- tryCatch (wilcox.test(my_test_vector[as.numeric(independent_variable) == idx[1,i]],my_test_vector[as.numeric(independent_variable) == idx[2,i]]),error = function(i) {fail <<- TRUE})
          }
          
          # Function to assign corrected and not-corrected pvalues  
          if (fail) {
            ppval_res[i] <- NaN
            fail <- FALSE
          } else {
            # Round the p-value down to four decimals 
            ppval_res[i] <- round(fit$p.value,8)
          }
          
          # Set the values of the pair and the corresponding p-value
          pair_name <- paste (idx_name[1,i],"-", idx_name[2,i],  sep = "")
          pair_num <- paste (idx[1,i],"-", idx[2,i],  sep = "")
          ppval <- round(ppval_res[i],8)
          
          # Create and add a new column to the plot table and the overall pairwise comparison table
          newRow <- data.frame( measure = my_name,name = pair_num,Group1 = idx_name[1,i],Group2 = idx_name[2,i],pvalue = ppval, corrected=0)
          pair_pval_table <- rbind(pair_pval_table,newRow)
        } 
        
        # Add the corrected p-values column to the dataframe
        if (isTRUE(KW)) pair_pval_table$corrected[-1] <- round(p.adjust(pair_pval_table$pvalue[-1], method = "BH"),8)
        Pforplot_table <- pair_pval_table[,c(-3,-4)]
        
        # Add the table with the corrected p-values to the complete list of pairwise p-values for all tests
        all_pair_pval_table <- rbind(all_pair_pval_table,pair_pval_table) 
        
        # Determine which groups are significantly different based on the results of the Wilcoxon test
        signif_pairs <- Pforplot_table[(Pforplot_table$pvalue < sig.cutoff) & !(is.na(Pforplot_table$pvalue)),]
        
      }
     

      # For the Fisher test:
      # If zeros are replaced by missing values, the prevalence and the points to be plotted are the same (number of samples with a value >0)
      # If zeros are not replaced and considered as true values, the prevalence is as above, but zeros are going to be plotted in the graphs
      if (ReplaceZero == "YES") {
        plot_df <- cbind.data.frame(abundance = my_test_vector,variable = independent_variable)
        plot_df$samplekaname <- row.names(input_table)
        plot_df_prevalence <- plot_df
      }else {
        plot_df<- cbind.data.frame(abundance = my_test_vector,variable = independent_variable)
        plot_df$samplekaname <- row.names(input_table)
        my_test_vector <- fill_zero.NA(my_test_vector,ReplaceZero = "YES")
        plot_df_prevalence <- cbind.data.frame(abundance = my_test_vector,variable = independent_variable)
        plot_df_prevalence$samplekaname <- row.names(input_table)
      }
      
      # Calculate prevalence by counting presence or absence of the given variable in each group
      prevalence <- table(plot_df_prevalence[!is.na(plot_df_prevalence[,1]),2])
      
      # Count how many samples are absent from total number of counts
      not_found <- total - prevalence 
      pre_table <- cbind(prevalence,not_found,total)
      
      # Calculate a two-sided Fisher's test 
      Fishtest <- fisher.test(pre_table[,-3],alternative = "two.sided",workspace=2e8)
      
      # Save the p-value of the Fisher's test
      fish_pvalue <- round(Fishtest$p.value,8)
      FnewRow <- data.frame(name = my_name,pvalue = fish_pvalue)
      Fdf <- rbind(Fdf,FnewRow)
      fppval_res <- numeric(ncol(idx))
      
      # Create an empty dataframe to hold the results of pairwise comparison for pairwise Fisher's Test
      pair_fpval_table <- data.frame(measure = character(0),name = character(0),Group1 = character(0),Group2 = character(0),pvalue = numeric(0),corrected = numeric(0))
      
      # Compute p-values from Wilcoxon test for all comparisons
      for (i in 1:ncol(idx))
      {
        pret <- as.data.frame(pre_table)
        pretrowbind <-rbind(pret[idx_name[1,i],-3],pret[idx_name[2,i],-3])
        
        # Compute two-sided Fisher's test
        fppval_res[i] <- fisher.test(pretrowbind,alternative = "two.sided",workspace=2e8)$p.value
        
        # Set values of the pair and corresponding p-value
        # Take variable name as "A-B"
        pair_name <- paste (idx_name[1,i],"-", idx_name[2,i],  sep = "")
        
        # Take variable name as "1-2"
        pair_num <- paste (idx[1,i],"-", idx[2,i],  sep = "")
        
        # Round the p-value to four decimals
        fppval <- round(fppval_res[i],8)
        
        # Create and add a new line to the plot table and the overall pairwise comparison table
        newRow <- data.frame(measure = my_name,name = pair_num,Group1 = idx_name[1,i],Group2 = idx_name[2,i],pvalue = fppval)
        pair_fpval_table <- rbind(pair_fpval_table,newRow)
      } 
      
      # Applying Benjamini-Hochberg (1995) correction 
      # Add a column with the corrected p-values 
      pair_fpval_table$corrected <- round(p.adjust(pair_fpval_table$pvalue, method = "BH"),8)
      
      # Add the table with corrected p-values to the complete list of pairwise p-values for all tests
      all_pair_fpval_table <- rbind(all_pair_fpval_table,pair_fpval_table) 
      
      # Make an object of "measure", "name", "p-value", "corrected" to print in PDF
      Fforplot_table <- pair_fpval_table[,c(-3,-4)]
      
      # Test whether the test is significant or not
      # If significant, then get value in "signif_fpairs"
      signif_fpairs <-Fforplot_table[(Fforplot_table$pvalue <= sig.cutoff) & !(is.na(Fforplot_table$pvalue)),]
      
      # Test for the case that there are no significant pairs
      x = x + 1
      # Check whether Fisher test is significant or not
      if (fish_pvalue <= sig.cutoff) {
        if (!is.na(signif_fpairs[1,1])) {
          signif_fpairs$measure <- as.character(signif_fpairs$measure)
          signif_fpairs$name <- as.character(signif_fpairs$name)
          colnames(signif_fpairs) <-c("Species","Groups","p-value","Adj. p-value")
          fpvaltable[[x]] <- list()
          
          # Pavlue table for significant pairs
          fpvaltable[[x]] <- tableGrob(signif_fpairs[2:4],rows = NULL,theme = mytheme)
          title <- textGrob("Fisher's Exact Test - pairwise",gp = gpar(fontsize = 9))
          padding <- unit(3,"mm")
          fpvaltable[[x]] <- gtable_add_rows(fpvaltable[[x]], heights = grobHeight(title) + padding,pos = 0)
          fpvaltable[[x]] <- gtable_add_grob(fpvaltable[[x]], title, 1, 1, 1, ncol(fpvaltable[[x]]))
        }else {
          FnewRow$name <- as.character(FnewRow$name)
          FnewRow$pvalue <- as.character(FnewRow$pvalue)
          FnewRow <- cbind(FnewRow$name,"-",FnewRow$pvalue,0)
          colnames(FnewRow) <-c("Species","Groups","p-value","Adj. p-value")
          allfpvaltable[[x]] <- list()
          allfpvaltable[[x]] <- tableGrob(FnewRow,rows = NULL,theme = mytheme)
        }
      }
      # Check whether Kruskal-Wallis test is significant or not
      if (!is.na(signif_pairs[1,1])) {
        signif_pairs$measure <- as.character(signif_pairs$measure)
        signif_pairs$name <- as.character(signif_pairs$name)
        colnames(signif_pairs) <-c("Species","Groups","p-value","Adj. p-value")
        pvaltable[[x]] <- list()
        pvaltableAll[[x]] <- list()
        
        # Pvalue table for significant pairs
        signif_all <- signif_pairs[1,c(1,3,4)]
        #if(dim(signif_pairs)[1] > 1) {
        signif_pairs <- signif_pairs[2:dim(signif_pairs)[1],]
        pvaltable[[x]] <-tableGrob(signif_pairs[2:4],rows = NULL,theme = mytheme)
        
        # Title of tables in the PDF
        title <- textGrob("Wilcoxon Rank Sum Test - pairwise",gp = gpar(fontsize = 9))
        padding <- unit(2,"mm")
        pvaltable[[x]] <- gtable_add_rows(pvaltable[[x]], heights = grobHeight(title) + padding,pos = 0)
        pvaltable[[x]] <- gtable_add_grob(pvaltable[[x]], title, 1, 1, 1, ncol(pvaltable[[x]]))
        colnames(signif_all) <-c(" ","p-value","Adj. p-value")
        pvaltableAll[[x]] <- tableGrob(signif_all,rows = NULL,theme = mytheme)
        title <- textGrob("Kruskal-Wallis Rank Sum Test - all groups ",gp = gpar(fontsize = 9))
        padding <- unit(2,"mm")
        pvaltableAll[[x]] <- gtable_add_rows(pvaltableAll[[x]], heights = grobHeight(title) + padding,pos = 0)
        pvaltableAll[[x]] <- gtable_add_grob(pvaltableAll[[x]], title, 1, 1, 1, ncol(pvaltableAll[[x]]))
      }
        }
    
    
    # Apply Benjamini-Hochberg (1995) correction to all p-values if only Wilcoxon and no Kruskall-Wallis was run (KW have been already corrected)
    # Add a column with the corrected p-values 
    if(isFALSE(KW)) all_pair_pval_table$corrected <- round(p.adjust(all_pair_pval_table$pvalue, method = "BH"),8)
    sig_all_pair_pval_table <- subset(all_pair_pval_table,all_pair_pval_table$pvalue<=sig.cutoff)
    Fdf$corrected <- round(p.adjust(Fdf$pvalue, method = "BH"),8)
    sig_Fdf <- subset(Fdf,Fdf$pvalue<=sig.cutoff)
    counter = 1
    
    ################################################################################
    ######                        Write Output Files in Separate folder       ######
    #################################################################################
    
    # Take the name of the independent variable to name the folder
    prefix = paste(params[[q]]$name, taxon, sep="_")
    
    # Make a directory name with independent variable name and date
    if(isTRUE(paired)){ dir = "paired" }else{ dir = "unpaired" }
    newdir <- file.path(rhea.dir, dir, paste(prefix, sep = "_"))
    
    # Create a directory 
    dir.create(newdir, recursive = T)
    
    # Main file generated after all pre-processing and used for statistical analysis
    # Serial-Group-Comparisons input Table
    write.table(input_table, file.path(newdir, "modified.txt"), sep = "\t", col.names = NA, quote = FALSE)
    
    # Table with the results of Kruskal-Wallis test 
    if(isTRUE(KW)) write.table(df[,c(1,2,4)], file.path(newdir, "pvalues.tab"), sep = "\t", col.names = NA, quote = FALSE)
    
    # Table with the results of Wilcoxon Rank Sum Test

    # Table with the results of Wilcoxon Rank Sum Test
    if(taxon == "ASV"){
      tax_tab <- read.xlsx(file.path(res.dir, "2.blast/tax_tab_SILVA_blast_species.xlsx"), rowNames = T)
      rownames(tax_tab) <- gsub("_", "", rownames(tax_tab))
      all_pair_pval_table <- merge(all_pair_pval_table, tax_tab, by.x = "measure", by.y = 0, all.x = T, all.y = F)
    }
    write.table(all_pair_pval_table[,-2], file.path(newdir, "sign.tab"), sep = "\t", col.names = NA, quote = FALSE)
    write.xlsx(all_pair_pval_table[,-2], file.path(newdir, "sign_pairs.xlsx"))
    
    # Table with the results of Fisher's Exact Test
    write.table(Fdf, file.path(newdir, "FisherTestAll.tab"), sep = "\t", col.names = NA, quote = FALSE)
    
    # Table with the results of pairwise Fisher's Exact Test
    write.table(all_pair_fpval_table[,-2], file.path(newdir, "FisherTestPairWise.tab"), sep = "\t", col.names = NA, quote = FALSE)
    
    plot_abundance = function(data, ylabn = "Relative abundance (%)",
                              Facet = "Taxa",
                              Color = "Diagnosis",
                              palette,
                             # n = NULL,
                              legend_title = NULL){
      
      # Select significant taxa to plot
      all_pair_pval_table_ordered = all_pair_pval_table[order(all_pair_pval_table$pvalue),]
      significant_pval_table_ordered = all_pair_pval_table_ordered[all_pair_pval_table_ordered$corrected < sig.cutoff &
                                  !is.na(all_pair_pval_table_ordered$corrected),]
      
      if(isTRUE(KW)){
        sig = significant_pval_table_ordered[significant_pval_table_ordered$Group1 == " -", "measure"] 
      }else{sig = significant_pval_table_ordered[, "measure"] }
      
      # if(is.null(n)){
      #   if(length(sig) <= 5){
      #     n = 1
      #   }else if(length(sig) <= 10){
      #     n = 2
      #   }else if(length(sig) <= 15){
      #     n = 3
      #   }else if(length(sig) <= 20){
      #     n = 4
      #   }
      # }
      nfeatures <- length(sig)
      if(nfeatures > 0){
        n = ceiling(nfeatures/5)
        ncol = min(5, nfeatures)
        
        if(is.null(palette)){
          palette <- scales::hue_pal()(nrow(unique(sample_data(physeq)[,condition_var])))
        }
        
        # Select metadata columns
        fixed_columns <- names(data)[1:(taxonomic_variables_start-1)]
        
        # Keep only metadata columns and significant taxa
        data <- data[,c(fixed_columns, sig)]
        
        # Transform data frame into long format for plotting
        data_long <- data %>%
          pivot_longer(
            cols = -all_of(fixed_columns),  # Select all columns EXCEPT the fixed ones
            names_to = "Taxa",              # New column for bacterial taxonomic names
            values_to = "Abundance"         # New column for abundance values
          )
        
        # Function for transforming y axis scale
        marks_no_sci <- function(x) format(x, big.mark = ".", decimal.mark = ",", scientific = F)
        
        P <- ggplot(data = data_long,
                    mapping = aes_string(x = Color, y = "Abundance",
                                         color = Color, fill = Color)) + thm +
          geom_boxplot(fill = NA) +
          facet_wrap(facets = Facet, nrow = n, ncol = ncol) +
          ylab(ylabn) +
          # stat_compare_means(method = "wilcox") +
          scale_fill_manual(values = palette) + scale_color_manual(values = palette) +
          #  scale_y_log10()
          scale_y_log10(labels = marks_no_sci)
        
        if (!is.null(legend_title)) {
          P <- P + labs(color = legend_title, fill = legend_title)
        } 
        print(P)
        
        if(nfeatures == 1){w = 4}else if(nfeatures < 5){w = 3.2*nfeatures}else{w = 16}
        ggsave(file.path(newdir, paste(prefix,"-sign.svg", sep = "")), width = w, height = 6*n, dpi = 300, limitsize = FALSE)
        # ggsave(file.path(newdir, paste(prefix,"-sign.tiff", sep = "")), width = 16, height = 8, dpi = 300)
        
      }else{
        print(paste("No significant taxa to plot for ", q, " and ", taxon))
      }
    }
    
    plot_abundance(original_table, Color = independant_variable_name, palette = params[[q]]$col, legend_title = params[[q]]$lab_plot) 
    
    # Adding log file in analysis
    sink(file = file.path(newdir, "my_analysis_log.txt"))
    cat ("***************************************","\n")
    cat ("Parameters Used for Analysis","\n")
    cat ("***************************************","\n","\n")
    cat ("study.name:OTUsCombined","\n","\n")
    cat ("input_file:",input_file,"\n","\n")
    cat ("independant_variable_name:",independant_variable_name,"\n","\n")
    cat ("dependant_variables_start:",dependant_variables_start,"\n","\n")
    cat ("taxonomic_variables_start:",taxonomic_variables_start,"\n","\n")
    cat ("perform_Kruskall-Wallis/Friedman test:",KW,"\n","\n")
    cat ("paired",paired,"\n","\n")
    cat ("abundance_cutoff:", abundance_cutoff,"\n","\n")
    cat ("prevalence_cutoff:",prevalence_cutoff,"\n","\n")
    cat ("perform_max_median_cutoff:",perform_max_median_cutoff,"\n","\n")
    cat ("max_median_cutoff:",max_median_cutoff,"\n","\n")
    cat ("PlotOption:",PlotOption,"\n","\n")
    cat ("ReplaceZero:",ReplaceZero,"\n","\n")
    cat ("sig.cutoff:",sig.cutoff,"\n","\n")
    sink()

    if(!flag) { stop("
                 It was not possible to install all required R libraries properly.
                 Please check the installation of all required libraries manually.\n
                 Required libaries:plotrix,PerformanceAnalytics,reshape,ggplot2,gridExtra,grid,ggrepel,gtable,Matrix,cowplot")
    }
    
    ########################################################
    ##    Script Ended !!!!
    #########################################################
  
  
  
  
  

