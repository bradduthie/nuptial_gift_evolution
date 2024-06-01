female_IBM_criteria <- function(MGf, MLf, lambda = 1, TF = 2){
  MM   <- MGf + MLf;
  gamF <- lambda / (TF * (MM - MLf));
  return(gamF);
}

male_IBM_criteria <- function(Male_M, alpha = 1, lambda = 1){
  return(lambda * alpha * Male_M);
}

simpleboot <- function(freqs, repli = 1000, alpha = 0.05){  
  vals  <- NULL;                             
  i     <- 0;                                    
  while(i < repli){                             
    boot  <- sample(x = freqs, size = length(freqs), replace = TRUE);     
    strap <- mean(boot);                      
    vals  <- c(vals,strap);                     
    i     <- i + 1;                          
  }                                           
  vals   <- sort(x = vals, decreasing = FALSE);    
  lowCI  <- vals[round( (alpha * 0.5) * repli )];    
  highCI <- vals[round( (1-(alpha*0.5) ) * repli)];     
  CIs    <- c(lowCI, highCI);        
  return(CIs);         
} 

IBM_output_to_summary <- function(res_file, out = "summary.csv"){
  dat  <- read.csv(file = res_file);
  # Remove any titles if CSVs are concatenated
  mims <- which(dat[,1] == "mim");
  dat  <- dat[-mims,];
  cols <- dim(dat)[2];
  for(i in 1:cols){
    dat[,i] <- as.numeric(dat[,i]);
  }
  
  # Pull out summary statistics
  res  <- unique(cbind(dat$a1, dat$gam));
  SUMM <- NULL;
  for(i in 1:dim(res)[1]){
    a1_val  <- res[i, 1];
    gam_val <- res[i, 2];
    temp_r  <- which(dat$a1 == a1_val & dat$gam == gam_val);
    mean_Tm <- mean(dat[temp_r,]$Tm);
    CIs_Tm  <- simpleboot(freqs = dat[temp_r,]$Tm);
    mean_Rj <- mean(dat[temp_r,]$RejPr);
    CIs_Rj  <- simpleboot(freqs = dat[temp_r,]$RejPr);
    mean_Gf <- mean(dat[temp_r,]$Gift);
    mean_Nu <- mean(dat[temp_r,]$N_m);
    CIs_Nu  <- simpleboot(freqs = dat[temp_r,]$N_m);
    mean_M  <- mean(dat[temp_r,]$M_m);
    mean_B  <- mean(dat[temp_r,]$Beta);
    mean_Mm <- mean(dat[temp_r,]$M_males);
    mean_Fm <- mean(dat[temp_r,]$M_females);
    mean_MG <- mean(dat[temp_r,]$F_Mnupt_enc);
    mean_MN <- mean(dat[temp_r,]$F_Mnonupt_enc);
    reps    <- length(temp_r);
    new_row <- c(a1_val, gam_val, mean_Tm, CIs_Tm, mean_Rj, CIs_Rj, mean_Gf,
                 mean_Nu, CIs_Nu, mean_M, mean_B, mean_Mm, mean_Fm, mean_MG,
                 mean_MN, reps);
    SUMM    <- rbind(SUMM, new_row);
    print(i);
  }
  
  colnames(SUMM) <- c("a1_val", "gam_val", "mean_TM", "LCI_Tm", "UCI_Tm",
                      "mean_Rj", "LCI_Rj", "UCI_Rj", "Gift", "Neutral",
                      "LCI_Neutral", "UCI_Neutral", "M", "Beta", "M_m", "F_m",
                      "Mf_nupt", "Mf_no_nput", "Replicates");
  
  write.csv(x = SUMM, file = out, row.names = FALSE);
}



plot_SI_mortality <- function(SUMM, outfile = "mortality_sensitivity.eps"){
  alpha_vals   <- sort(unique(SUMM[,3]));
  mort_combn   <- paste(SUMM[,1], SUMM[,2]);
  mort_label   <- unique(paste(SUMM[,1], ", ", SUMM[,2], sep = ""));
  TM_vals      <- tapply(X = SUMM[,4],  INDEX = mort_combn, FUN = mean);
  m_vals       <- tapply(X = SUMM[,14], INDEX = mort_combn, FUN = mean);
  beta_vals    <- tapply(X = SUMM[,15], INDEX = mort_combn, FUN = mean);
  Male_M       <- tapply(X = SUMM[,16], INDEX = mort_combn, FUN = mean);
  Female_M     <- tapply(X = SUMM[,17], INDEX = mort_combn, FUN = mean);
  Mf_nupt      <- tapply(X = SUMM[,18], INDEX = mort_combn, FUN = mean);
  Mf_nonupt    <- tapply(X = SUMM[,19], INDEX = mort_combn, FUN = mean);
  male_search  <- as.numeric(which(SUMM[,5] > 0 & SUMM[,6] > 0));
  fem_choose   <- as.numeric(which(SUMM[,8] > 0 & SUMM[,9] > 0));
  mimmom       <- unique(cbind(dat$mim, dat$mom));
  
  mm <- SUMM[male_search,]
  mp <- NULL;
  for(i in 1:dim(mm)[1]){
    rval <- which(mm[i, 1] == mimmom[,1] & mm[i, 2] == mimmom[,2]);
    mp   <- rbind(mp, c(rval, mm[i, 3]));
  }

  ff <- SUMM[fem_choose,]
  fp <- NULL;
  for(i in 1:dim(ff)[1]){
    rval <- which(ff[i, 1] == mimmom[,1] & ff[i, 2] == mimmom[,2]);
    fp   <- rbind(fp, c(rval, ff[i, 3]));
  }
  
  setEPS();
  postscript(outfile);
  par(mar = c(8, 5, 1, 1));
  plot(x = mp[,1], y = mp[,2], pch = 20, ylim = c(0, 2), col = "blue",
       xlab = "", ylab = "", xaxt = "n", cex.lab = 1.25);
  points(x = fp[,1], y = fp[,2], cex = 2.5, col = "red", pch = 20);
  points(x = mp[,1], y = mp[,2], cex = 1.5, col = "blue", pch = 20);
  axis(side = 1, labels = mort_label, at = 1:25, las = 2);
  mtext(text = expression(paste("Mortality combination (",mu["in"],", ", 
                                mu["out"],")")), side = 1, line = 6.5, 
        cex = 1.25); 
  mtext(text = expression(paste("Fitness threshold (", gamma[f], " & ", 
                                gamma[m], ")")), side = 2, line = 2.5, 
        cex = 1.25); 
  dev.off();
}



