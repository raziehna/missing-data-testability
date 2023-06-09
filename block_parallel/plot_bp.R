

process_df <- function(df_1, df_2, name="binary_boxplots.png", title_name='Block-parallel model holds'){

  
  dat = data.frame(rbind(df_1[, 1], df_1[, 2], df_1[, 3], df_1[, 4], df_1[, 5], 
                         df_2[, 1], df_2[, 2], df_2[, 3], df_2[, 4], df_2[, 5]))
  
  dat$type = c("1k", "3k", "5k", "7k", "9k", 
               "1k", "3k", "5k", "7k", "9k")
  
  dat$case = c("bin", "bin", "bin", "bin", "bin",  
               "cont", "cont", "cont", "cont", "cont")
  
  stacked.dat = melt(dat, id = c('type', 'case'))
  stacked.dat = stacked.dat[, -3] 
  
  png(filename= paste0(name), 
      units="in", 
      width=6,
      height=8,
      pointsize=12,
      res=200)
  
  boxplots.triple = boxplot(value~type + case, 
                            data = stacked.dat, 
                            las = 2,
                            names = dat$type, 
                            at = c(1:5, 7:11), 
                            outline=FALSE, 
                            # xaxt='n',
                            # yaxt='n', 
                            width=rep(1, 10),
                            # ylim = c(0, 6.0),
                            ylim = c(0, 2.0),
                            col = c('gray', 'gray', 'gray', 'gray', 'gray'), 
                            xlab="", ylab="")
  
  axis(side=2, at=c(-0.75, -0.5)) 
  title(title_name)
  title(ylab = "Odds ratio vlaues", line = 2.5) 
  
  text(c(3, 9.5),
       # c(5.9, 5.9),
       c(1.9, 1.9),
       c('Binary Data', 'Gaussian Data'),
       font = 1)
  
  # dashed lines
  # y = seq(0.2, 6.0, by = 0.1)
  y = seq(0.1, 2, by = 0.1)
  x = rep(6, length(y))
  lines(x, y, type="l", lty=2) 
  
  # true value 
  x = seq(0, 12, by=1)
  y = rep(1, length(x))
  lines(x, y, type="l", lty=1, col="red")
  
  dev.off()

}

  
df_1 = read.csv(paste0("BP_binary.csv"), sep = ",")
df_2 = read.csv(paste0("BP_continuous.csv"), sep = ",")
df_3 = read.csv(paste0("NoSelf_binary.csv"), sep = ",")
df_4 = read.csv(paste0("NoSelf_continuous.csv"), sep = ",")


process_df(df_1, df_2, name="BP_boxplots.png", title_name='Block-parallel model holds')
  
# process_df(df_3, df_4, name="NoSelf_boxplots.png", title_name='Block-parallel model does not hold')

