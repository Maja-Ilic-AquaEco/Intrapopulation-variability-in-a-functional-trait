ggHistNorm <- function(data, x, variable.lab, binwidth) {
  
  norm.p <- shapiro.test(x)$p.value
  
  if(norm.p < 0.001){
    p.lab <- "p < 0.001"
  }
  
  if(norm.p >= 0.001 & norm.p < 0.01){
    p.lab <- "p < 0.01"
  }
  
  if(norm.p >= 0.01 & norm.p < 0.05){
    p.lab <- "p < 0.05"
  }
  
  if(norm.p >= 0.05 ){
    p.lab <- parse(text = paste0('p = ', round(norm.p, digits = 3)))
  }
  
  Median <- median(x, na.rm = T)
  Mean <- mean(x, na.rm = T)
  SD <- sd(x, na.rm = T)
  Var <- var(x, na.rm = T)
  
  data_Mean_SD <- length(x[x >= Mean - SD & x <= Mean + SD]) / length(!is.na(x)) * 100
  
  plot.norm <- ggplot(data, aes(x = x)) +
    geom_histogram(aes(y = ..density..), binwidth = binwidth, alpha = 0.6, fill = "lightblue", color = "grey70") + 
    geom_density(alpha = 0.5, fill = "lightblue", color = "grey40") +
    theme_minimal() + 
    stat_function(fun = dnorm, color = rgb(22, 160, 133, max = 255), size = 1, 
                  args = list(mean = Mean, sd = SD)) +
    geom_vline(xintercept = Mean, color = "blue", size = 1) +
    geom_vline(xintercept = Mean*0.9, color = "blue", linetype = "dashed", size = 0.6) +
    geom_vline(xintercept = Mean*1.1, color = "blue", linetype = "dashed", size = 0.6) +
    geom_vline(xintercept = Median, color = "red", size = 1) +
    geom_vline(xintercept = Median*0.9, color = "red", linetype = "dashed", size = 0.6) +
    geom_vline(xintercept = Median*1.1, color = "red", linetype = "dashed", size = 0.6) +
    geom_vline(xintercept = Mean - SD, color = "green", linetype = "dashed", size = 0.8) +
    geom_vline(xintercept = Mean + SD, color = "green", linetype = "dashed", size = 0.8) +
    labs(x = variable.lab, 
         y = "Density",
         title = paste0("Normality of ", variable.lab),
         subtitle = paste0("Shapiro's-Test: ",p.lab,
                           "\nMedian (red line): ",round(Median, digits = 4),
                           "\nMean (blue line): ",round(Mean, digits = 4),
                           "\nVariance: ",round(Var, digits = 4),
                           "\nSD (green line): ",round(SD, digits = 4),
                           "\nData within Mean \u00B1 1 SD: ",round(data_Mean_SD, digits = 2)," %"))
  
  return(list(plot.norm, norm.p))
}
