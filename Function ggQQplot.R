ggQQplot <- function(data, x, variable.lab) {
  
  plot.qq <- ggplot(data, aes(sample = x)) +
    stat_qq(color = rgb(44, 62, 80, max = 255), size = 2) + 
    stat_qq_line(size = 1, color = rgb(22, 160, 133, max = 255)) +
    theme_minimal() +
    labs(x = "Theoretical Quantiles", 
         y = "Sample Quantiles", 
         title = "Normal Q-Q Plot",
         subtitle = paste0("Variable: ", variable.lab))
  
  return(plot.qq)
  
}
