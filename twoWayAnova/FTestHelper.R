FTest <- function(FStat, df1, df2, alpha) {
  FQuant <- qf((1-alpha), df1, df2);
  pVal <- (1 - pf(FStat, df1, df2));
  
  if(FStat > FQuant)  {
    result <- paste("Significant. P Value: ", toString(pVal), sep = "");
  } else  {
    result <- paste("not significant. P Value: ", toString(pVal), sep = "");
  }
  
  return(result);
}
