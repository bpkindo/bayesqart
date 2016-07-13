## weighted mean absoluted deviation
#'A function to calculate weigthed mean absolute deviation.  
#'@param y_true actual response,
#'@param y_hat conditional quantile prediction,
#'@param quant quantile value.
#'@return weighted mean absoulute deviation, \code{wmad}.
#'@useDynLib BayesQArt
wmad = function(y_true,y_hat,quant){
  n=length(y_true);
  wmad = sum(quant*abs(y_true - y_hat)*(y_true > y_hat) + (1-quant)*abs(y_true - y_hat)*(y_true <= y_hat))/n;
  return(wmad);
}
