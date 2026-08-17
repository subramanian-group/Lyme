#Title: General Functions
#Author: Clifford Rostomily


# Color Gradient from Values ----------------------------------------------

val2gradient = function(x, gradient = c("blue", "red")){
  rgb(colorRamp(gradient)( (x + abs(min(x)))/(max(x + abs(min(x))) ))/255)
}


# Return Significance Stars -----------------------------------------------

returnSigStars = function(x){
  x2 = x
  x2[is.numeric(x) & x > 0.05] = "ns"
  x2[is.numeric(x) & x < 0.05] = "*"
  x2[is.numeric(x) & x < 0.01] = "**"
  x2[is.numeric(x) & x < 0.001] = "***"
  x2[is.numeric(x) & x < 0.0001] = "****"
  x2
}
