#remove redundancy by nodes

remove.node.duplicates <- function (x,prior) {
  x <- x[order(x[,2]), ]
  x <- x[order(x[,4]), ]
  i <- 1
  if (prior=="match"){
  while (i < nrow(x)) { 
    if ( x[i,2]==x[i+1,2] &&x[i,4]==x[i+1,4] ){ 
      if (x$match[i] < x$match[i+1]) {
        x=x[-c(i),]
        i=i-1
      } else if (x$match[i] == x$match[i+1]) {
        if (x$intensity[i] < x$intensity[i+1]) {
          x=x[-c(i),]
          i=i-1
        }else{
          x=x[-c(i+1),]
          i=i-1 
        }
      }else{
        x=x[-c(i+1),]
        i=i-1
      } 
    }
    
    i=i+1
}
} else if (prior=="intensity"){
  while (i < nrow(x)) { 
    if ( x[i,2]==x[i+1,2] &&x[i,4]==x[i+1,4] ){
      if (x$intensity[i] < x$intensity[i+1]) {
        x=x[-c(i),]
        i=i-1
    }else{
        x=x[-c(i+1),]
        i=i-1
      } 
    }
    i=i+1
}
}
  return(x)
}
