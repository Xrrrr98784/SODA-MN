#remove what has been found in previous round
remove.previous.found <- function(pre,new){
  i=1
  while (i < nrow(pre)+1) {
    j=1
    while (j < nrow(new)+1) {
      if (new[j,3] == pre[i,2] || (new [j,2]==pre[i,2] && new [j,3]==pre[i,3])){
        new=new[-c(j),]
        j=j-1
      }
      j=j+1
    }
    i=i+1
  }
  return(new)
}
