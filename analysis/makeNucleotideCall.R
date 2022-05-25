#function to call base of mutation
#Parameters:
# n = int from 1 to 3 (length of a codon)
# df = data frame, kras.mut.df
#Returns:
# var.return = string corresponding to called base, either "fracA", 
#  "fracC", "fracG", or "fracT"
makeNucleotideCall <- function(n, df){
  row.try <- df[n,3:6]
  second.highest <- sort(row.try)[length(row.try)-1]
  second.allele <- substr(colnames(second.highest), 5,5)
  if(second.highest>0.15 & second.highest<0.5){
    var.return <- colnames(second.highest)
  } else if(second.highest==0.5 & second.allele==df[n,7]){
    first.highest <- sort(row.try)[length(row.try)]
    var.return <- colnames(first.highest)
  }
  else{
    var.return <- paste0("frac",df[n,7])
  }
  return(var.return)
}