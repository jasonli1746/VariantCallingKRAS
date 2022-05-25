#function to return actual codon based on variant call results
#Parameters:
# startposition = int, starting position on chromosome
# refcodon = string of 3 characters, e.g. "GCC"
# kras.mut.df = result from calcVAFGeneMutList function
#Returns:
# codon = string of 3 chars, alternative codon based on mutation calls, e.g. "GGC"

checkMutation <- function(startposition, refcodon, kras.mut.df){
  firstpos <- kras.mut.df %>% filter(grepl(startposition+2,position)) %>% dplyr::select(nuc_call) 
  if(nrow(firstpos)>0){
    firstpos.res <- stringr::str_sub(firstpos, 5, 5)
  }else{
    firstpos.res <- stringr::str_sub(refcodon, 1, 1)
  }
  secondpos <- kras.mut.df %>% filter(grepl(startposition+1,position)) %>% dplyr::select(nuc_call) 
  if(nrow(secondpos)>0){
    secondpos.res <- stringr::str_sub(secondpos, 5, 5)
  }else{
    secondpos.res <- stringr::str_sub(refcodon, 2,2)
  }
  thirdpos <- kras.mut.df %>% filter(grepl(startposition,position)) %>% dplyr::select(nuc_call) 
  if(nrow(thirdpos)>0){
    thirdpos.res <- stringr::str_sub(thirdpos, 5, 5)
  }else{
    thirdpos.res <- stringr::str_sub(refcodon, 3,3)
  }
  codon <- paste0(firstpos.res, secondpos.res, thirdpos.res)
  return(codon)
}