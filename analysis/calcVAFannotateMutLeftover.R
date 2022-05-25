#same function as calcVAFannotateMut, but for leftover VCF results that
#have multiple mutation calls
#Parameters:
# n = int, number from 1 to length of res.list
# orig.AA = char, original amino acid, e.g. "G"
# res.list = list of dataframes, result of kras-allelesfromprot-aalocs.sh
#Returns:
# kras.mut.df = dataframe, annotated VCF with columns chr, startpos, ref_allele, alt_allele,
#  vaf_rnaseq, nuc_call, codon, mutation, bam_idx
#  chr, startpos, ref_allele, alt_allele are 4 columns of VCF
#  vaf_rnaseq is calculated from VCF, frac
#  nuc_call is called base (but now there can be multiple calls)
#  codon is called codon
#  mutation is HGVSp_short mutation
#  bam_idx is index of input BAM file
calcVAFannotateMutLeftover <- function(n, orig.AA, res.list){
  
  aa.row <- amino.acid.table[which(amino.acid.table$orig_AA==orig.AA),]
  ref.AA <- aa.row$AA.1
  refcodon=aa.row$codon
  startpos=aa.row$genomic_position
  #orig.AA = aminoAcidTable(refcodon)
  kras.mut.df <- res.list[[n]]
  row.sel <- names(res.list)[[n]]
  kras.mut.df <- kras.mut.df %>% mutate(vaf_rnaseq=case_when(refAllele=="A" ~ 1-fracA,
                                                             refAllele=="C" ~ 1-fracC,
                                                             refAllele=="G" ~ 1-fracG,
                                                             refAllele=="T" ~ 1-fracT))
  nuc_calls <- sapply(c(1:3), makeNucleotideCall, df=kras.mut.df)
  kras.mut.df$nuc_call <- nuc_calls
  codon <- checkMutation(startpos, refcodon, kras.mut.df)
  
  kras.mut.df$codon <- codon
  kras.mut.df <- kras.mut.df %>% mutate(AA=aminoAcidTable(codon))
  
  kras.mut.df <- kras.mut.df %>% mutate(mutation=ifelse(AA!=ref.AA, paste0("p.", orig.AA, AA), "no.chg"))
  kras.mut.df$bam_idx <- row.sel
  return(kras.mut.df)
}