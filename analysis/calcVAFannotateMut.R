#function to calculate Variant Allele Frequency and annotate a filtered VCF
#Parameters:
# n = int, number from 1 to length of res.list
# orig.AA = char, original amino acid, e.g. "G"
# res.list = list of dataframes, result of kras-allelesfromprot-aalocs.sh
#Returns:
# kras.mut.df = dataframe, annotated VCF with columns chr, startpos, ref_allele, alt_allele,
#  vaf_rnaseq, nuc_call, codon, mutation, bam_idx
#  chr, startpos, ref_allele, alt_allele are 4 columns of VCF
#  vaf_rnaseq is calculated from VCF, frac
#  nuc_call is called base
#  codon is called codon
#  mutation is HGVSp_short mutation
#  bam_idx is index of input BAM file
calcVAFannotateMut <- function(n, orig.AA, res.list){
  
  #orig.AA = aminoAcidTable(refcodon)
  aa.row <- amino.acid.table[which(amino.acid.table$orig_AA==orig.AA),]
  ref.AA <- aa.row$AA.1
  refcodon=aa.row$codon
  startpos=aa.row$genomic_position
  kras.mut.df <- res.list[[n]]
  row.sel <- names(res.list)[[n]]
  kras.mut.df <- kras.mut.df %>% mutate(vaf_rnaseq=case_when(refAllele=="A" ~ 1-fracA,
                                                             refAllele=="C" ~ 1-fracC,
                                                             refAllele=="G" ~ 1-fracG,
                                                             refAllele=="T" ~ 1-fracT))
  kras.mut.df$nuc_call <- names(kras.mut.df)[3:6][max.col(kras.mut.df[3:6], "first")]
  codon <- checkMutation(startpos, refcodon, kras.mut.df)
  
  kras.mut.df$codon <- codon
  kras.mut.df <- kras.mut.df %>% mutate(AA=aminoAcidTable(codon))
  
  kras.mut.df <- kras.mut.df %>% mutate(mutation=ifelse(AA!=ref.AA, paste0("p.", orig.AA, AA), "no.chg"))
  kras.mut.df$bam_idx <- row.sel
  return(kras.mut.df)
}