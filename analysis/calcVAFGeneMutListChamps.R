#function to calculate Variant Allele Frequency from selected rows of VCFs
#Parameters:
# n = integer from 1 to length of res.list
# aa.pos = amino acid position, amino acid of mutation (e.g. G12C, 12 or Q61H, 61)
# res.list = list of results from KRAS-allelesfromprot.sh, list of dataframes, each which is a
#  filtered rows from VCF
#Returns: 
# kras.mut.df = data frame with columns chr, startpos, ref_allele, alt_allele,
#  vaf_rnaseq, nuc_call, codon, mutation
#  chr, startpos, ref_allele, alt_allele are 4 columns of VCF
#  vaf_rnaseq is calculated from VCF, frac
#  nuc_call is called base
#  codon is called codon
#  mutation is HGVSp_short mutation

calcVAFGeneMutListChamps <- function(n, aa.pos, res.list){
  kras.mut.df <- res.list[[n]]
  row.sel <- names(res.list)[[n]]
  kras.mut.df <- kras.mut.df %>% mutate(vaf_rnaseq=case_when(refAllele=="A" ~ 1-fracA,
                                                             refAllele=="C" ~ 1-fracC,
                                                             refAllele=="G" ~ 1-fracG,
                                                             refAllele=="T" ~ 1-fracT))
  kras.mut.df$nuc_call <- names(kras.mut.df)[3:6][max.col(kras.mut.df[3:6], "first")]
  firstpos <- kras.mut.df %>% filter(grepl("25227340",position)) %>% dplyr::select(nuc_call) 
  if(nrow(firstpos)>0){
    firstpos.res <- stringr::str_sub(firstpos, 5, 5)
  }else{
    firstpos.res <- "C"
  }
  
  secondpos <- kras.mut.df %>% filter(grepl("25227339",position)) %>% dplyr::select(nuc_call) 
  if(nrow(secondpos)>0){
    secondpos.res <- stringr::str_sub(secondpos, 5, 5)
  }else{
    secondpos.res <- "T"
  }
  thirdpos <- kras.mut.df %>% filter(grepl("25227338",position)) %>% dplyr::select(nuc_call) 
  if(nrow(thirdpos)>0){
    thirdpos.res <- stringr::str_sub(thirdpos, 5, 5)
  }else{
    thirdpos.res <- "C"
  }
  codon <- paste0(firstpos.res, secondpos.res, thirdpos.res)
  kras.mut.df$codon <- codon
  kras.mut.df <- kras.mut.df %>% mutate(AA=case_when(codon %in% c("CCA", "CCC", "CCG", "CCT") ~ "G",
                                                     codon %in% c("CTA", "CTG") ~ "D",
                                                     codon %in% c("CAA", "CAC", "CAG", "CAT") ~ "V",
                                                     codon %in% c("ACA", "ACG") ~ "C",
                                                     codon %in% c("CGA", "CGC", "CGG", "CGT") ~ "A",
                                                     codon %in% c("TGA", "TGC", "TGG", "TGT") ~ "T",
                                                     codon %in% c("GTC", "GTT") ~ "Q",
                                                     codon %in% c("GTA", "GTG") ~ "H",
                                                     codon %in% c("GCA", "GCC", "GCG", "GCT", "TCT", "TCC") ~ "R",
                                                     codon %in% c("TCG", "TCA", "AGA", "AGC", "AGG", "AGT") ~ "S",
                                                     codon %in% c("GAA", "GAC", "GAG", "GAT", "AAT", "AAC") ~ "L",
                                                     codon %in% c("ATA", "ATG") ~ "Y",
                                                     codon %in% c("CTC", "CTT") ~ "E",
                                                     codon %in% c("TAC") ~ "M",
                                                     codon %in% c("TAA", "TAG", "TAT") ~ "I",
                                                     codon %in% c("AAA", "AAG") ~ "F",
                                                     codon %in% c("ACC") ~ "W",
                                                     codon %in% c("GGA", "GGC", "GGG", "GGT") ~ "P",
                                                     codon %in% c("TTG", "TTA") ~ "N",
                                                     codon %in% c("TTC", "TTT") ~ "K",
                                                     codon %in% c("ATT", "ATC", "ACT") ~ "stop"))
  orig.AA="E"
  kras.mut.df <- kras.mut.df %>% mutate(mutation=ifelse(AA!=orig.AA, paste0("p.", orig.AA, aa.pos, AA), "no.chg"))
  #G, D, V, C, R, S, A, I, L, M, Q, T, Y, E
  #G (Gly) CCA, CCC, CCG, CCT
  #D (Asp) CTA, CTG
  #V (Val) CAA, CAC, CAG, CAT
  #C (Cys) ACA, ACG
  #A (Ala) CGA, CGC, CGG, CGT
  #T (Thr) TGA, TGC, TGG, TGT
  #Q (Gln) GTC, GTT
  #H (His) GTA, GTG
  #R (Arg) GCA, GCC, GCG, GCt, TCT, TCC
  #S (Ser) TCG, TCA, AGA, AGC, AGG, AGT
  #L (Leu) GAA, GAC, GAG, GAT, AAT, AAC
  #Y (Tyr) ATA, ATG
  #E (Glu) CTC, CTT
  #M (Met) TAC
  #I (Ile) TAA, TAG, TAT
  #F (Phe) AAA, AAG
  #W (Trp) ACC
  #P (Pro) GGA, GGC, GGG, GGT
  #N (Asn) TTG, TTA 
  #K (Lys) TTC, TTT
  #stop ATT, ATC, ACT
  kras.mut.df$bam_idx <- row.sel
  return(kras.mut.df)
}