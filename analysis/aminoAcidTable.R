#function to return amino acid code based on codon
#Parameters:
# codon = 3 char string, e.g. "CCA"
#Returns:
# AA = 1 char, e.g. "G"

aminoAcidTable <- function(codon){
  AA <- case_when(codon %in% c("CCA", "CCC", "CCG", "CCT") ~ "G",
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
                  codon %in% c("ATT", "ATC", "ACT") ~ "stop")
  return(AA)
}