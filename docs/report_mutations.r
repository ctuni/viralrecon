vcffile <- paste("mutations/", sample, ".snpSift.table.modified.txt", sep= "")
mutations_df <- read.table(vcffile, sep = "\t", header = T)
ref_alt_nuc <- mutations_df[,4]
mutations_db <- read.table(mutations_tsv, header = F, sep = "\t")
mutations_column <- as.character(mutations_db$V1)
mutations_df$EFF....AA = gsub("\\.", "",mutations_df$EFF....AA)
mutations_df$EFF....AA = gsub(",", "",mutations_df$EFF....AA)

i=1
for (i in 1:length(mutations_df$EFF....AA)) {
  if (mutations_df$EFF....AA[i]=="") {
    mutations_df$EFF....AA[i] <- "NNo aminoacid change"
  } 
}

mutations_df = mutations_df[mutations_df$EFF....AA !="",]
pvalues_column <- as.character(formatC(mutations_df$PVALUE,format = "e",digits=1))

i=1
for (i in 1:length(pvalues_column)) {
  if (pvalues_column[i]=="0.000e+00") {
    pvalues_column[i]<-"0"
  }
}

pvalues_column <- as.character(formatC(mutations_df$PVALUE,format = "e",digits=3))
refcodon_column <- as.character(mutations_df$REFCODON)
refcodon_column[is.na(refcodon_column)]<-"0"
altcodon_column <- as.character(mutations_df$ALTCODON)
altcodon_column[is.na(altcodon_column)]<-"0"

sample_mutations= mutations_df$EFF....AA
sample_mutations = substring(sample_mutations, 2)
mutation_position = mutations_df$POS
sub_sample_mutations <- stri_extract(sample_mutations,regex="[a-zA-Z][a-zA-Z-][a-zA-Z][0-9]*")

i=1
codon=""
no_changes=0
changed_aa=""
for (i in 1:length(sample_mutations)) {
  if(sub_sample_mutations[i]=="ami"){
    no_changes=no_changes+1
  } else {
    #two consecutive mutations
    if (sub_sample_mutations[i]==sub_sample_mutations[i+1] & mutation_position[i]+1==mutation_position[i+1] ) {
      if (sub_sample_mutations[i+1]==sub_sample_mutations[i+2] & mutation_position[i+1]+1 == mutation_position[i+2]) {
        #changes in three nucleotides, same codon.
        codon <- paste(as.character(ref_alt_nuc[i]),as.character(ref_alt_nuc[i+1]),as.character(ref_alt_nuc[i+2]),sep = "")
        changed_aa<-GENETIC_CODE[codon][[1]]
        sample_mutations[i]<-paste(sub_sample_mutations[i],changed_aa,sep="")
        sample_mutations[i+1]<-paste(sub_sample_mutations[i+1],changed_aa,sep="")
        sample_mutations[i+2]<-paste(sub_sample_mutations[i+2],changed_aa,sep="")
      } else if (sub_sample_mutations[i]==sub_sample_mutations[i+1] & sub_sample_mutations[i]!=sub_sample_mutations[i-1] & mutation_position[i]+2 != mutation_position[i+2]){
        #changes in nuc1 and nuc2 but not nuc3
        codon<- paste(substring(altcodon_column[i],1,1),substring(altcodon_column[i+1],2,2),substring(refcodon_column[i],3,3),sep="")
        changed_aa<-GENETIC_CODE[codon][[1]]
        sample_mutations[i]<-paste(sub_sample_mutations[i],changed_aa,sep="")
        sample_mutations[i+1]<-paste(sub_sample_mutations[i+1],changed_aa,sep="")
      } else if (mutation_position[i]-1 != mutation_position[i-1]) {
        #changes in nuc2 and nuc3 but not nuc1
        codon<- paste(substring(refcodon_column[i],1,1),substring(altcodon_column[i+1],2,2),substring(altcodon_column[i],3,3),sep="")
        changed_aa<-GENETIC_CODE[codon][[1]]
        sample_mutations[i]<-paste(sub_sample_mutations[i],changed_aa,sep="")
        sample_mutations[i+1]<-paste(sub_sample_mutations[i+1],changed_aa,sep="")
      }
    } else if (sub_sample_mutations[i]==sub_sample_mutations[i+1] & mutation_position[i]+2==mutation_position[i+1]) {
      #changes in nuc1 and nuc3 but not nuc2
      codon<- paste(substring(altcodon_column[i],1,1),substring(refcodon_column[i],2,2),substring(altcodon_column[i+1],3,3),sep="")
      changed_aa<-GENETIC_CODE[codon][[1]]
      sample_mutations[i]<-paste(sub_sample_mutations[i],changed_aa,sep="")
      sample_mutations[i+1]<-paste(sub_sample_mutations[i+1],changed_aa,sep="")
    }
  }
}

sample_mutations