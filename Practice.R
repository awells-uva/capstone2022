library(plotly)
library(cluster)
library(stringr)
library(caret)
library(broom)
library(regclass)
library(CustomerScoringMetrics)
library(aod)
library(OneR)
library(smotefamily)
library(dplyr)
library(scales)
# library(factoextra)
# library(FactoMineR)
# library(nFactors)
# library(GPArotation)
library(clValid)
library(stats)
library(MLmetrics)
library(yardstick)
library(maditr)
library(ape)
library(seqinr)
library(kmer)

seq_1_DNAbin <- read.GenBank("JF806202")
attr(seq_1_DNAbin, "species") 
seq_1_DNAbin$JF806202
str(seq_1_DNAbin)


#save as character object:
seq_1_character <- read.GenBank("JF806202", as.character = TRUE)
seq_1_character #this is not a very nice format


#vector of lizard accession numbers
lizards_accession_numbers <- c("JF806202", "HM161150", "FJ356743", "JF806205",
                               "JQ073190", "GU457971", "FJ356741", "JF806207",
                               "JF806210", "AY662592", "AY662591", "FJ356748",
                               "JN112660", "AY662594", "JN112661", "HQ876437",
                               "HQ876434", "AY662590", "FJ356740", "JF806214",
                               "JQ073188", "FJ356749", "JQ073189", "JF806216",
                               "AY662598", "JN112653", "JF806204", "FJ356747",
                               "FJ356744", "HQ876440", "JN112651", "JF806215",
                               "JF806209")

#read sequences and place them in a DNAbin object
lizards_sequences <- read.GenBank(lizards_accession_numbers) 
#a brief summary of what is in the object, including base composition
lizards_sequences 
#a list of the DNAbin elements with length of the sequences
#notice the one of the attributes is the species names
str(lizards_sequences) 


#see the list of attributes and contents
attributes(lizards_sequences)
#the accession numbers
names(lizards_sequences) 
# we get the species list. Notice this
attr(lizards_sequences, "species") 
# attr is slightly different function


lizards_sequences_GenBank_IDs <- paste(attr(lizards_sequences, "species"), names
                                       (lizards_sequences), sep ="_RAG1_") 

lizards_sequences_GenBank_IDs


setwd(dir = "C:/Users/MSachs.MSACHS-DELL/Documents/UVA MSDS/DS 6011/FASTA Files/")
### we are going to write in fasta format
write.dna(lizards_sequences, file ="lizard_fasta_1.fasta", format = "fasta", append =
            FALSE, nbcol = 6, colsep = " ", colw = 10)

lizard_seq_seqinr_format <- read.fasta(file = "lizard_fasta_1.fasta", seqtype = "DNA",
                                       as.string = TRUE, forceDNAtolower = FALSE)
lizard_seq_seqinr_format #this shows dif


write.fasta(sequences = lizard_seq_seqinr_format, names = lizards_sequences_GenBank_IDs,
            nbchar = 10, file.out = "lizard_seq_seqinr_format.fasta")

lizard_seq_seqinr_format <- read.fasta(file = "lizard_seq_seqinr_format.fasta", seqtype =
                                               "DNA", as.string = TRUE, forceDNAtolower = FALSE)
lizard_seq_seqinr_format


data(woodmouse, package = "ape")
ape::as.character.DNAbin(woodmouse[1:5, 1:5])

woodmouse <- woodmouse[, apply(woodmouse, 2, function(v) !any(v == 0xf0))]
woodmouse.kdist <- kdistance(woodmouse, k = 6)
print(as.matrix(woodmouse.kdist)[1:7, 1:7], digits = 2)


suppressWarnings(RNGversion("3.5.0"))
set.seed(999)
seeds <- sample(1:15, size = 3)
woodmouse.mbed <- mbed(woodmouse, seeds = seeds, k = 6)
print(woodmouse.mbed[,], digits = 2)



## compute pairwise distance matrices
dist1 <- ape::dist.dna(woodmouse, model = "K80") 
dist2 <- kdistance(woodmouse, k = 7) 

## build neighbor-joining trees
phy1 <- ape::nj(dist1)
phy2 <- ape::nj(dist2)

## rearrange trees in ladderized fashion
phy1 <- ape::ladderize(phy1)
phy2 <- ape::ladderize(phy2)

## convert phylo objects to dendrograms
dnd1 <- as.dendrogram(phy1)
dnd2 <- as.dendrogram(phy2)

## plot the tanglegram
dndlist <- dendextend::dendlist(dnd1, dnd2)
dendextend::tanglegram(dndlist, fast = TRUE, margin_inner = 5)
