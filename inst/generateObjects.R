library(crisprBase)
library(devtools)

nnn_pams <- c("AA", "AC", "AG", "AT",
              "CA", "CC", "CG", "CT",
              "GA", "GC", "GG", "GT",
              "TA", "TC", "TG", "TT")
nnn_pams <- paste0("(3/3)N", nnn_pams)

cas9_pam_weights <- c("0.000", "0.000", "0.259", "0.000",
                      "0.000", "0.000", "0.107", "0.000",
                      "0.069", "0.022", "1.000", "0.016",
                      "0.000", "0.000", "0.039", "0.000")

cas9ng_pam_weights <- c("0.138", "0.092", "0.415", "0.305",
                        "0.000", "0.000", "0.007", "0.000",
                        "0.870", "0.694", "0.952", "1.000",
                        "0.017", "0.007", "0.092", "0.007")

spry_pam_weights <- c("0.514", "0.472", "1.000", "0.611",
                      "0.208", "0.347", "0.653", "0.778",
                      "0.653", "0.764", "0.806", "0.319",
                      "0.014", "0.347", "0.319", "0.042")

spryc_pam_weights <- c("0.470", "0.614", "0.699", "0.434",
                       "0.964", "0.614", "0.602", "0.542",
                       "0.614", "0.675", "0.542", "1.000",
                       "0.578", "0.229", "0.470", "0.277")




SpCas9 <- CrisprNuclease("SpCas9",
                         pams=nnn_pams,
                         weights=cas9_pam_weights,
                         metadata=list(description="Wildtype Streptococcus pyogenes Cas9 (SpCas9) nuclease"),
                         pam_side="3prime",
                         spacer_length=20)

SpCas9NG <- CrisprNuclease("SpCas9NG",
                         pams=nnn_pams,
                         weights=cas9ng_pam_weights,
                         metadata=list(description="relaxed-PAM engineered SpCas9 variant (Nishimasu 2018)."),
                         pam_side="3prime",
                         spacer_length=20)

SpRY <- CrisprNuclease("SpRY",
                            pams=nnn_pams,
                            weights=spry_pam_weights,
                            metadata=list(description="Near-PAMless engineered SpCas9 variant (Walton 2020)."),
                            pam_side="3prime",
                            spacer_length=20)

SpRYc <- CrisprNuclease("SpRYc",
                            pams=nnn_pams,
                            weights=spryc_pam_weights,
                            metadata=list(description="SpRY PID with Sc++ N-terminus (Zhao 2023)."),
                            pam_side="3prime",
                            spacer_length=20)


SaCas9 <- CrisprNuclease("SaCas9",
                         pams=c("(3/3)NNGRRA","(3/3)NNGRRC", "(3/3)NNGRRG", "(3/3)NNGRRT"),
                         weights=c(0.1666, 0.1666, 0.1666, 0.5),
                         metadata=list(description="Wildtype Staphylococcus aureus Cas9 (SaCas9) nuclease"),
                         pam_side="3prime",
                         spacer_length=21)


AsCas12a <- CrisprNuclease("AsCas12a",
                           pams="TTTV(18/23)",
                           metadata=list(description="Wildtype Acidaminococcus Cas12a (AsCas12a) nuclease."),
                           pam_side="5prime",
                           spacer_length=23)

MAD7 <- CrisprNuclease("MAD7",
                       pams="YTTV(18/23)",
                       metadata=list(description="MAD7 nuclease (Cas12a-like family)."),
                       pam_side="5prime",
                       spacer_length=23)


### RNA-targeting nucleases:
CasRx <- CrisprNuclease("CasRx",
                        targetType="RNA",
                        pams="N",
                        metadata=list(description="Cas13d-NLS from Ruminococcus flavefaciens strain XPD3002.",
                                      doi="10.1016/j.cell.2018.02.033"),
                        pam_side="3prime",
                        spacer_length=23)

Csm <- CrisprNuclease("Csm",
                      targetType="RNA",
                      pams="N",
                      metadata=list(description="RNA-targeting Csm complex from Streptococcus thermophilus",
                                    doi="10.1038/s41587-022-01649-9"),
                      pam_side="3prime",
                      spacer_length=32)


# SpG nuclease:
load("pams/spg/spg.pams.rda")
pams <- spg.pams[, c("PAM", "Score")]
pams <- pams[order(-pams$Score),,]
motifs <- paste0("(3/3)", pams$PAM)
SpG <- CrisprNuclease("SpG",
                             pams=motifs,
                             weights=pams$Score,
                             metadata=list(description="Engineered Streptococcus pyogenes Cas9 nuclease SpCas9NG (SpG)."),
                             pam_side="3prime",
                             spacer_length=20)


load("pams/cas12a/cas12a.pams.rda")
pams <- cas12a.pams[, c("PAM", "Score_Doench")]
pams <- pams[!is.na(pams$Score_Doench),]
pams$Score_Doench[pams$PAM %in% c("TTTC","TTTA", "TTTG")] <- 1
pams <- pams[order(-pams$Score_Doench),,]
motifs <- paste0(pams$PAM, "(18/23)")
enAsCas12a <- CrisprNuclease("enAsCas12a",
                             pams=motifs,
                             weights=pams$Score_Doench,
                             metadata=list(description="Enhanced Acidaminococcus Cas12a (AsCas12a) nuclease."),
                             pam_side="5prime",
                             spacer_length=23)


#Generate base editor
ws <- t(read.csv("../inst/be/b4max.csv"))
colnames(ws) <- ws["Position",]
ws <- ws[-c(match("Position", rownames(ws))),]
load("../data/SpCas9.rda")
BE4max <- BaseEditor(SpCas9,
                     baseEditorName="BE4max",
                     editingStrand="original",
                     editingWeights=ws)
metadata(BE4max)$description_base_editor <- "BE4max cytosine base editor."


#dir.create("../data")
use_data(SpCas9,
         SpCas9NG,
         SpRY,
         SpRYc,
         SaCas9,
         SpG,
         AsCas12a,
         enAsCas12a,
         MAD7,
         CasRx,
         Csm,
         BE4max,
         compress="xz", internal=FALSE, overwrite=TRUE)





# Enzymes:
# EcoRI <- Nuclease("EcoRI",
#                   motifs=c("G^AATTC"),
#                   metadata="EcoRI restriction enzyme")

# SmaI <- Nuclease("SmaI",
#                   motifs=c("CCC^GGG"),
#                   metadata="SmaI restriction enzyme")

# HgaI <- Nuclease("HgaI",
#                  motifs=c("GACGC(5/10)"),
#                  metadata="HgaI restriction enzyme")

# PfaAI <- Nuclease("PfaAI",
#                   motifs=c("G^GYRCC"),
#                   metadata="PfaAI restriction enzyme")








