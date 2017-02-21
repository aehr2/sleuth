wd  <- '/Volumes/work512/data/staff_agbio/'
setwd(wd)

########## ########## ########## ########## ########## 
########## ########## ########## ########## ########## 

library("RBPs", lib.loc="~/Rlibs")
DEGSETS      <- get.DEGsets()

SamplesRNA   <- DEGSETS$SamplesRNA
DEGsetsRNA   <- DEGSETS$DEGsetsRNA
DEGsetsRNATP <- DEGSETS$DEGsetsRNATP
SamplesRIP   <- DEGSETS$SamplesRIP
DEGsetsRIP   <- DEGSETS$DEGsetsRIP

##########
SUB='RNA'
##########

SamplesSUB <- switch(SUB,
                     "RNA"     = SamplesRNA  ,
                     "RNATP"   = SamplesRNA ,
                     "RIP"     = SamplesRIP
)

DEGsets <- switch(SUB,
                  "RNA"     = DEGsetsRNA ,
                  "RNATP"   = DEGsetsRNATP ,
                  "RIP"     = DEGsetsRIP 
)


########## ########## ########## ########## ########## 
TYPE='salmon'
print(TYPE)
#######

library("RBPs", lib.loc="~/Rlibs")
Samples <- get.Samples()

folder  <- switch(TYPE,
                  "salmon"   = paste0(wd,'RIP/salmon/STAR_TAIR10.atRTD.',as.character(Samples$ID) ) , 
                  "kallisto" = paste0(wd,'RIP/kallisto/TAIR10.atRTD.',as.character(Samples$ID) )
)

end     <- switch(TYPE,
                  "salmon"   =  c( rep('_R1single.Q.30.transcripts_quant0.6.1',sum(Samples$Experiment=='RNASeq')), rep(".Q.20.transcripts_quant0.6.1",sum(Samples$Experiment=='RIPSeq')) ) ,
                  "kallisto" =  c( rep('_R1single.Q.30.transcripts_quant',sum(Samples$Experiment=='RNASeq')), rep(".Q.20.transcripts_quant",sum(Samples$Experiment=='RIPSeq')) )
)

DIR <- switch(TYPE,
              "salmon"   = paste0(folder,end),
              "kallisto" = paste0(folder,end)
)
names(DIR) <- as.character(Samples$ID) 



sets <-"LL36__Col-2__vs__D1c"
SET <- DEGsets[,sets]

SamplesTwoSet  <- switch(SUB,
                         "RNA"     = SamplesSUB[SamplesSUB$TP  == SET["TP"]   & ( SamplesSUB$SID == SET["SID1"] | SamplesSUB$SID == SET["SID2"]  ),] ,
                         "RNATP"   = SamplesSUB[SamplesSUB$SID == SET["SID"]  & ( SamplesSUB$TP  == SET["TP1"]  | SamplesSUB$TP  == SET["TP2"]   ),] ,
                         "RIP"     = SamplesSUB[SamplesSUB$TP  == SET["TP"]   & ( SamplesSUB$SID == SET["SID1"] | SamplesSUB$SID == SET["SID2"]  ),]  
)


s3c <- data.frame( sample    = SamplesTwoSet$ID,
                   condition = SamplesTwoSet$Group)

s3c <- dplyr::mutate(s3c, path = DIR[SamplesTwoSet$ID] )


TXI <- readRDS(paste0("/Volumes/Macintosh HD/Users/weinhol/Promotion/ServerTMP/TXI.RDS") )
TXIrn <- rownames(TXI$Transcript$abundance)

library("sleuthDev", lib.loc="~/Rlibs")
so3 <- sleuthDev::sleuth_prep(s3c, ~ condition,max_bootstrap = 100,TXIrn = TXIrn )#,filter_fun = TXIbasic_filter)
so3 <- sleuth_fit(so3)
so.wt2.C <- sleuth_wt(so3, which_beta = paste0('condition','D1c_LL36') )
so.res.wt2.C<-sleuth_results(so.wt2.C,test = paste0('condition','D1c_LL36'),  test_type = "wt", which_model = "full", rename_cols = TRUE, show_all = TRUE)
rownames(so.res.wt2.C) <- as.character(so.res.wt2.C$target_id)
sum(!is.na(so.res.wt2.C[f.input2(as.character(so.res.wt2.C$target_id),TXIrn)$diffAB,]$pval))



