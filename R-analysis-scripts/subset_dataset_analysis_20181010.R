################################
## Import libraries ############
library(lme4)
library(lmerTest)
library(plyr)
library(effsize)

################################
## Prep the datasets ###########
# Import data and clean up
mydata <- read.csv('merge_subset.txt')
mydata$group[mydata$subject_id<2000] <- 'native'
mydata$group[mydata$subject_id>=2000] <- 'nonnative'
mydata$condition <- relevel(mydata$condition,'SSN')

# Expand and rebind for morphosyntactic and content word errors
morphdat <- mydata[mydata$DNH!=1,]
morphdat$error_type <- 'morph'
morphdat$error_rate <- morphdat$total_morph_errors/morphdat$total.num.words
worddat <- mydata[mydata$DNH!=1,]
worddat$error_type <- 'content'
worddat$error_rate <- worddat$total_cont_errors/worddat$total.num.words
mydata_errors <- rbind(morphdat,worddat)

# Build the subject-level datasets
mydata_subj <- ddply(mydata_errors,c('group','subject_id','condition','error_type'),
                     summarise,
                     error_rate = mean(error_rate))

################################
## Perform model analyses ######

# Model of DNH errors
dnh_model <- glmer(DNH~(1|subject_id)+group*condition,mydata,family='binomial')
summary(dnh_model)
anova(dnh_model)

# Omnibus Model with Error Type
omni_model <- lmer(error_rate~(1|subject_id)+group*condition*error_type,mydata_errors[mydata_errors$error_type!='dnh',])
summary(omni_model)
anova(omni_model)

# Word-Level Analysis
contentword <- lmer(error_rate~(1|subject_id)+group*condition,worddat)
summary(contentword)
anova(contentword)

# Morpheme-Level Analysis
morphemes <- lmer(error_rate~(1|subject_id)+group*condition,morphdat)
summary(morphemes)
anova(morphemes)


################################
## Planned and post-hoc tests ##

# Subject-level posthocs on error rates for native vs. non-native
# Average across the mask type (condition)-level averages to get subject-level average
mydata_subj_bal <- ddply(mydata_subj, c('group','subject_id','error_type'),
                         summarise, 
                         error_rate = mean(error_rate))

t.test(
  mydata_subj_bal$error_rate[mydata_subj_bal$group=='native' & mydata_subj_bal$error_type=='content'],
  mydata_subj_bal$error_rate[mydata_subj_bal$group=='nonnative' & mydata_subj_bal$error_type=='content']
)
cohen.d(
  mydata_subj_bal$error_rate[mydata_subj_bal$group=='native' & mydata_subj_bal$error_type=='content'],
  mydata_subj_bal$error_rate[mydata_subj_bal$group=='nonnative' & mydata_subj_bal$error_type=='content'],
  na.rm=TRUE
)

t.test(
  mydata_subj_bal$error_rate[mydata_subj_bal$group=='native' & mydata_subj_bal$error_type=='morph'],
  mydata_subj_bal$error_rate[mydata_subj_bal$group=='nonnative' & mydata_subj_bal$error_type=='morph']
)
cohen.d(
  mydata_subj_bal$error_rate[mydata_subj_bal$group=='native' & mydata_subj_bal$error_type=='morph'],
  mydata_subj_bal$error_rate[mydata_subj_bal$group=='nonnative' & mydata_subj_bal$error_type=='morph'],
  na.rm=TRUE
)

# Subject-level posthocs for mask types
mydata_paired = data.frame(subject_id = unique(mydata_subj$subject_id),Mask1T=NA,Mask2T=NA,Mask8T=NA,MaskSSN=NA)

for(ss in 1:length(mydata_paired$subject_id)){
  mydata_paired[ss,]$Mask1T <- mydata_subj$error_rate[mydata_subj$subject_id==mydata_paired$subject_id[ss] & mydata_subj$error_type=='content' & mydata_subj$condition=='1T'][1]
  mydata_paired[ss,]$Mask2T <- mydata_subj$error_rate[mydata_subj$subject_id==mydata_paired$subject_id[ss] & mydata_subj$error_type=='content' & mydata_subj$condition=='2Talker'][1]
  mydata_paired[ss,]$Mask8T <- mydata_subj$error_rate[mydata_subj$subject_id==mydata_paired$subject_id[ss] & mydata_subj$error_type=='content' & mydata_subj$condition=='8Talker'][1]
  mydata_paired[ss,]$MaskSSN <- mydata_subj$error_rate[mydata_subj$subject_id==mydata_paired$subject_id[ss] & mydata_subj$error_type=='content' & mydata_subj$condition=='SSN'][1]
}
t.test(mydata_paired$Mask1T,
       mydata_paired$MaskSSN,
       paired=TRUE)
t.test(mydata_paired$Mask2T,
       mydata_paired$MaskSSN,
       paired=TRUE)
t.test(mydata_paired$Mask8T,
       mydata_paired$MaskSSN,
       paired=TRUE)

for(ss in 1:length(mydata_paired$subject_id)){
  mydata_paired[ss,]$Mask1T <- mydata_subj$error_rate[mydata_subj$subject_id==mydata_paired$subject_id[ss] & mydata_subj$error_type=='morph' & mydata_subj$condition=='1T'][1]
  mydata_paired[ss,]$Mask2T <- mydata_subj$error_rate[mydata_subj$subject_id==mydata_paired$subject_id[ss] & mydata_subj$error_type=='morph' & mydata_subj$condition=='2Talker'][1]
  mydata_paired[ss,]$Mask8T <- mydata_subj$error_rate[mydata_subj$subject_id==mydata_paired$subject_id[ss] & mydata_subj$error_type=='morph' & mydata_subj$condition=='8Talker'][1]
  mydata_paired[ss,]$MaskSSN <- mydata_subj$error_rate[mydata_subj$subject_id==mydata_paired$subject_id[ss] & mydata_subj$error_type=='morph' & mydata_subj$condition=='SSN'][1]
}
t.test(mydata_paired$Mask1T,
       mydata_paired$MaskSSN,
       paired=TRUE)
t.test(mydata_paired$Mask2T,
       mydata_paired$MaskSSN,
       paired=TRUE)
t.test(mydata_paired$Mask8T,
       mydata_paired$MaskSSN,
       paired=TRUE)
