## Import each subject's results file and concatenate into two data frames:
# 1. A full data frame for all trials (uneven distribution: 16 1Talker, 16 2Talker, 24 SSN, 8 8Talker)
# 2. A balanced data frame with the first n of all conditions, where n is the smallest condition

path = "~/Documents/Projects/EAB/outputs_with_block/"
setwd(path)

full.file<-""
sub.file<-""

file.names <- dir(path, pattern ="*.csv")
for(i in 1:length(file.names)){
  file <- read.csv(file.names[i],header=TRUE, sep=",", stringsAsFactors=TRUE)
  conditions <- levels(file$condition)
  
  # put in block order (trial #, temporally ordered) just in case it isn't
  file <- file[order(file$block),]
  
  # strip out nan rows
  file <- file[!is.na(file$condition),]
  
  # minimum number of trials / cond
  num_trials = min(xtabs(~condition,file))
  
  file.subset <-""
  for(j in 1:length(conditions)){
    this.cond <- file[file$condition==conditions[j],]
    file.subset <- rbind(file.subset,this.cond[c(1:num_trials),])
  }
  file.subset <- file.subset[order(file.subset$block),]
  # strip out nan rows
  file.subset <- file.subset[!is.na(file.subset$condition),]
  print(file.names[i])
  print(xtabs(~condition,file))
  print(xtabs(~condition,file.subset))
  
  full.file <- rbind(full.file, file)
  sub.file <- rbind(sub.file,file.subset)
}

# Write out tables for posterity
write.table(full.file, file = "merge_full.txt",sep=",", 
            row.names = FALSE, qmethod = "double",fileEncoding="windows-1252")
write.table(sub.file, file = "merge_subset.txt",sep=",", 
          row.names = FALSE, qmethod = "double",fileEncoding="windows-1252")

## Perform analyses
