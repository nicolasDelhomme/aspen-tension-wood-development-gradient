library(dplyr)
library(here)
files<-sort(list.files(here("data/raw"), pattern="*.fastq.gz", full.names=TRUE))
md5 <- readr::read_table("data/raw/MD5.txt",
                         col_names=c("MD5","file"),
                         show_col_types=FALSE) %>%
  mutate(file=basename(file))

readr::write_csv(tibble::tibble(
  ExperimentTitle="Aspen tension wood development gradient",
  SampleName=sub("_[1,2].fastq.gz","",basename(files)),
  SampleDescription="",
  SequencingDate=sprintf("20%s-%s-%sT10:00:00",
                         substr(basename(files),3,4),
                         substr(basename(files),5,6),
                         substr(basename(files),7,8)),
  FileName=basename(files),
  FileLocation=dirname(files),
  MD5checksum=unlist(md5[match(basename(files),md5$file),"MD5"],use.names=FALSE)),
file="doc/samples.csv")
  
