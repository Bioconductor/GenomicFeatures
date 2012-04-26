transcripts_deprecated <- function(genes, proximal = 500L, distal = 10000L) {
  msg <- c("  transcripts_deprecated() is defunct.\n",
           "  Please use transcripts() on a TranscriptDb object instead ",
           "(see ?transcripts).")
  .Defunct(msg=msg)
}


.countPos <- function(pos) {
   strsplit(as.character(pos), ",", fixed=TRUE)
}

.splitPos <- function(pos) {
    as.integer(unlist(.countPos(pos), use.names=FALSE))
}


## Introns are stored in the local DB ... for now

exons_deprecated <- function(genes) {
  msg <- c("  exons_deprecated() is defunct.\n",
           "  Please use exons() on a TranscriptDb object instead ",
           "(see ?exons).")
  .Defunct(msg=msg)
}


introns_deprecated <- function(genes) {
  msg <- c("  introns_deprecated() is defunct.\n",
           "  Please use intronsByTranscript() on a TranscriptDb object instead ",
           "(see ?intronsByTranscript).")
  .Defunct(msg=msg)
}
