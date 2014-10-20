### =========================================================================
### sortExonsByRank()
### -------------------------------------------------------------------------
###
### FIXME: Current implementation is a quick-and-dirty one that leverages the
### work done in .makeUCSCTxListFromGRangesList(). Therefore it looses the
### metadata columns. We need something better, faster, and documented (and it
### should probably go somewhere else).
### TODO: Also, maybe exonsBy(... , by="tx") should get a
### 'by.decreasing.rank.on.minus.strand' arg or something like that.
###

sortExonsByRank <- function(x, decreasing.rank.on.minus.strand=FALSE)
{
    .Defunct()
}

