### =========================================================================
### Extend exons by a given number of bases into their adjacent introns
### -------------------------------------------------------------------------
###

extendExonsIntoIntrons <- function(ex_by_tx, extent=2)
{
    if (!is(ex_by_tx, "GRangesList"))
        stop(wmsg("'ex_by_tx' must be a GRangesList object"))
    if (!isSingleNumber(extent))
        stop(wmsg("'extent' must be a single number"))
    if (!is.integer(extent))
        extent <- as.integer(extent)

    resize_idx <- which(lengths(ex_by_tx) >= 2L)
    ex_to_resize <- ex_by_tx[resize_idx]

    ## Resize first exons.
    first_ex <- heads(ex_to_resize, n=1L)
    unlisted <- unlist(first_ex, use.names=FALSE)
    unlisted <- resize(unlisted, width(unlisted) + extent,
                       fix="start", use.names=FALSE)
    first_ex <- relist(unlisted, first_ex)

    ## Resize last exons.
    last_ex <- tails(ex_to_resize, n=1L)
    unlisted <- unlist(last_ex, use.names=FALSE)
    unlisted <- resize(unlisted, width(unlisted) + extent,
                       fix="end", use.names=FALSE)
    last_ex <- relist(unlisted, last_ex)

    ## Resize intermediate exons.
    mid_ex <- tails(heads(ex_to_resize, n=-1L), n=-1L)
    unlisted <- unlist(mid_ex, use.names=FALSE)
    unlisted <- resize(unlisted, width=width(unlisted) + 2L*extent,
                       fix="center", use.names=FALSE)
    mid_ex <- relist(unlisted, mid_ex)

    ## Put exons back together.
    ex_by_tx[resize_idx] <- pc(first_ex, mid_ex, last_ex)
    ex_by_tx
}

