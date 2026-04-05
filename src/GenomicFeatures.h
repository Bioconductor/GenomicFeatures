#include <Rdefines.h>


/* transcriptLocs2refLocs.c */

SEXP C_transcript_widths(
	SEXP exonStarts,
	SEXP exonEnds
);

SEXP C_tlocs2rlocs(
	SEXP tlocs,
	SEXP exonStarts,
	SEXP exonEnds,
	SEXP strand,
	SEXP decreasing_rank_on_minus_strand,
	SEXP error_if_out_of_bounds
);

