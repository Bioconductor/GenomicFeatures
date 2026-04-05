#include "GenomicFeatures.h"

#include <R_ext/Rdynload.h>

#define CALLMETHOD_DEF(fun, numArgs) {#fun, (DL_FUNC) &fun, numArgs}

static const R_CallMethodDef callMethods[] = {

/* transcriptLocs2refLocs.c */
	CALLMETHOD_DEF(C_transcript_widths, 2),
	CALLMETHOD_DEF(C_tlocs2rlocs, 6),

	{NULL, NULL, 0}
};

void R_init_GenomicFeatures(DllInfo *info)
{
	R_registerRoutines(info, NULL, callMethods, NULL, NULL);
	return;
}

