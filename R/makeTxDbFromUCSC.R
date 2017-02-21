### =========================================================================
### makeTxDbFromUCSC()
### -------------------------------------------------------------------------

## take any of the supported genomes at UCSC, remove the version number, 
## and then lookup the supported organism name.
UCSCGenomeToOrganism <- function(genome){
  genome <- gsub("\\d+$","",genome)
  genome2org <- c("hg"="Homo sapiens",
                 "felCat"="Felis catus",
                    "galGal"="Gallus gallus",
                    "panTro"="Pan troglodytes",
                    "bosTau"="Bos taurus",
                    "canFam"="Canis familiaris",
                    "loxAfr"="Loxodonta africana",
                    "fr"="Fugu rubripes",
                    "cavPor"="Cavia porcellus",
                    "equCab"="Equus caballus",
                    "petMar"="Petromyzon marinus",
                    "anoCar"="Anolis carolinensis",
                    "calJac"="Callithrix jacchus",
                    "oryLat"="Oryzias latipes",
                    "mm"="Mus musculus",
                    "monDom"="Monodelphis domestica",
                    "ponAbe"="Pongo abelii",
                    "ailMel"="Ailuropoda melanoleuca",
                    "susScr"="Sus scrofa",
                    "ornAna"="Ornithorhynchus anatinus",
                    "oryCun"="Oryctolagus cuniculus",
                    "rn"="Rattus norvegicus",
                    "rheMac"="Macaca mulatta",
                    "oviAri"="Ovis aries",
                    "gasAcu"="Gasterosteus aculeatus",
                    "tetNig"="Tetraodon nigroviridis",
                    "xenTro"="Xenopus tropicalis",
                    "taeGut"="Taeniopygia guttata",
                    "danRer"="Danio rerio",
                    "ci"="Ciona intestinalis",
                    "braFlo"="Branchiostoma floridae",
                    "strPur"="Strongylocentrotus purpuratus",
                    "apiMel"="Apis mellifera",
                    "anoGam"="Anopheles gambiae",
                    "droAna"="Drosophila ananassae",
                    "droEre"="Drosophila erecta",
                    "droGri"="Drosophila grimshawi",
                    "dm"="Drosophila melanogaster",
                    "droMoj"="Drosophila mojavensis",
                    "droPer"="Drosophila persimilis",
                    "dp"="Drosophila pseudoobscura",
                    "droSec"="Drosophila sechellia",
                    "droSim"="Drosophila simulans",
                    "droVir"="Drosophila virilis",
                    "droYak"="Drosophila yakuba",
                    "caePb"="Caenorhabditis brenneri",
                    "cb"="Caenorhabditis briggsae",
                    "ce"="Caenorhabditis elegans",
                    "caeJap"="Caenorhabditis japonica",
                    "caeRem"="Caenorhabditis remanei",
                    "priPac"="Pristionchus pacificus",
                    "aplCal"="Aplysia californica",
                    "sacCer"="Saccharomyces cerevisiae",
                    "papAnu"="Papio anubis",
                    "vicPac"= "Vicugna pacos",
                    "dasNov"= "Dasypus novemcinctus",
                    "otoGar"= "Otolemur garnettii",
                    "papHam"= "Papio hamadryas",
                    "papAnu"= "Papio anubis",
                    "turTru"= "Tursiops truncatus",
                    "nomLeu"= "Nomascus leucogenys",
                    "gorGor"= "Gorilla gorilla",
                    "eriEur"= "Erinaceus europaeus",
                    "dipOrd"= "Dipodomys ordii",
                    "triMan"= "Trichechus manatus",
                    "pteVam"= "Pteropus vampyrus",
                    "myoLuc"= "Myotis lucifugus",
                    "micMur"= "Microcebus murinus",
                    "hetGla"= "Heterocephalus glaber",
                    "ochPri"= "Ochotona princeps",
                    "proCap"= "Procavia capensis",
                    "sorAra"= "Sorex  araneus",
                    "choHof"= "Choloepus hoffmanni",
                    "speTri"= "Spermophilus tridecemlineatus",
                    "saiBol"= "Saimiri sciureus",
                    "sorAra"= "Sorex araneus",
                    "sarHar"= "Sarcophilus harrisii",
                    "echTel"= "Echinops telfairi",
                    "tupBel"= "Tupaia belangeri",
                    "macEug"= "Macropus eugenii",
                    "cerSim"= "Ceratotherium simum",
                    "gadMor"= "Gadus morhua",
                    "melUnd"= "Melopsittacus undulatus",
                    "latCha"= "Latimeria chalumnae",
                    "geoFor"= "Geospiza fortis",
                    "oreNil"= "Oreochromis niloticus",
                    "chrPic"= "Chrysemys picta",
                    "melGal"= "Meleagris gallopavo",
                    "panPan"= "Pan paniscus",
                    "aptMan"= "Apteryx australis",
                    "criGri"= "Cricetulus griseus",
                    "macFas"= "Macaca fascicularis",
                    "musFur"= "Mustela putorius",
                    "chlSab"= "Chlorocebus sabaeus",
                    "galVar"= "Galeopterus variegatus",
                    "balAcu"= "Balaenoptera acutorostrata",
                    "tarSyr"= "Tarsius syrichta",
                    "allMis"= "Alligator mississippiensis",
                    "calMil"= "Callorhinchus milii",
                    "eboVir"= "Filoviridae ebolavirus"
                  )
  genome2org[genome]
}

### makeTxDbFromUCSC() expects a UCSC transcript table to have at least
### the following columns:
.UCSC_TXCOL2CLASS <- c(
    name="character",
    chrom="factor",
    strand="factor",
    txStart="integer",
    txEnd="integer",
    cdsStart="integer",
    cdsEnd="integer",
    exonCount="integer",
    exonStarts="character",
    exonEnds="character"
)
### Note that, from a strictly technical point of view, the 'name' and
### 'exonCount' cols are not required and .makeTxDbFromUCSCTxTable()
### could easily be modified to accept tables where they are missing.

### Lookup between UCSC transcript tables and their associated track.
.SUPPORTED_UCSC_TABLES <- c(
  ## tablename (unique key)           track                subtrack

  ## Tables/tracks shared by hg18/hg19.
  ## All the tables/tracks listed in this section belong to the "Genes and
  ## Gene Prediction" group of tracks for hg18 and hg19.
  ## On Aug 13 2010, makeTxDbFromUCSC() was successfully tested by hand on
  ## all of them for hg18 (i.e. with 'genome="hg18"').
  ## Note: the "acembly" table contains more than 250000 transcripts!
  "knownGene",                        "UCSC Genes",        NA,
  "knownGeneOld8",                    "Old UCSC Genes",    NA,
  "knownGeneOld7",                    "Old UCSC Genes",    NA,
  "knownGeneOld6",                    "Old UCSC Genes",    NA,
  "knownGeneOld4",                    "Old UCSC Genes",    NA,
  "knownGeneOld3",                    "Old UCSC Genes",    NA,
  "knownGenePrevious",                "Old Known Genes",   NA,
  "ccdsGene",                         "CCDS",              NA,
  "refGene",                          "RefSeq Genes",      NA,
  "xenoRefGene",                      "Other RefSeq",      NA,
  "vegaGene",                         "Vega Genes",        "Vega Protein Genes",
  "vegaPseudoGene",                   "Vega Genes",        "Vega Pseudogenes",
  "ensGene",                          "Ensembl Genes",     NA,
  "acembly",                          "AceView Genes",     NA,
  "sibGene",                          "SIB Genes",         NA,
  "nscanPasaGene",                    "N-SCAN",            "N-SCAN PASA-EST",
  "nscanGene",                        "N-SCAN",            "N-SCAN",
  "sgpGene",                          "SGP Genes",         NA,
  "geneid",                           "Geneid Genes",      NA,
  "genscan",                          "Genscan Genes",     NA,
  "exoniphy",                         "Exoniphy",          NA,
  "augustusGene",                     "Augustus",          NA,
  "augustusHints",                    "Augustus",          "Augustus Hints",
  "augustusXRA",                      "Augustus",          "Augustus De Novo",
  "augustusAbinitio",                 "Augustus",          "Augustus Ab Initio",
  "acescan",                          "ACEScan",           NA,
  "lincRNAsTranscripts",              "lincRNAsTranscripts", NA,

  ## Tables/tracks specific to hg18.
  "wgEncodeGencodeManualV3",          "Gencode Genes",     "Gencode Manual",
  "wgEncodeGencodeAutoV3",            "Gencode Genes",     "Gencode Auto",
  "wgEncodeGencodePolyaV3",           "Gencode Genes",     "Gencode PolyA",

  ## Tables/tracks specific to hg19.
  "wgEncodeGencodeBasicV19",          "GENCODE Genes V19", NA,
  "wgEncodeGencodeCompV19",           "GENCODE Genes V19", NA,
  "wgEncodeGencodePseudoGeneV19",     "GENCODE Genes V19", NA,
  "wgEncodeGencode2wayConsPseudoV19", "GENCODE Genes V19", NA,
  "wgEncodeGencodePolyaV19",          "GENCODE Genes V19", NA,

  "wgEncodeGencodeBasicV17",          "GENCODE Genes V17", NA,
  "wgEncodeGencodeCompV17",           "GENCODE Genes V17", NA,
  "wgEncodeGencodePseudoGeneV17",     "GENCODE Genes V17", NA,
  "wgEncodeGencode2wayConsPseudoV17", "GENCODE Genes V17", NA,
  "wgEncodeGencodePolyaV17",          "GENCODE Genes V17", NA,

  "wgEncodeGencodeBasicV14",          "GENCODE Genes V14", NA,
  "wgEncodeGencodeCompV14",           "GENCODE Genes V14", NA,
  "wgEncodeGencodePseudoGeneV14",     "GENCODE Genes V14", NA,
  "wgEncodeGencode2wayConsPseudoV14", "GENCODE Genes V14", NA,
  "wgEncodeGencodePolyaV14",          "GENCODE Genes V14", NA,

  "wgEncodeGencodeBasicV7",           "GENCODE Genes V7",  NA,
  "wgEncodeGencodeCompV7",            "GENCODE Genes V7",  NA,
  "wgEncodeGencodePseudoGeneV7" ,     "GENCODE Genes V7",  NA,
  "wgEncodeGencode2wayConsPseudoV7",  "GENCODE Genes V7",  NA,
  "wgEncodeGencodePolyaV7",           "GENCODE Genes V7",  NA,

  ## Tables/tracks specific to D. melanogaster.
  "flyBaseGene",                      "FlyBase Genes",     NA,

  ## Tables/tracks specific to sacCer2.
  ## makeTxDbFromUCSC(genome="sacCer2", tablename="sgdGene")
  ## successfully tested on On Aug 13 2010.
  "sgdGene",                          "SGD Genes",         NA
)

### Return a data.frame with 3 columns (tablename, track, and subtrack) and
### 1 row per UCSC table known to work with makeTxDbFromUCSC().
### A note about the current implementation:
### Current implementation uses hard-coded .SUPPORTED_UCSC_TABLES matrix
### above which is not satisfying in the long run (the matrix needs to be
### manually updated from times to times, a long and boring and error-prone
### process, and is probably out-of-sync at the moment). Ideally we'd like
### to be able to generate the 3-column data.frame programmatically in
### reasonable time. For this we need to be able to retrieve all the "central
### tables" for all the transcript-centric tracks available for a given
### organism. Using a combination of calls to rtracklayer::trackNames(session)
### and rtracklayer::tableNames(ucscTableQuery(session, track=track)) would
### partly achieve this but is unfortunately very slow.
supportedUCSCtables <- function(genome="hg19",
                                url="http://genome.ucsc.edu/cgi-bin/")
{
    if (is(genome, "UCSCSession")) {
        if (!missing(url))
            warning("'url' is ignored when 'genome' is a UCSCSession object")
    } else {
        if (!isSingleStringOrNA(genome))
            stop("'genome' must be a single string or NA")
        if (!isSingleString(url))
            stop("'url' must be a single string")
    }
    mat <- matrix(.SUPPORTED_UCSC_TABLES, ncol=3, byrow=TRUE)
    colnames(mat) <- c("tablename", "track", "subtrack")
    ans_tablename <- mat[ , "tablename"]
    ans_track <- mat[ , "track"]
    ans_subtrack <- mat[ , "subtrack"]
    ans <- data.frame(tablename=ans_tablename,
                      track=ans_track,
                      subtrack=ans_subtrack,
                      stringsAsFactors=FALSE)
    if (isSingleStringOrNA(genome) && is.na(genome)) {
        ans$track <- factor(ans_track, levels=unique(ans_track))
        return(ans)
    }
    if (is(genome, "UCSCSession")) {
        session <- genome
        genome <- genome(session)
    } else {
        session <- browserSession(url=url)
        genome(session) <- genome
    }
    if (genome %in% c("hg17", "hg16", "mm8", "mm7", "rn3")) {
        ans_track[ans$tablename == "knownGene"] <- "Known Genes"
        ans$track <- ans_track
    } else if (genome %in% "hg38") {
        ans_track[ans$tablename == "knownGene"] <- "GENCODE v24"
        ans$track <- ans_track
    }
    ## trackNames() returns a mapping from track names to "central table" names
    ## in the form of a named character vector where the names are the track
    ## names and the values the "central table" names (more than 1 table can
    ## be connected to a given track via joins thru the "central table").
    ## Unfortunately such mapping cannot handle the situation where a track is
    ## mapped to more than 1 "central table". This happens for example when a
    ## track has subtracks (e.g. the Augustus track for hg18 has 3 subtracks),
    ## in which case there is 1 "central table" per subtrack. So trackNames()
    ## alone cannot be used to get the one-to-many mapping from tracks to
    ## "central tables". Calling tableNames(ucscTableQuery(session,
    ## track=track)) in a loop on all the tracks returned by trackNames()
    ## would work but is very slow :-/
    genome_tracknames <- trackNames(session)
    ans <- ans[ans$track %in% names(genome_tracknames), , drop=FALSE]
    ans$track <- factor(ans$track, levels=unique(ans$track))
    rownames(ans) <- NULL
    ans
}

### Can be used to quickly check that a combination of genome/tablename
### actually exists.
browseUCSCtrack <- function(genome="hg19",
                            tablename="knownGene",
                            url="http://genome.ucsc.edu/cgi-bin/")
{
    if (!isSingleString(genome))
        stop("'genome' must be a single string")
    if (!isSingleString(tablename))
        stop("'tablename' must be a single string")
    if (!isSingleString(url))
        stop("'url' must be a single string")
    url <- sprintf("%s/hgTrackUi?db=%s&g=%s", url, genome, tablename)
    browseURL(url)
}

.tablename2track <- function(tablename, session)
{
    if (!isSingleString(tablename))
        stop("'tablename' must be a single string")
    supported_tables <- supportedUCSCtables(session)
    idx <- which(supported_tables$tablename == tablename)
    if (length(idx) == 0L)
        stop("UCSC table \"", tablename, "\" is not supported")
    ## Sanity check.
    stopifnot(length(idx) == 1L)  # should never happen
    track <- as.character(supported_tables$track[idx])
    track_tables <- tableNames(ucscTableQuery(session, track=track))
    if (!(tablename %in% track_tables))
        stop("UCSC table \"", tablename, "\" does not exist ",
             "for genome \"", genome(session), "\", sorry")
    track
}

### The table names above (unique key) must be used to name the top-level
### elements of the list below. If no suitable tx_name-to-gene_id mapping is
### available in the UCSC database for a supported table, then there is no
### entry in the list below for this table and makeTxDbFromUCSC() will leave
### the gene table empty.
.UCSC_TXNAME2GENEID_MAPDEFS <- list(
    knownGene=list(
        L2Rchain=list(
            c(tablename="knownToLocusLink",
              Lcolname="name",
              Rcolname="value")
        ),
        gene_id_type="Entrez Gene ID"
    ),
    refGene=list(
        L2Rchain=list(
            c(tablename="hgFixed.refLink",
              Lcolname="mrnaAcc",
              Rcolname="locusLinkId")
        ),
        gene_id_type="Entrez Gene ID"
    ),
    vegaGene=list(
        L2Rchain=list(
            c(tablename="vegaGtp",
              Lcolname="transcript",
              Rcolname="gene")
        ),
        gene_id_type="HAVANA Pseudogene ID"
    ),
    vegaPseudoGene=list(
        L2Rchain=list(
            c(tablename="vegaGtp",
              Lcolname="transcript",
              Rcolname="gene")
        ),
        gene_id_type="HAVANA Pseudogene ID"
    ),
    ## UCSC changed its db schema in September 2011 and the ensGtp table became
    ## unavailable for some genomes (sacCer2, hg18, etc..., apparently for
    ## those assemblies that are not the latest). But the (new?) name2 column
    ## in the ensGene table seems to contain the Ensembl gene ids so the join
    ## with the ensGtp table is not needed anymore.
    #ensGene=list(
    #    L2Rchain=list(
    #        c(tablename="ensGtp",
    #          Lcolname="transcript",
    #          Rcolname="gene")
    #    ),
    #)
    ensGene=c(
        colname="name2",
        gene_id_type="Ensembl gene ID"
    ),
    lincRNAsTranscripts=c(
        colname="name",
        gene_id_type="Name of gene"
    ),
    wgEncodeGencodeManualV3=list(
        L2Rchain=list(
            c(tablename="wgEncodeGencodeClassesV3",
              Lcolname="name",
              Rcolname="geneId")
        ),
        gene_id_type="Ensembl gene ID"
    ),
    wgEncodeGencodeAutoV3=list(
        L2Rchain=list(
            c(tablename="wgEncodeGencodeClassesV3",
              Lcolname="name",
              Rcolname="geneId")
        ),
        gene_id_type="Ensembl gene ID"
    ),
    wgEncodeGencodePolyaV3=list(
        L2Rchain=list(
            c(tablename="wgEncodeGencodeClassesV3",
              Lcolname="name",
              Rcolname="geneId")
        ),
        gene_id_type="Ensembl gene ID"
    ),
    wgEncodeGencodeBasicV17=list(
        L2Rchain=list(
            c(tablename="wgEncodeGencodeAttrsV17",
              Lcolname="transcriptId",
              Rcolname="geneId")
        ),
        gene_id_type="Ensembl gene ID"
    ),
    wgEncodeGencodeCompV17=list(
        L2Rchain=list(
            c(tablename="wgEncodeGencodeAttrsV17",
              Lcolname="transcriptId",
              Rcolname="geneId")
        ),
        gene_id_type="Ensembl gene ID"
    ),
    wgEncodeGencodePseudoGeneV17=list(
        L2Rchain=list(
            c(tablename="wgEncodeGencodeAttrsV17",
              Lcolname="transcriptId",
              Rcolname="geneId")
        ),
        gene_id_type="Ensembl gene ID"
    ),
    wgEncodeGencode2wayConsPseudoV17=list(
        L2Rchain=list(
            c(tablename="wgEncodeGencodeAttrsV17",
              Lcolname="transcriptId",
              Rcolname="geneId")
        ),
        gene_id_type="Ensembl gene ID"
    ),
    wgEncodeGencodePolyaV17=list(
        L2Rchain=list(
            c(tablename="wgEncodeGencodeAttrsV17",
              Lcolname="transcriptId",
              Rcolname="geneId")
        ),
        gene_id_type="Ensembl gene ID"
    ),
    wgEncodeGencodeBasicV14=list(
        L2Rchain=list(
            c(tablename="wgEncodeGencodeAttrsV14",
              Lcolname="transcriptId",
              Rcolname="geneId")
        ),
        gene_id_type="Ensembl gene ID"
    ),
    wgEncodeGencodeCompV14=list(
        L2Rchain=list(
            c(tablename="wgEncodeGencodeAttrsV14",
              Lcolname="transcriptId",
              Rcolname="geneId")
        ),
        gene_id_type="Ensembl gene ID"
    ),
    wgEncodeGencodePseudoGeneV14=list(
        L2Rchain=list(
            c(tablename="wgEncodeGencodeAttrsV14",
              Lcolname="transcriptId",
              Rcolname="geneId")
        ),
        gene_id_type="Ensembl gene ID"
    ),
    wgEncodeGencode2wayConsPseudoV14=list(
        L2Rchain=list(
            c(tablename="wgEncodeGencodeAttrsV14",
              Lcolname="transcriptId",
              Rcolname="geneId")
        ),
        gene_id_type="Ensembl gene ID"
    ),
    wgEncodeGencodePolyaV14=list(
        L2Rchain=list(
            c(tablename="wgEncodeGencodeAttrsV14",
              Lcolname="transcriptId",
              Rcolname="geneId")
        ),
        gene_id_type="Ensembl gene ID"
    ),
    wgEncodeGencodeBasicV7=list(
        L2Rchain=list(
            c(tablename="wgEncodeGencodeAttrsV7",
              Lcolname="transcriptId",
              Rcolname="geneId")
        ),
        gene_id_type="Ensembl gene ID"
    ),
    wgEncodeGencodeCompV7=list(
        L2Rchain=list(
            c(tablename="wgEncodeGencodeAttrsV7",
              Lcolname="transcriptId",
              Rcolname="geneId")
        ),
        gene_id_type="Ensembl gene ID"
    ),
    wgEncodeGencodePseudoGeneV7=list(
        L2Rchain=list(
            c(tablename="wgEncodeGencodeAttrsV7",
              Lcolname="transcriptId",
              Rcolname="geneId")
        ),
        gene_id_type="Ensembl gene ID"
    ),
    wgEncodeGencode2wayConsPseudoV7=list(
        L2Rchain=list(
            c(tablename="wgEncodeGencodeAttrsV7",
              Lcolname="transcriptId",
              Rcolname="geneId")
        ),
        gene_id_type="Ensembl gene ID"
    ),
    wgEncodeGencodePolyaV7=list(
        L2Rchain=list(
            c(tablename="wgEncodeGencodeAttrsV7",
              Lcolname="transcriptId",
              Rcolname="geneId")
        ),
        gene_id_type="Ensembl gene ID"
    ),
    flyBaseGene=list(
        L2Rchain=list(
            c(tablename="flyBaseIsoforms",
              Lcolname="transcript",
              Rcolname="clusterId"),
            c(tablename="flyBaseCanonical",
              Lcolname="clusterId",
              Rcolname="transcript")
        ),
        gene_id_type="Name of canonical transcript in cluster"
    ),
    sgdGene=list(
        L2Rchain=list(
            c(tablename="sgdIsoforms",
              Lcolname="transcript",
              Rcolname="clusterId"),
            c(tablename="sgdCanonical",
              Lcolname="clusterId",
              Rcolname="transcript")
        ),
        gene_id_type="Name of canonical transcript in cluster"
    )
)

.howToGetTxName2GeneIdMapping <- function(tablename)
    .UCSC_TXNAME2GENEID_MAPDEFS[[tablename]]

### The 2 functions below must return a named list with 2 elements:
###   $genes: data.frame with tx_name and gene_id cols;
###   $gene_id_type: single string.
.fetchTxName2GeneIdMappingFromUCSC <- function(session, track,
                                               Ltablename, mapdef)
{
    nlink <- length(mapdef$L2Rchain)
    for (i in seq_len(nlink)) {
        L2Rlink <- mapdef$L2Rchain[[i]]
        tablename <- L2Rlink[["tablename"]]
        Lcolname <- L2Rlink[["Lcolname"]]
        Rcolname <- L2Rlink[["Rcolname"]]
        message("Download the ", tablename, " table ... ", appendLF=FALSE)
        if (tablename == "hgFixed.refLink") {
            query <- ucscTableQuery(session, track, table=tablename)
        } else {
            ## The tables involved in the "left join" don't necessarily belong
            ## to the track of the leftmost table (e.g.
            ## "wgEncodeGencodeAttrsV17" table does NOT belong to the "GENCODE
            ## Genes V17" track).
            #query <- ucscTableQuery(session, track, table=tablename)
            query <- ucscTableQuery(session, table=tablename)
        }
        ucsc_table <- getTable(query)
        message("OK")
        if (!all(has_col(ucsc_table, c(Lcolname, Rcolname))))
            stop("expected cols \"", Lcolname, "\" or/and \"",
                 Rcolname, "\" not found in table ", tablename)
        Lcol <- ucsc_table[[Lcolname]]
        Rcol <- ucsc_table[[Rcolname]]
        if (!is.character(Lcol))
            Lcol <- as.character(Lcol)
        if (!is.character(Rcol))
            Rcol <- as.character(Rcol)
        if (i == 1L) {
            tmp <- data.frame(Lcol=Lcol, Rcol=Rcol)
        } else {
            name2val <- Rcol
            names(name2val) <- Lcol
            tmp <- joinDataFrameWithName2Val(tmp, "Rcol", name2val, "Rcol")
        }
    }
    genes <- data.frame(tx_name=tmp$Lcol, gene_id=tmp$Rcol)
    gene_id_type <- mapdef$gene_id_type
    list(genes=genes, gene_id_type=gene_id_type)
}

.extractTxName2GeneIdMappingFromUCSCTxTable <- function(ucsc_txtable, mapdef)
{
    genes <- data.frame(tx_name=ucsc_txtable[["name"]],
                        gene_id=ucsc_txtable[[mapdef["colname"]]])
    gene_id_type <- mapdef["gene_id_type"]
    list(genes=genes, gene_id_type=gene_id_type)
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Extract the 'transcripts' data frame from UCSC table.
###

.extractTranscriptsFromUCSCTxTable <- function(ucsc_txtable)
{
    message("Extract the 'transcripts' data frame ... ",
            appendLF=FALSE)
    tx_id <- seq_len(nrow(ucsc_txtable))
    tx_name <- ucsc_txtable$name
    tx_chrom <- ucsc_txtable$chrom
    tx_strand <- ucsc_txtable$strand
    tx_start <- ucsc_txtable$txStart + 1L
    tx_end <- ucsc_txtable$txEnd
    transcripts <- data.frame(
        tx_id=tx_id,
        tx_name=tx_name,
        tx_chrom=tx_chrom,
        tx_strand=tx_strand,
        tx_start=tx_start,
        tx_end=tx_end
    )
    message("OK")
    transcripts
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Extract the 'splicings' data frame from UCSC table.
###

### Returns a named list with 2 elements. Each element is itself a list of
### integer vectors with no NAs. The 2 elements have the same "shape".
.extractExonLocsFromUCSCTxTable <- function(ucsc_txtable, check.exonCount=FALSE)
{ 
    exon_count <- ucsc_txtable$exonCount
    if (is.null(exon_count) && check.exonCount)
        stop("UCSC data anomaly: 'ucsc_txtable' has no \"exonCount\" column")
    exon_start <- strsplitAsListOfIntegerVectors(ucsc_txtable$exonStarts)
    exon_end <- strsplitAsListOfIntegerVectors(ucsc_txtable$exonEnds)
    if (is.null(exon_count)) {
        if (!identical(elementNROWS(exon_start),
                       elementNROWS(exon_end)))
            stop("UCSC data anomaly: shape of 'ucsc_txtable$exonStarts' ",
                 "and 'ucsc_txtable$exonEnds' differ")
    } else {
        if (!identical(elementNROWS(exon_start), exon_count))
            stop("UCSC data anomaly: 'ucsc_txtable$exonStarts' ",
                 "inconsistent with 'ucsc_txtable$exonCount'")
        if (!identical(elementNROWS(exon_end), exon_count))
            stop("UCSC data anomaly: 'ucsc_txtable$exonEnds' ",
                 "inconsistent with 'ucsc_txtable$exonCount'")
    }
    list(start=exon_start, end=exon_end)
}

### 'cdsStart0', 'cdsEnd1': single integers (resp. 0-based and 1-based).
### 'exon_start0', 'exon_end1': integer vectors of equal lengths (resp.
### 0-based and 1-based) and with no NAs.
### Returns a list with 2 elements, each of them being an integer vector of
### the same length as 'exon_start0' (or 'exon_end1') that may contain NAs.
### Notes (using genome="hg18"):
###   (1) In refGene table, transcript NM_001146685: cds cumulative length is
###       not a multiple of 3:
###                 name chrom strand txStart   txEnd cdsStart  cdsEnd
###         NM_001146685  chr1      + 1351370 1353029  1351370 1353029
###         exonCount       exonStarts         exonEnds id   name2
###                 2 1351370,1352796, 1351628,1353029,  0 TMEM88B     
###         cdsStartStat cdsEndStat exonFrames
###                 cmpl     incmpl       0,0,
###       --> cds lengths: 1351628 - 1351370 -> 258
###                        1353029 - 1352796 -> 233
###       --> cds cum length: 491
###       Note that the cds end is marked as "incomplete" (see the cdsEndStat
###       col) which, according to UCSC, means that "the CDS is NOT completely
###       contained in the alignment at this end". See this post on the Genome
###       mailing list for more information:
###       https://lists.soe.ucsc.edu/pipermail/genome/2005-December/009184.html
###       Note that the post is about the Gencode Genes. Is it reasonable to
###       assume that this applies to RefSeq Genes too?
###   (2) Same thing in ensGene table, transcript ENST00000371841.
###   (3) All transcripts in the knownGene and ccdsGene tables have a cds
###       cumulative length that is a multiple of 3 (the former doesn't even
###       have the cdsStartStat/cdsStartEnd columns). For hg19 ccdsGene, this
###       is not true anymore :-/
### TODO: Investigate (1) and (2).
.extractUCSCCdsStartEnd <- function(cdsStart0, cdsEnd1,
                                    exon_start0, exon_end1, tx_name)
{
    cds_start0 <- cds_end1 <- integer(length(exon_start0))
    cds_start0[] <- NA_integer_
    cds_end1[] <- NA_integer_
    if (cdsStart0 >= cdsEnd1)
        return(list(cds_start0, cds_end1, FALSE))
    first_exon_with_cds <- which(exon_start0 <= cdsStart0
                                 & cdsStart0 < exon_end1)
    if (length(first_exon_with_cds) != 1L)
        stop("UCSC data ambiguity in transcript ", tx_name,
             ": cannot determine first exon with cds ('cdsStart' ",
             "falls in 0 or more than 1 exon)")
    last_exon_with_cds <- which(exon_start0 < cdsEnd1
                                & cdsEnd1 <= exon_end1)
    if (length(last_exon_with_cds) != 1L)
        stop("UCSC data ambiguity in transcript ", tx_name,
             ": cannot determine last exon with cds ('cdsEnd' ",
             "falls in 0 or more than 1 exon)")
    if (last_exon_with_cds < first_exon_with_cds)
        stop("UCSC data anomaly in transcript ", tx_name,
             ": last exon with cds occurs before first exon with cds")
    exons_with_cds <- first_exon_with_cds:last_exon_with_cds
    cds_start0[exons_with_cds] <- exon_start0[exons_with_cds]
    cds_end1[exons_with_cds] <- exon_end1[exons_with_cds]
    cds_start0[first_exon_with_cds] <- cdsStart0
    cds_end1[last_exon_with_cds] <- cdsEnd1
    ## NAs are OK in here since they indicate the absence of any CDS
    ## (which is common and nothing to write home about)
    ## changed from 50K to 19.3K ...   but the 'bad' ones are not present?
    bad <- sum(cds_end1 - cds_start0, na.rm=TRUE) %% 3L != 0L 
    list(cds_start0, cds_end1, bad) 
}

### 'exon_locs' must be the list of 2 lists returned by
### .extractExonLocsFromUCSCTxTable().
### Returns a named list with 2 elements. Each element is itself a list of
### integer vectors with eventually NAs. The 2 elements have the same
### "shape" as the elements in 'exon_locs' and the NAs occur at the same
### places in the 2 elements.
.extractCdsLocsFromUCSCTxTable <- function(ucsc_txtable, exon_locs)
{
    cdsStart <- ucsc_txtable$cdsStart
    cdsEnd <- ucsc_txtable$cdsEnd
    cds_start <- cds_end <- vector(mode="list", length=nrow(ucsc_txtable))

    startend <- Map(.extractUCSCCdsStartEnd, cdsStart, cdsEnd,
                    exon_locs$start, exon_locs$end, ucsc_txtable$name)

    bad <- sapply(startend, "[[", 3)  
    if (any(bad)) {
        bad_cds <- ucsc_txtable$name[bad]
        msg <- sprintf("UCSC data anomaly in %d transcript(s):
            the cds cumulative length is not a multiple of 3
            for transcripts %s", length(bad_cds),
            paste(sQuote(bad_cds), collapse=" "))
        warning(paste(strwrap(msg, exdent=2), collapse="\n"))
    }

    list(start=sapply(startend, "[[", 1),
         end=sapply(startend, "[[", 2))
}

.extractSplicingsFromUCSCTxTable <- function(ucsc_txtable, transcripts_tx_id)
{
    message("Extract the 'splicings' data frame ... ",
            appendLF=FALSE)
    exon_count <- ucsc_txtable$exonCount
    splicings_tx_id <- rep.int(transcripts_tx_id, exon_count)
    if (length(exon_count) == 0L) {
        exon_rank <- exon_start <- exon_end <-
            cds_start <- cds_end <- integer(0)
    } else {
        if (min(exon_count) <= 0L)
            stop("UCSC data anomaly: 'ucsc_txtable$exonCount' contains ",
                 "non-positive values")
        exon_rank <- makeExonRankCol(exon_count, ucsc_txtable$strand)
        exon_locs <- .extractExonLocsFromUCSCTxTable(ucsc_txtable,
                                                     check.exonCount=TRUE)
        cds_locs <- .extractCdsLocsFromUCSCTxTable(ucsc_txtable, exon_locs)
        exon_start <- unlist(exon_locs$start) + 1L
        exon_end <- unlist(exon_locs$end)
        cds_start <- unlist(cds_locs$start) + 1L
        cds_end <- unlist(cds_locs$end)
    }
    splicings <- data.frame(
        tx_id=splicings_tx_id,
        exon_rank=exon_rank,
        exon_start=exon_start,
        exon_end=exon_end,
        cds_start=cds_start,
        cds_end=cds_end
    )
    message("OK")
    splicings
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Preprocess the 'genes' data frame.
###

.makeUCSCGenes <- function(genes, ucsc_txtable)
{
    #genes <- S4Vectors:::extract_data_frame_rows(genes,
    #                             genes$tx_name %in% ucsc_txtable$name)
    genes
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Download and preprocess the 'chrominfo' data frame.
###

.makeUCSCChrominfo <- function(genome, circ_seqs=DEFAULT_CIRC_SEQS,
        goldenPath_url="http://hgdownload.cse.ucsc.edu/goldenPath")
{
    message("Download and preprocess the 'chrominfo' data frame ... ",
            appendLF=FALSE)
    ucsc_chrominfotable <- GenomeInfoDb:::fetch_ChromInfo_from_UCSC(genome,
                                                goldenPath_url)
    COL2CLASS <- c(
        chrom="character",
        size="integer"
    )
    ucsc_chrominfotable <- setDataFrameColClass(ucsc_chrominfotable, COL2CLASS,
                                                drop.extra.cols=TRUE)
    chrominfo <- data.frame(
        chrom=ucsc_chrominfotable$chrom,
        length=ucsc_chrominfotable$size,
        is_circular=make_circ_flags_from_circ_seqs(ucsc_chrominfotable$chrom,
                                                   circ_seqs)
    )
    oo <- order(rankSeqlevels(chrominfo[ , "chrom"]))
    chrominfo <- S4Vectors:::extract_data_frame_rows(chrominfo, oo)
    message("OK")
    chrominfo
}

## User-friendly wrapper to .makeUCSCChrominfo().
getChromInfoFromUCSC <- function(genome,
          goldenPath_url="http://hgdownload.cse.ucsc.edu/goldenPath")
{
  chrominfo <-.makeUCSCChrominfo(genome, circ_seqs=character(),
                                 goldenPath_url=goldenPath_url)
  chrominfo[ , 1:2]
}



### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Prepare the 'metadata' data frame.
###

.prepareUCSCMetadata <- function(genome, tablename, track, gene_id_type,
                                 full_dataset,
                                 taxonomyId=NA, miRBaseBuild=NA)
{
    message("Prepare the 'metadata' data frame ... ",
            appendLF=FALSE)
    if (!isSingleStringOrNA(miRBaseBuild))
        stop("'miRBaseBuild' must be a a single string or NA")
    if(is.na(taxonomyId)){
        taxonomyId <- GenomeInfoDb:::.taxonomyId(UCSCGenomeToOrganism(genome))
    }else{
        GenomeInfoDb:::.checkForAValidTaxonomyId(taxonomyId)
    }
        
    metadata <- data.frame(
        name=c("Data source", "Genome", "Organism", "Taxonomy ID",
               "UCSC Table", "UCSC Track",
               "Resource URL", "Type of Gene ID",
               "Full dataset",
               "miRBase build ID"),
        value=c("UCSC", genome, UCSCGenomeToOrganism(genome), taxonomyId,
                tablename, track,
                "http://genome.ucsc.edu/", gene_id_type,
                ifelse(full_dataset, "yes", "no"),
                miRBaseBuild)
    )
    message("OK")
    metadata
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### makeTxDbFromUCSC()
###

.makeTxDbFromUCSCTxTable <- function(ucsc_txtable, genes,
        genome, tablename, track, gene_id_type,
        full_dataset,
        circ_seqs=DEFAULT_CIRC_SEQS,
        goldenPath_url="http://hgdownload.cse.ucsc.edu/goldenPath",
        taxonomyId=NA,
        miRBaseBuild=NA)
{
    ucsc_txtable <- setDataFrameColClass(ucsc_txtable, .UCSC_TXCOL2CLASS,
                                         drop.extra.cols=TRUE)

    strand_is_dot <- ucsc_txtable$strand == "."
    if (any(strand_is_dot)) {
        msg <- sprintf("dropped %d transcript(s) for which strand
            was not set (i.e. was set to '.')", sum(strand_is_dot))
        warning(paste(strwrap(msg, exdent=2), collapse="\n"))
        keep_idx <- which(!strand_is_dot)
        ucsc_txtable <- S4Vectors:::extract_data_frame_rows(ucsc_txtable,
                                                            keep_idx)
    }

    transcripts <- .extractTranscriptsFromUCSCTxTable(ucsc_txtable)
    splicings <- .extractSplicingsFromUCSCTxTable(ucsc_txtable,
                                                  transcripts$tx_id)
    genes <- .makeUCSCGenes(genes, ucsc_txtable)
    chrominfo <- .makeUCSCChrominfo(genome, circ_seqs, goldenPath_url)
    metadata <- .prepareUCSCMetadata(genome, tablename, track, gene_id_type,
                                     full_dataset,
                                     taxonomyId,  miRBaseBuild)

    message("Make the TxDb object ... ", appendLF=FALSE)
    txdb <- makeTxDb(transcripts, splicings, genes=genes,
                     chrominfo=chrominfo, metadata=metadata,
                     reassign.ids=TRUE)
    message("OK")
    txdb
}

### The 2 main tasks that makeTxDbFromUCSC() performs are:
###   (1) download the data from UCSC into a data.frame (the getTable() call);
###   (2) store that data.frame in an SQLite db (the
###       .makeTxDbFromUCSCTxTable() call).
### Speed:
###   - for genome="hg18" and tablename="knownGene":
###       (1) download takes about 40-50 sec.
###       (2) db creation takes about 30-35 sec.
makeTxDbFromUCSC <- function(genome="hg19",
        tablename="knownGene",
        transcript_ids=NULL,
        circ_seqs=DEFAULT_CIRC_SEQS,
        url="http://genome.ucsc.edu/cgi-bin/",
        goldenPath_url="http://hgdownload.cse.ucsc.edu/goldenPath",
        taxonomyId=NA,
        miRBaseBuild=NA)
{
    if (!is.null(transcript_ids)) {
        if (!is.character(transcript_ids) || any(is.na(transcript_ids)))
            stop("'transcript_ids' must be a character vector with no NAs")
    }
    if (!isSingleString(url))
        stop("'url' must be a single string")
    if (!isSingleString(goldenPath_url))
        stop("'goldenPath_url' must be a single string")

    ## Create an UCSC Genome Browser session.
    session <- browserSession(url=url)
    genome(session) <- genome
    track <- .tablename2track(tablename, session)

    ## Download the transcript table.
    message("Download the ", tablename, " table ... ", appendLF=FALSE)
    query <- ucscTableQuery(session, track, table=tablename,
                            names=transcript_ids)
    ucsc_txtable <- getTable(query)
    if (ncol(ucsc_txtable) < length(.UCSC_TXCOL2CLASS))
        stop("GenomicFeatures internal error: ", tablename, " table doesn't ",
             "exist, was corrupted during download, or doesn't contain ",
             "transcript information. ",
             "Thank you for reporting this to the GenomicFeatures maintainer ",
             "or to the Bioconductor mailing list, and sorry for the ",
             "inconvenience.")
    message("OK")

    ## Get the tx_name-to-gene_id mapping.
    mapdef <- .howToGetTxName2GeneIdMapping(tablename)
    if (is.null(mapdef)) {
        txname2geneid <- list(genes=NULL, gene_id_type="no gene ids")
    } else if (is.list(mapdef)) {
        txname2geneid <- .fetchTxName2GeneIdMappingFromUCSC(session, track,
                                 tablename, mapdef)
    } else if (is.character(mapdef)) {
        txname2geneid <- .extractTxName2GeneIdMappingFromUCSCTxTable(
                                 ucsc_txtable, mapdef)
    } else {
        stop("GenomicFeatures internal error: invalid 'mapdef'")
    }
    .makeTxDbFromUCSCTxTable(ucsc_txtable, txname2geneid$genes,
                             genome, tablename, track,
                             txname2geneid$gene_id_type,
                             full_dataset=is.null(transcript_ids),
                             circ_seqs=circ_seqs,
                             goldenPath_url=goldenPath_url,
                             taxonomyId=taxonomyId,
                             miRBaseBuild=miRBaseBuild)
}

