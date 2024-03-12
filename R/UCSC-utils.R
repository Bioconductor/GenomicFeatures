### =========================================================================
### Utilities for fetching data from UCSC
### -------------------------------------------------------------------------
###


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### list_UCSC_primary_tables()
###

.cached_genome_tracks <- new.env(parent=emptyenv())

.get_genome_tracks <- function(genome, recache=FALSE)
{
    stopifnot(isSingleString(genome), isTRUEorFALSE(recache))
    ans <- .cached_genome_tracks[[genome]]
    if (is.null(ans) || recache) {
        url <- "https://api.genome.ucsc.edu"
        response <- GET(url, path="list/tracks", query=list(genome=genome))
        if (response$status_code != 200L)
            stop(wmsg(genome, ": unknown genome (or ", url, " is down?)"))
        json <- content(response, as="text", encoding="UTF-8")
        ans <- fromJSON(json)[[genome]]
        .cached_genome_tracks[[genome]] <- ans
    }
    ans
}

### Typical usage: list_UCSC_primary_tables("mm9", group="genes")
### Returns a data.frame with 1 row per primary table and 5 columns:
### primary_table, track, type, group, composite_track.
### Note that the "group" and "composite_track" columns can contain NAs.
list_UCSC_primary_tables <- function(genome, group=NULL, recache=FALSE)
{
    stopifnot(is.null(group) || isSingleString(group))
    genome_tracks <- .get_genome_tracks(genome, recache=recache)

    track_groups <- vapply(genome_tracks,
        function(track) {
            group <- track$group
            if (is.null(group)) NA_character_ else group
        },
        character(1), USE.NAMES=FALSE
    )
    if (!is.null(group)) {
        keep_idx <- which(track_groups %in% group)
        genome_tracks <- genome_tracks[keep_idx]
        track_groups <- track_groups[keep_idx]
    }
    track_names <- vapply(genome_tracks,
        function(track) track$shortLabel,
        character(1), USE.NAMES=FALSE
    )
    track_types <- vapply(genome_tracks,
        function(track) track$type,
        character(1), USE.NAMES=FALSE
    )

    ## Extract tracks nested in composite tracks.
    is_composite <- vapply(genome_tracks,
        function(track) identical(track$compositeTrack, "on"),
        logical(1), USE.NAMES=FALSE
    )
    nested_tracks <- lapply(genome_tracks[is_composite],
        function(track) {
            track[vapply(track, is.list, logical(1), USE.NAMES=FALSE)]
        }
    )
    nested_primary_tables <- lapply(nested_tracks, names)
    nested_track_names <- lapply(nested_tracks,
        function(tracks) vapply(tracks,
                                function(track) track$shortLabel,
                                character(1), USE.NAMES=FALSE)
    )
    nested_track_types <- lapply(nested_tracks,
        function(tracks) vapply(tracks,
                                function(track) track$type,
                                character(1), USE.NAMES=FALSE)
    )
    nested_tracks_count <- lengths(nested_tracks)

    ## Sanity checks.
    stopifnot(
        identical(lengths(nested_primary_tables), nested_tracks_count),
        identical(lengths(nested_track_names), nested_tracks_count),
        identical(lengths(nested_track_types), nested_tracks_count)
    )

    ## Prepare columns of final data frame.
    times <- rep.int(1L, length(genome_tracks))
    times[is_composite] <-  nested_tracks_count
    ans_is_composite <- rep.int(is_composite, times)
    ans_primary_table <- rep.int(names(genome_tracks), times)
    ans_primary_table[ans_is_composite] <-
        unlist(nested_primary_tables, use.names=FALSE)
    stopifnot(anyDuplicated(ans_primary_table) == 0L)
    ans_track <- ans_composite_track <- rep.int(track_names, times)
    ans_track[ans_is_composite] <-
        unlist(nested_track_names, use.names=FALSE)
    ans_type <- rep.int(track_types, times)
    ans_type[ans_is_composite] <-
        unlist(nested_track_types, use.names=FALSE)
    ans_group <- rep.int(track_groups, times)
    ans_composite_track[!ans_is_composite] <- NA_character_

    data.frame(
        primary_table=ans_primary_table,
        track=ans_track,
        type=ans_type,
        group=ans_group,
        composite_track=ans_composite_track
    )
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### lookup_organism_by_UCSC_genome()
###

lookup_organism_by_UCSC_genome <- function(genome)
{
    genome <- gsub("\\d+$", "", genome)

    ## Fetch all UCSC genomes with:
    ##   library(RMariaDB)
    ##   dbconn <- dbConnect(MariaDB(), username="genome",
    ##                       host="genome-mysql.soe.ucsc.edu", port=3306)
    ##   genomes <- sort(dbGetQuery(dbconn, "SHOW DATABASES")[[1L]])
    ##   unique(gsub("\\d+$", "", genomes))
    genome2org <- c(
        ailMel="Ailuropoda melanoleuca",
        allMis="Alligator mississippiensis",
        anoCar="Anolis carolinensis",
        anoGam="Anopheles gambiae",
        apiMel="Apis mellifera",
        aplCal="Aplysia californica",
        aptMan="Apteryx australis",
        aquChr="Aquila chrysaetos canadensis",

        balAcu="Balaenoptera acutorostrata",
        bisBis="Bison bison",
        bosTau="Bos taurus",
        braFlo="Branchiostoma floridae",

        caeJap="Caenorhabditis japonica",
        caePb="Caenorhabditis brenneri",
        caeRem="Caenorhabditis remanei",
        calJac="Callithrix jacchus",
        calMil="Callorhinchus milii",
        canFam="Canis familiaris",
        cavPor="Cavia porcellus",
        cb="Caenorhabditis briggsae",
        ce="Caenorhabditis elegans",
        cerSim="Ceratotherium simum",
        chlSab="Chlorocebus sabaeus",
        choHof="Choloepus hoffmanni",
        chrPic="Chrysemys picta",
        ci="Ciona intestinalis",
        criGri= "Cricetulus griseus",

        danRer="Danio rerio",
        dasNov="Dasypus novemcinctus",
        dipOrd="Dipodomys ordii",
        dm="Drosophila melanogaster",
        dp="Drosophila pseudoobscura",
        droAna="Drosophila ananassae",
        droEre="Drosophila erecta",
        droGri="Drosophila grimshawi",
        droMoj="Drosophila mojavensis",
        droPer="Drosophila persimilis",
        droSec="Drosophila sechellia",
        droSim="Drosophila simulans",
        droVir="Drosophila virilis",
        droYak="Drosophila yakuba",

        eboVir="Filoviridae ebolavirus",
        echTel="Echinops telfairi",
        equCab="Equus caballus",
        eriEur= "Erinaceus europaeus",

        felCat="Felis catus",
        fr="Fugu rubripes",

        gadMor="Gadus morhua",
        galGal="Gallus gallus",
        galVar="Galeopterus variegatus",
        gasAcu="Gasterosteus aculeatus",
        geoFor="Geospiza fortis",
        gorGor="Gorilla gorilla",

        hetGla="Heterocephalus glaber",
        hg="Homo sapiens",

        latCha="Latimeria chalumnae",
        loxAfr="Loxodonta africana",

        macEug="Macropus eugenii",
        macFas="Macaca fascicularis",
        manPen="Manis pentadactyla",
        melGal="Meleagris gallopavo",
        melUnd="Melopsittacus undulatus",
        micMur="Microcebus murinus",
        mm="Mus musculus",
        monDom="Monodelphis domestica",
        musFur="Mustela putorius",
        myoLuc="Myotis lucifugus",

        nanPar="Nanorana parkeri",
        nasLar="Nasalis larvatus",
        nomLeu="Nomascus leucogenys",

        ochPri="Ochotona princeps",
        oreNil="Oreochromis niloticus",
        ornAna="Ornithorhynchus anatinus",
        oryCun="Oryctolagus cuniculus",
        oryLat="Oryzias latipes",
        otoGar="Otolemur garnettii",
        oviAri="Ovis aries",

        panPan="Pan paniscus",
        panTro="Pan troglodytes",
        papAnu="Papio anubis",
        papHam="Papio hamadryas",
        petMar="Petromyzon marinus",
        ponAbe="Pongo abelii",
        priPac="Pristionchus pacificus",
        proCap="Procavia capensis",
        pteVam="Pteropus vampyrus",

        rheMac="Macaca mulatta",
        rhiRox="Rhinopithecus roxellana",
        rn="Rattus norvegicus",

        sacCer="Saccharomyces cerevisiae",
        saiBol="Saimiri sciureus",
        sarHar="Sarcophilus harrisii",
        sorAra="Sorex araneus",
        speTri="Spermophilus tridecemlineatus",
        strPur="Strongylocentrotus purpuratus",
        susScr="Sus scrofa",

        taeGut="Taeniopygia guttata",
        tarSyr="Tarsius syrichta",
        tetNig="Tetraodon nigroviridis",
        triMan="Trichechus manatus",
        tupBel="Tupaia belangeri",
        turTru="Tursiops truncatus",

        vicPac="Vicugna pacos",

        xenLae="Xenopus laevis",
        xenTro="Xenopus tropicalis"
    )
    genome2org[genome]
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### UCSC_dbselect()
###

### See https://genome.ucsc.edu/goldenpath/help/mysql.html for how to connect
### to a MySQL server at UCSC.
### Here is an example of how to query the server on the US west coast from
### the Unix command line:
###
###   mysql --user=genome --host=genome-mysql.soe.ucsc.edu mm10 -e "select count(*) from knownToLocusLink;"
###
### By default UCSC_dbselect() uses the server located on the US west coast.
UCSC_dbselect <- function(dbname, from, columns=NULL, where=NULL,
                          server="genome-mysql.soe.ucsc.edu")
{
    columns <- if (is.null(columns)) "*" else paste0(columns, collapse=",")
    SQL <- sprintf("SELECT %s FROM %s", columns, from)
    if (!is.null(where)) {
        stopifnot(isSingleString(where))
        SQL <- paste(SQL, "WHERE", where)
    }
    dbconn <- dbConnect(RMariaDB::MariaDB(), dbname=dbname,
                                             username="genome",
                                             host=server,
                                             port=3306)
    on.exit(dbDisconnect(dbconn))
    dbGetQuery(dbconn, SQL)
}

