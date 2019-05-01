### =========================================================================
### Utilities for fetching data from UCSC
### -------------------------------------------------------------------------
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

