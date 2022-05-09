#!/usr/bin/env nextflow
 
/*
 * Defines the pipeline inputs parameters (giving a default value for each for them)
 * Each of the following parameters can be specified as command line options
 */
/*params.query = "$/data/sample.fa"
*params.db = "$baseDir/blast-db/pdb/tiny"
*/


 
db_name = file(params.db).name //#these come from the config files 
db_dir = file(params.db).parent 

includeConfig "$baseDir/config/compute_resources.config"

 
/*
 * Given the query parameter creates a channel emitting the query fasta file(s),
 * the file is split in chunks containing as many sequences as defined by the parameter 'chunkSize'.
 * Finally assign the result channel to the variable 'fasta_ch'
 */
Channel
    .fromPath(params.query)
    .splitFasta(by: params.chunkSize, file:true)
    .set { fasta_ch }
 
/*
 * Executes a BLAST job for each chunk emitted by the 'fasta_ch' channel
 * and creates as output a channel named 'top_hits' emitting the resulting
 * BLAST matches 
 */
process ipr {
    input:
    path 'query.fa' from fasta_ch
    path db from db_dir
    
    output:

//        file 'out.tsv'
    file 'top_hits' into hits_ch
 
    """
    bash interproscan.sh  ${db} --input query.fa    --iprlookup --goterms -pa -dp -t p  --formats  tsv   --tempdir /temp 


    """

}
 

process extract {
    
    input:
        path 'all_hits' from hits_ch
    
    output: 
        file 'sequences' into sequences_ch

        """
        cat all_hits > ipr.results.tsv 
        """
    
                }


params.fiout = launchDir + "/" + params.out 

sequences_ch
    .collectFile(name: params.fiout)
    //.view { file -> "matching sequences:\n ${file.text}" }





