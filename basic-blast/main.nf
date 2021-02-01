#!/usr/bin/env nextflow
 
/*
 * Defines the pipeline inputs parameters (giving a default value for each for them)
 * Each of the following parameters can be specified as command line options
 */
/*params.query = "$/data/sample.fa"
*params.db = "$baseDir/blast-db/pdb/tiny"
*/

//qry = "mgram.biotch.fa"

//print launchDir + "\n\n\n"

params.query = launchDir +  '/' +  params.qry 
params.db = launchDir + '/' + params.database

 
db_name = file(params.db).name
db_dir = file(params.db).parent 

       

 
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
process blast {
    input:
    path 'query.fa' from fasta_ch
    path db from db_dir
    
    output:

        file 'all_hits' into hits_ch       
 //   file 'top_hits' into hits_ch
 
    """
    blastp -db $db_dir/$db_name -query query.fa -outfmt 6 > all_hits 
    """

}
 
/*
 * Each time a file emitted by the 'top_hits' channel an extract job is executed
 * producing a file containing the matching sequences
 *
*
*process extract {
*    input:
*    path 'top_hits' from hits_ch
*    path db from db_dir
* 
*    output:
*    file 'sequences' into sequences_ch
* 
*    """
*    blastdbcmd -db $db/$db_name -entry_batch top_hits | head -n 10 > sequences
*
*    """
*}
*/
 
/*
 * Collects all the sequences files into a single file
 * and prints the resulting file content when complete
 */

process extract {
    
    input:
        path 'all_hits' from hits_ch
    
    output: 
        file 'sequences' into sequences_ch

        """
        cat all_hits >  sequences 
        """
    
                }


params.fiout = launchDir + "/" + params.out 

sequences_ch
    .collectFile(name: params.fiout)
    //.view { file -> "matching sequences:\n ${file.text}" }





