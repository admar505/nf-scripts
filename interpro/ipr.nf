#! /usr/bin/env nextflow

nextflow.preview.dsl=2

/*
 * Default pipeline parameters. They can be overriden on the command line eg.
 * given `params.foo` specify on the run command line `--foo some_value`.
 */

params.outdir = "results"


params.records_per_file = 100

// blastp parameters

params.interproscan_db = ''

// agat_sp_manage_functional_annotation.pl parameters
//params.id_prefix = 'NBIS'
//params.protein_existence = '5'

log.info("""
NBIS
  _   _ ____ _____  _____
 | \\ | |  _ \\_   _|/ ____|
 |  \\| | |_) || | | (___
 | . ` |  _ < | |  \\___ \\
 | |\\  | |_) || |_ ____) |
 |_| \\_|____/_____|_____/  Annotation Service
 Functional annotation workflow
 ===================================
     outdir                         : ${params.outdir}
 Parallelisation parameters
     records_per_file               : ${params.records_per_file}
 Interproscan parameters
     interproscan_db                : ${params.interproscan_db}
 Merge functional annotation parameters
     id_prefix                      : ${params.id_prefix}
     protein_existence              : ${params.protein_existence}
 """)

workflow {

    main:
        Channel.fromPath(params.gff_annotation, checkIfExists: true)
            .ifEmpty { exit 1, "Cannot find gff file matching ${params.gff_annotation}!\n" }
            .set {annotation}
        Channel.fromPath(params.genome, checkIfExists: true)
            .ifEmpty { exit 1, "Cannot find genome matching ${params.genome}!\n" }
            .set {genome}
        Channel.fromPath("${params.blast_db_fasta}" + (params.blast_db_fasta =~ /^(ht|f)tps?:/ ? '': "*"), checkIfExists: true)
            .ifEmpty { exit 1, "Cannot find blast database files matching ${params.blast_db_fasta}?(.p*)" }
            .set {blastdb}
        functional_annotation(annotation,genome,blastdb)

}

workflow functional_annotation {

    take:
        gff_file
        genome
        blastdb

    main:
        makeblastdb(blastdb,blastdb.filter { it =~ /.p(hr|in|sq)$/ }.ifEmpty('DBFILES_ABSENT'))
        makeblastdb.out.mix(blastdb.filter { !(it =~ /[^.]f(ast|n|)a$/) }).unique().collect().set{blastfiles}
        gff2protein(gff_file,genome.collect())
        blastp(gff2protein.out.splitFasta(by: params.records_per_file, file: true),
            blastfiles)
        interproscan(gff2protein.out.splitFasta(by: params.records_per_file, file: true))
        merge_functional_annotation(gff_file,
            blastp.out.collectFile(name:'blast_merged.tsv').collect(),
            interproscan.out.collectFile(name:'interproscan_merged.tsv').collect(),
            blastdb.collect())

    emit:
        annotation = merge_functional_annotation.out[0]
}

process gff2protein {//ok, so this is emmits the proteins, I _THINK_
                     //grab the protein channel from the blastand attempt to
    label 'protein'  //haxx this bitch  

    input:
    path proteins 

    output:
    path "${gff_file.baseName}_proteins.fasta"

    script://like, split or something? not sure here.
    //"""
    //agat_sp_extract_sequences.pl -o ${gff_file.baseName}_proteins.fasta -f $genome_fasta \\
    //    -p -cfs -cis -ct ${params.codon_table} --gff $gff_file
    //"""
    // agat_sp_extract_sequences.pl is a script from AGAT
    """
    splitfasta.pl <something>
    """

}


process interproscan {

    input:
    path protein_fasta

    output:
    path '*.tsv'

    script:
    applications = params.interproscan_db ? "-appl ${params.interproscan_db}" : ''
    tmpdir = task.scratch ? "-T ${task.scratch}" : ''
  
    """
    export PATH="/opt/interproscan:\$PATH"
    interproscan.sh ${applications} -i $protein_fasta -o ${protein_fasta.baseName}.tsv \\
        -f TSV --iprlookup --goterms -pa -dp -t p ${tmpdir}
    """

}

process merge_functional_annotation {

    publishDir "${params.outdir}/blast_tsv", mode: 'copy', pattern: 'blast_merged.tsv'
    publishDir "${params.outdir}/interproscan_tsv", mode: 'copy', pattern: 'interproscan_merged.tsv'
    publishDir "${params.outdir}/final_annotation", mode: 'copy', pattern: "${gff_annotation.baseName}_plus-functional-annotation.gff"
    label 'AGAT'

    input:
    path gff_annotation
    path merged_blast_results
    path merged_interproscan_results
    path blast_files

    output:
    path "${gff_annotation.baseName}_plus-functional-annotation.gff"
    path "*.tsv", includeInputs:true

    script:
    fasta = blast_files.find { it =~ /\.f(ast|n)?a$/ }
    """
    agat_sp_manage_functional_annotation.pl -f ${gff_annotation} \\
        -b ${merged_blast_results} -i ${merged_interproscan_results} \\
        -db ${params.blast_db_fasta} -id ${params.id_prefix} \\
        -pe ${params.protein_existence} \\
        -o ${gff_annotation.baseName}_plus-functional-annotation.gff
    """
    // agat_sp_manage_functional_annotation.pl is a script from AGAT

}

workflow.onComplete {
    log.info ( workflow.success ? "\nFunctional annotation input preparation complete!\n" : "Oops .. something went wrong\n" )
}

