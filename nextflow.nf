params.referenceGenome = file("/umut/MahmutKarabulut/referans_fasta/fox_reference_genome.fasta")
params.assemblyDir = "/umut/MahmutKarabulut/nextflow_results/assembly_results"
params.mappingDir = "/umut/MahmutKarabulut/nextflow_results/mapping_results"
params.statsDir = "/umut/MahmutKarabulut/nextflow_results/statistics"
params.inputDataDir = "/umut/MahmutKarabulut/input_files"

process spadesAssembly {
    input:
    path dirPath from params.inputDataDir
    file(illuminaFile) from fileTree(dirPath).include("*.fastq")

    output:
    path "${params.assemblyDir}/${illuminaFile.baseName}" into assemblyResults

    script:
    """
    spades.py -o ${params.assemblyDir}/${illuminaFile.baseName} -1 ${illuminaFile} -2 ${illuminaFile}
    """
}



process minimapMapping {
    input:
    path(sampleAssembly) from assemblyResults
    file(referenceGenome)
    output:
    path "${params.mappingDir}/${sampleAssembly.baseName}" into mappingResults
    script:
    """
    minimap2 -ax map-ont ${referenceGenome} ${sampleAssembly}/scaffolds.fasta | samtools view -bS - | samtools sort -o ${params.mappingDir}/${sampleAssembly.baseName}.sorted.bam -
    """
}

process calculateN50N75 {
    input:
    path(mappedBam) from mappingResults
    output:
    path "${params.statsDir}/${mappedBam.baseName}.stats" into n50n75Stats
    script:
    """
    samtools stats ${mappedBam} > ${params.statsDir}/${mappedBam.baseName}.stats
    """
}

process calculateChromosomeCoverage {
    input:
    path(mappedBam) from mappingResults
    output:
    path "${params.statsDir}/${mappedBam.baseName}.chromosome_coverage_summary.txt" into chromosomeCoverage
    script:
    """
    python percentFinder.py
    """
}



workflow {
    // Updated both processes to work with at least 2 files
    spadesAssembly()
    minimapMapping(referenceGenome)
    calculateN50N75()
    calculateChromosomeCoverage()
}



