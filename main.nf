#!/usr/bin/env nextflow

def helpMessage() {
    log.info """
    Usage:
    The typical command for running the pipeline is as follows:
    nextflow run main.nf --bams sample.bam [Options]
    
    Inputs Options:
    --input         Input file

    Resource Options:
    --max_cpus      Maximum number of CPUs (int)
                    (default: $params.max_cpus)  
    --max_memory    Maximum memory (memory unit)
                    (default: $params.max_memory)
    --max_time      Maximum time (time unit)
                    (default: $params.max_time)
    See here for more info: https://github.com/lifebit-ai/hla/blob/master/docs/usage.md
    """.stripIndent()
}

// Show help message
if (params.help) {
  helpMessage()
  exit 0
}

// Define channels from repository files
projectDir = workflow.projectDir
ch_run_sh_script = Channel.fromPath("${projectDir}/bin/run.sh")

ch_fa = Channel.value(file(params.genome_fasta))
ch_fai = Channel.value(file("${params.genome_fasta}.fai"))
ch_alt = Channel.value(file("${params.genome_fasta}.alt"))

outfmt_uc = params.format.toUpperCase()
outfmt_lc = params.format.toLowerCase()

Channel
  .fromPath(params.input)
  .splitCsv(header:true)
  .map{ row-> tuple(row.rgId, row.rgLb, row.rgPl, row.rgSm, row.rgPu, file(row.read1), file(row.read2)) }
  .set { fqpair_ch }

process bwamem2_postalt {
  input:
    set rgId, rgLb, rgPl, rgSm, rgPu, file(read1), file(read2) from fqpair_ch
    file(fasta) from ch_fa
    file(alt) from ch_alt
    
  output:
    // stdout result
    file("sorted.bam") into sorted_ch

  script:
    def cpus = params.cpus
    def bwamem = "bwa-mem2 mem -K 10000000 -R '@RG\\tID:$rgId\\tLB:$rgLb\\tPL:$rgPl\\tSM:$rgSm\\tPU:$rgPu' -t $cpus $fasta $read1 $read2"
    def postalt = "k8 bwa-postalt.js $alt"
    def fixmate = "samtools fixmate -m --output-fmt bam,level=0 -@ $cpus - -"
    def sort = "samtools sort -m 2G --output-fmt bam,level=1 -T sorttmp -@ $cpus -"

    """
    rm -f sorttmp*
    set -o pipefail
    echo "$bwamem | $postalt | $fixmate | $sort" > sorted.bam
    """
}

sorted_ch = sorted_ch.view { print it }

process samtools_markdup {
  input:
    file '*.bam' from sorted_ch.collect()
  
  output:
    stdout result
  
  script:
    def cpus = params.cpus
    def outidx = outfmt_lc == 'cram' ? 'cram.crai' : 'bam.bai'
    """
    rm -f marked.*
    CLEANED=`ls -m *.bam | tr -d ','`
    echo "samtools markdup --write-index --mode s --output-fmt $outfmt_uc -S --include-fails -T marktmp -@ $cpus -f markdup.met \$CLEANED marked.${outfmt_lc}##idx##marked.${outidx}"
    """
}

result.view { it.trim() }


