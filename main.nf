#!/usr/bin/env nextflow

def helpMessage() {
    log.info """
    Usage:
    The typical command for running the pipeline is as follows:

    nextflow run main.nf --genome_fasta genome.fasta --input data.csv [Options]
    
    Required Arguments:
    --input         Input CSV, see end for format
    --genome_fasta  Reference genome fasta, other files based on dropping extension or appending:
                     - genome.fa - provide path to this file
                     - genome.dict
                     - genome.fa.{0123,alt,amb,ann,pac}
                     - genome.fa.{bwt.2bit.64,pos_packed}
                     - genome.fa.suffixarray_uint64_L{1,2}_PARAMETERS

    Options:
    --outdir        Directory for final results [results]
    --format        Final output as bam or cram [bam]
    --bp_in_batch   Input bases in each batch see bwa/bwa-mem2/bwa-meme -K
                    (for reproducibility) [100000000]

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

faFile = file(params.genome_fasta)
ch_dict = Channel.value(file("${faFile.getParent()}/${faFile.baseName}.dict")) // reheader

ch_fa = Channel.value(file(params.genome_fasta))

ch_0123 = Channel.value(file("${params.genome_fasta}.0123")) // bwa-mem2, bwa-meme
ch_alt = Channel.value(file("${params.genome_fasta}.alt")) // bwa-mem2, bwa-meme
ch_amb = Channel.value(file("${params.genome_fasta}.amb")) // bwa-mem2, bwa-meme
ch_ann = Channel.value(file("${params.genome_fasta}.ann")) // bwa-mem2, bwa-meme
ch_bwt2bit = Channel.value(file("${params.genome_fasta}.bwt.2bit.64")) // bwa-mem2, bwa-meme
ch_pac = Channel.value(file("${params.genome_fasta}.pac")) // bwa-mem2, bwa-meme
ch_pos = Channel.value(file("${params.genome_fasta}.pos_packed")) // bwa-meme
ch_L1 = Channel.value(file("${params.genome_fasta}.suffixarray_uint64_L1_PARAMETERS")) // bwa-meme
ch_L2 = Channel.value(file("${params.genome_fasta}.suffixarray_uint64_L2_PARAMETERS")) // bwa-meme

outfmt_uc = params.format.toUpperCase()
outfmt_lc = params.format.toLowerCase()

Channel
  .fromPath(params.input)
  .splitCsv(header:true)
  .map{ row-> tuple(row.rgId, row.rgLb, row.rgPl, row.rgSm, row.rgPu, file(row.read1), file(row.read2)) }
  .set { fqpair_ch }


process bwamem {
  input:
    set rgId, rgLb, rgPl, rgSm, rgPu, file(read1), file(read2) from fqpair_ch
    file(fasta) from ch_fa
    file(zero123) from ch_0123
    file(alt) from ch_alt
    file(amb) from ch_amb
    file(ann) from ch_ann
    file(bwt2bit) from ch_bwt2bit
    file(pac) from ch_pac
    file(pospack) from ch_pos
    file(l_one) from  ch_L1
    file(l_two) from ch_L2
    file(dict) from ch_dict

  output:
    file 'sorted.bam' into sorted_ch
    file '.command.*'
  
  publishDir path: "${params.outdir}/logs/${task.process}/${task.index}",
    mode: 'copy',
    overwrite: true,
    enabled: true,
    pattern: '.command.*',
    saveAs: { "command.${file(it).getExtension()}" }
  
  script:
    // 'bwa-meme mem' without '-7' behaves like 'bwa-mem2 mem'
    def bwamem = "bwa-meme mem -7 -K ${params.bp_in_batch} -R '@RG\\tID:$rgId\\tLB:$rgLb\\tPL:$rgPl\\tSM:$rgSm\\tPU:$rgPu' -t ${task.cpus} $fasta $read1 $read2"
    def fixmate = "samtools fixmate -m --output-fmt bam,level=0 -@ 1 - -"
    def reheader = "samtools reheader -c 'grep -v ^@SQ > tmp-head && cat $dict tmp-head' -"
    def sort = "samtools sort --output-fmt bam,level=1 -T ./sorttmp -@ ${task.cpus} -"

    def command = "$bwamem | $fixmate | $reheader | $sort > sorted.bam"

    """
    rm -f sorttmp*
    set -o pipefail
    $command
    """
}

process samtools_markdup {
  input:
    // ? is due to odd behaviour in collect when only 1 file received, you get a file of just `.bam`
    file 'sorted_?.bam' from sorted_ch.collect()
    file('genome.fa') from ch_fa
  
  output:
    file 'marked.*' into result

  publishDir path: "${params.outdir}/logs/${task.name}/${task.index}", mode: 'copy', overwrite: true, enabled: true, pattern: '.command*'
  publishDir path: "${params.outdir}", mode: 'move', overwrite: true, enabled: true
  
  script:
    // @TODO merge needs to only be applied if more than one file
    
    def outidx = outfmt_lc == 'cram' ? 'cram.crai' : 'bam.bai'
    def cram_fa = outfmt_lc == 'cram' ? ',reference=genome.fa' : ''
    def merge = "samtools merge -u -@ ${task.cpus} - sorted_*.bam"
    def markdup = "samtools markdup --write-index --mode s --output-fmt ${outfmt_uc}${cram_fa} -S --include-fails -T ./marktmp -@ ${task.cpus} -s -f marked.met \$MARK_IN marked.${outfmt_lc}##idx##marked.${outidx}"
    """
    rm -f marked.* marktmp*
    set -o pipefail
    INPUTS=`ls -1 sorted_*.bam | wc -l`
    if [[ \$INPUTS -eq 1 ]]; then
      MARK_IN='sorted_1.bam'
      $markdup
    else
      MARK_IN='-'
      $merge | $markdup
    fi
    """
}


