/*
 * Authors:
 *      Mamana Mbiyavanga
 *
 *  On behalf of the H3ABionet Consortium
 *  2017
 *
 *
 * Description  : Nextflow pipeline for ...
 *
*/

//---- General definitions --------------------------------------------------//

chromosomes = params.chromosomes.split(',')

println "Project : $workflow.projectDir"
println "Cmd line: $workflow.commandLine"
println "User's populations: ${params.datasets.values().join(", ")}"
println "Chromosomes used: ${chromosomes.join(', ')}"

// Check if files exist
// dbSNP vcf
if(!file(params.dbsnp_vcf).exists()){
    System.err.println "File ${params.dbsnp_vcf} not found. Please check your config file."
    exit 1
}

// Create a channel for initial data which is in chromosomes
vcf_datas = Channel.create()
params.datasets.each { data ->
    if (file(data.value).exists()) {
        vcf_datas << [data.key, file(data.value)]
    }
    else{
        System.err.println "File ${file(data.value)} not found. Please check your config file."
        exit 1
    }
}
vcf_datas.close()

'''
Step 1.1: Split VCF into chrmosomes
'''
process split_vcf_to_chrom {
    tag "split_vcf_to_chrom_${chrm}_${file(vcf_file.baseName).baseName}"
    input:
        set dataset, file(vcf_file) from vcf_datas
        each chrm from chromosomes
    output:
        set dataset, chrm, file(vcf_out) into split_vcf_to_chrom
    script:
        vcf_out = "${file(vcf_file.baseName).baseName}_chr${chrm}.vcf.gz"
        """
        bcftools index --tbi -f ${vcf_file}
        bcftools view \
            --regions ${chrm} \
            ${vcf_file} \
            -Oz -o ${vcf_out}
        """
}

'''
Step 1.2: Annotate dataset using snpEff database
'''
process snpeff_vcf {
    tag "snpeff_${chrm}_${file(vcf_file.baseName).baseName}"
    label "bigmem"
    input:
        set dataset, chrm, file(vcf_file) from split_vcf_to_chrom
    output:
        set dataset, chrm, file("${vcf_out}.gz") into snpeff_vcf
    script:
        vcf_out = "${file(vcf_file.baseName).baseName}_snpeff.vcf"
        """
        snpEff \
            ${params.snpEff_human_db} \
            -lof \
            -stats ${vcf_file.baseName}.html \
            -csvStats ${vcf_file.baseName}.csv \
            -dataDir ${params.snpEff_database} \
            ${vcf_file} > ${vcf_out} -v -nodownload
        bgzip -f ${vcf_out}
        bcftools index --tbi -f ${vcf_out}.gz
        """
}


'''
Step 1.3: Annotate whole dataset with snpEff using dbSNP database
'''
process annotate_dbsnp {
    tag "dbsnp_${chrm}_${file(vcf_file.baseName).baseName}"
    label "bigmem"
    input:
        set dataset, chrm, file(vcf_file) from snpeff_vcf
    output:
        set dataset, chrm, file("${vcf_out}.gz") into annotate_dbsnp
    script:
        vcf_out = "${file(vcf_file.baseName).baseName}_dbsnp.vcf"
        """
        SnpSift \
            annotate \
            -v \
            ${params.dbsnp_vcf} \
            ${vcf_file} > ${vcf_out} 
        bgzip -f ${vcf_out}
        bcftools index --tbi -f ${vcf_out}.gz
        """
}


'''
Step 1.4.1: Convert gwasCat to bed
'''
gwascatalog_cha = Channel.fromPath(params.gwascatalog)
process gwascatalog_to_bed {
    tag "b37tob38_${gwascatalog.baseName}"
    input:
        file(gwascatalog) from gwascatalog_cha
    output:
        set file(gwascatalog), file(outBed) into gwascatalog_to_bed
    script:
        toBed = gwascatalog
        outBed = "${gwascatalog.baseName}.bed"
        template "gwascat2bed.py"
}

'''
Step 1.4.2: LiftOver from build37 to build37 using CrossMap
'''
process bed_b38tob37 {
    tag "b37tob38_${gwascatalog.baseName}"
    label "bigmem"
    time{ 1.hour * task.attempt }
    input:
        set file(gwascatalog), file(inBed) from gwascatalog_to_bed
    output:
        set file(gwascatalog), file(bed_out) into bed_b38tob37
    script:
        file_out = "${gwascatalog.baseName}_b37.tsv"
        bed_out = "${gwascatalog.baseName}_b37.bed"
        """
        CrossMap.py bed \
            ${params.b38tob37_chain} \
            ${gwascatalog.baseName}.bed > \
            ${bed_out}
        """
}


'''
Step 1.4.2: LiftOver from build37 to build37 using CrossMap
'''
process gwascatalog_b37tob38 {
    tag "b37tob38_${gwascat_b38.baseName}"
    label "bigmem"
    time{ 1.hour * task.attempt }
    input:
        set file(gwascat_b38), file(mappedBed) from bed_b38tob37
    output:
        file(gwascat_b37) into gwascatalog_b37tob38
    script:
        gwascat_b37 = "${gwascat_b38.baseName}_b37.tsv"
        template "gwascatb38tob37.py"
}

'''
Step 1.4.2: Annotate whole dataset with snpEff using gwas catalog
'''
annotate_dbsnp_1 = annotate_dbsnp.combine(gwascatalog_b37tob38)
process annotate_gwascat {
    echo true
    tag "gwascat_${chrm}_${file(vcf_file.baseName).baseName}"
    label "bigmem"
    input:
        set dataset, chrm, file(vcf_file), file(gwascat_b37) from annotate_dbsnp_1
    output:
        set dataset, chrm, file("${vcf_out}.gz") into annotate_gwascat
    script:
        vcf_out = "${file(vcf_file.baseName).baseName}_gwascat.vcf"
        """
        SnpSift \
            gwasCat \
            -db ${gwascat_b37} \
            ${vcf_file} \
            -v \
            > ${vcf_out}
        bgzip -f ${vcf_out}
        bcftools index --tbi -f ${vcf_out}.gz
        """
}


'''
Step 1.5: Annotate whole dataset with snpEff using clinvar
'''
process annotate_clinvar {
    tag "clinvar_${chrm}_${file(vcf_file.baseName).baseName}"
    label "bigmem"
    input:
        set dataset, chrm, file(vcf_file) from annotate_gwascat
    output:
        set dataset, chrm, file("${vcf_out}.gz") into annotate_clinvar
    script:
        vcf_out = "${file(vcf_file.baseName).baseName}_clinvar.vcf"
        """
        SnpSift \
            annotate \
            ${params.clinvar} \
            ${vcf_file} \
            -v \
            > ${vcf_out}
        bgzip -f ${vcf_out}
        bcftools index --tbi -f ${vcf_out}.gz
        """
}


'''
Step 1.6: Annotate whole baylor database with snpEff using Cosmic
'''
process annotate_cosmic {
    tag "cosmic_${chrm}_${file(vcf_file.baseName).baseName}"
    label "bigmem"
    input:
        set dataset, chrm, file(vcf_file) from annotate_clinvar
    output:
        set dataset, chrm, file("${vcf_out}.gz") into annotate_cosmic
    script:
        vcf_out = "${file(vcf_file.baseName).baseName}_cosmic.vcf"
        """
        SnpSift \
            annotate \
            ${params.cosmic} \
            ${vcf_file} \
            -v \
            > ${vcf_out}
        bgzip -f ${vcf_out}
        bcftools index --tbi -f ${vcf_out}.gz
        """
}


'''
Step 1.7.1: combine frequency annotation files AGVP, SAHGP, TRYPANOGEN, gnomAD, ExAC per chromosome
'''
mafs_annotations = params.mafs_annotations
mafs_annotations_data = []
mafs_annotations.each{ dataset ->
    chromosomes.each{ chrm ->
        mafs_annotations_data << [dataset.key, chrm, file(dataset.value)]
    }
}
mafs_annotations_data_cha = Channel.from(mafs_annotations_data)
process generate_mafs_chr {
    tag "generate_mafs_${mafs}_${chrm}"
    label "bigmem"
    input:
        set mafs, chrm, file(mafs_file) from mafs_annotations_data_cha
    output:
        set mafs, chrm, file(outTSV) into mafs_annot mode flatten
    script:
        inTSV = mafs_file
        outTSV = "${mafs}_chr${chrm}_mafs.tsv"
        template "mergeMAFannotations.py"
}


'''
Step 1.7.2: Annotate whole dataset  with mafs from AGVP, SAHGP, TRYPANOGEN, gnomAD, ExAC
'''
mafs_annot.into { mafs_annot; mafs_annot_1 }
annotate_cosmic.into { annotate_cosmic; annotate_cosmic_1}
mafs_annot_list = mafs_annot_1.toSortedList().val

def cosmic_annot_list = { dataset, chrm, vcf_file ->
    annotate_cosmic = [:]
    annotate_cosmic[chrm] = [dataset, chrm, file(vcf_file)]
    annot = []
    datasets = []
    mafs_annot_list.each{ mafs, chr, mafs_file ->
        if (chrm == chr){
            annot << file(mafs_file)
            if( mafs.size() >= 6 ) {
                mafs = mafs[0..5]
            }
            datasets << mafs
        }
    }
    annotate_cosmic[chrm] << annot.join(';')
    annotate_cosmic[chrm] << datasets.toSorted().join('-')
    return annotate_cosmic.values()
}
annotate_cosmic_2_cha = annotate_cosmic_1.
        flatMap{ it -> cosmic_annot_list (it) }

process annotate_mafs {
    tag "mafs_${file(vcf_file.baseName).baseName}"
    label "bigmem"
    input:
        set dataset, chrm, file(vcf_file), val(mafs_file), val(mafs) from annotate_cosmic_2_cha
    output:
        set dataset, chrm, file(outVCF) into annotate_mafs
    script:
        outVCF = "${file(vcf_file.baseName).baseName}_mafs-${mafs}.vcf.gz"
        inTSV = mafs_file
        inVCF = vcf_file
        template "annotateVCFwithTSV.py"
}


'''
Step 1.9: Concatenate chromosome VCFs into one
'''
annotate_mafs_list = [:]
annotate_mafs_1 = annotate_mafs.toSortedList().val
annotate_mafs_1.each { dataset, chrm, vcf_chrm_pop ->
    if(dataset in annotate_mafs_list.keySet()) {
        annotate_mafs_list[dataset][1] = annotate_mafs_list[dataset][1] + ' ' + vcf_chrm_pop
    }
    else{
        annotate_mafs_list[dataset] = [dataset, vcf_chrm_pop]
    }
}
annotate_mafs_cha = Channel.from(annotate_mafs_list.values())

process concat_dataset {
    tag "concat_dataset_${dataset}_${chromosomes[0]}-${chromosomes[-1]}"
    label "bigmem"
    publishDir "${params.output_dir}/ANN", overwrite: true, mode:'copy'
    input:
        set dataset, vcf_files from annotate_mafs_cha
    output:
        set dataset, file(vcf_out) into vcf_annotated
    script:
        vcf_out = "${dataset}_annotated_chr${chromosomes[0]}-${chromosomes[-1]}.vcf.gz"
        """
        bcftools concat \
            ${vcf_files} \
            -Oz -o ${dataset}.vcf.gz
        bcftools +fill-tags \
            ${dataset}.vcf.gz \
            -Oz -o ${dataset}.tmp.vcf.gz
        bcftools sort \
            ${dataset}.tmp.vcf.gz \
            -Oz -o ${vcf_out}
        bcftools index --tbi -f ${vcf_out}
        rm -f ${dataset}.tmp* ${dataset}.vcf.gz
        """
}

workflow.onComplete {
    def subject = 'My pipeline execution'
    def recipient = 'mypandos@gmail.com'

    ['mail', '-s', subject, recipient].execute() << """

    Pipeline execution summary
    ---------------------------
    Completed at: ${workflow.complete}
    Duration    : ${workflow.duration}
    Success     : ${workflow.success}
    workDir     : ${workflow.workDir}
    exit status : ${workflow.exitStatus}
    Error report: ${workflow.errorReport ?: '-'}
    """
}