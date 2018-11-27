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
        bcftools index --tbi -f ${vcf_file}
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
        bcftools index --tbi -f ${vcf_file}
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
        bcftools index --tbi -f ${vcf_file}
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
        bcftools index --tbi -f ${vcf_file}
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


"""
1.8: Annotate for dbNSFP
"""
annotate_mafs.into{ annotate_mafs; annotate_mafs_1}
process annotate_dbnsfp{
    tag "snpeff_${file(vcf_file.baseName).baseName}"
    label "medium"
    publishDir "${params.work_dir}/data/${dataset}/ALL/VCF", overwrite: true, mode:'symlink'
    input:
        set dataset, chrm, file(vcf_file) from annotate_mafs_1
    output:
        set dataset, chrm, file("${vcf_out}.gz") into annotate_dbnsfp
    script:
        vcf_out = "${file(vcf_file.baseName).baseName}_dbnsfp.vcf"
        """
        SnpSift dbnsfp \
            -f Ancestral_allele,SIFT_score,SIFT_converted_rankscore,SIFT_pred,Uniprot_acc_Polyphen2,Uniprot_id_Polyphen2,Uniprot_aapos_Polyphen2,Polyphen2_HDIV_score,Polyphen2_HDIV_rankscore,Polyphen2_HDIV_pred,Polyphen2_HVAR_score,Polyphen2_HVAR_rankscore,Polyphen2_HVAR_pred,LRT_score,LRT_converted_rankscore,LRT_pred,LRT_Omega,MutationTaster_score,MutationTaster_converted_rankscore,MutationTaster_pred,MutationTaster_model,MutationTaster_AAE,MutationAssessor_UniprotID,MutationAssessor_variant,MutationAssessor_score,MutationAssessor_score_rankscore,MutationAssessor_pred,FATHMM_score,FATHMM_converted_rankscore,FATHMM_pred,PROVEAN_score,PROVEAN_converted_rankscore,PROVEAN_pred,Transcript_id_VEST3,Transcript_var_VEST3,VEST3_score,VEST3_rankscore,MetaSVM_score,MetaSVM_rankscore,MetaSVM_pred,MetaLR_score,MetaLR_rankscore,MetaLR_pred,Reliability_index,M-CAP_score,M-CAP_rankscore,M-CAP_pred,REVEL_score,REVEL_rankscore,MutPred_score,MutPred_rankscore,MutPred_protID,MutPred_AAchange,MutPred_Top5features,CADD_raw,CADD_raw_rankscore,CADD_phred,DANN_score,DANN_rankscore,fathmm-MKL_coding_score,fathmm-MKL_coding_rankscore,fathmm-MKL_coding_pred,fathmm-MKL_coding_group,Eigen_coding_or_noncoding,Eigen-raw,Eigen-phred,Eigen-PC-raw,Eigen-PC-phred,Eigen-PC-raw_rankscore,GenoCanyon_score,GenoCanyon_score_rankscore,integrated_fitCons_score,integrated_fitCons_score_rankscore,integrated_confidence_value,GM12878_fitCons_score,GM12878_fitCons_score_rankscore,GM12878_confidence_value,H1-hESC_fitCons_score,H1-hESC_fitCons_score_rankscore,H1-hESC_confidence_value,HUVEC_fitCons_score,HUVEC_fitCons_score_rankscore,HUVEC_confidence_value,GERP++_NR,GERP++_RS,GERP++_RS_rankscore,phyloP100way_vertebrate,phyloP100way_vertebrate_rankscore,phyloP20way_mammalian,phyloP20way_mammalian_rankscore,phastCons100way_vertebrate,phastCons100way_vertebrate_rankscore,phastCons20way_mammalian,phastCons20way_mammalian_rankscore,SiPhy_29way_pi,SiPhy_29way_logOdds,SiPhy_29way_logOdds_rankscore  \
            -v -db ${params.dbnsfp_db} \
            ${vcf_file} -v > ${vcf_out}
        bgzip -f ${vcf_out}
        bcftools index --tbi -f ${vcf_out}.gz
        """
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
    def recipient = 'mamana.mbiyavanga@uct.ac.za'

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