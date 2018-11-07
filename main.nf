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

CHRMS = params.chromosomes.split(',')

// All POP
def POPS_ALL = []
params.POPS.each { entry->
    POPS_ALL.addAll(entry.value.split(','))
}
println "Project : $workflow.projectDir"
//println "Git info: $workflow.repository - $workflow.revision [$workflow.commitId]"
println "Cmd line: $workflow.commandLine"
pop_in_sampleFile = file(params.sample_file).readLines().collect{ it.split('\t')[1] }.unique()
println "User's populations: ${POPS_ALL.join(", ")}"
println "Populations in sample file: ${pop_in_sampleFile.join(", ")}"
println "Chromosomes used: ${CHRMS.join(', ')}"


//// Help functions

def split_sample_list_per_pop(POP_, input_file, output_file){
    '''
    Read integrated_call_samples_v3.20130502.ALL.panel and split it by population
    '''
    data = []
    println("Extracting "+POP_+' to '+output_file)
    new File(input_file).eachLine { line, nline ->
        if(nline > 1){
            line = line.trim().split()
            POP = line[1].trim()
            if(POP == POP_){
                data.add(line.join('\t')+'\n')
            }
        }
    }
    data = data.join(' ')
    myFile = file(output_file)
    myFile.text = data
    return data
}

// Create a channel for initial data which is in chromosomes
datas = []
datas_merge = []
CHRMS.each { chromosome ->
    datas << [chromosome, file(sprintf(params.data, chromosome))]
    datas_merge << [chromosome, file(sprintf(params.data_merge, chromosome))]
}
vcf_data = Channel.from(datas)
vcf_data_merge = Channel.from(datas_merge)



'''
Step 1.1: Annotate merge dataset using snpEff database
'''
vcf_data_merge.into{ vcf_data_merge; vcf_data_merge_1}
process annotate_vcf_merge {
    //maxForks 15
    tag "annot_merge_${chrm}_${file(vcf_file.baseName).baseName}"
    memory { 5.GB * task.attempt }
    time{ 6.hour * task.attempt }
    publishDir "${params.work_dir}/VCF_ANN/MERGED/", overwrite: true, mode:'symlink'
    input:
        set val(chrm), file(vcf_file) from vcf_data_merge_1
    output:
        set val(chrm), file("${vcf_out}.gz") into annotate_vcf_merge
    script:
        vcf_out = "${file(vcf_file.baseName).baseName}_snpeff.vcf"
        """
        ${params.snpEff} \
            ${params.snpEff_human_db} \
            -lof \
            -stats ${vcf_file.baseName}.html \
            -csvStats ${vcf_file.baseName}.csv \
            -dataDir ${params.snpEff_database} \
            -c ${params.snpeff_config} \
            ${vcf_file} > ${vcf_out} -v -nodownload
        bgzip -f ${vcf_out}
        bcftools index --tbi -f ${vcf_out}.gz
        """
}


'''
Step 1.2: Annotate whole merge database with snpEff using dbSNP database
'''
annotate_vcf_merge.into { annotate_vcf_merge; annotate_vcf_merge_1}
process annotate_dbsnp_merge {
    echo true
    //maxForks 15
    tag "dbSNP_${chrm}_${file(vcf_file.baseName).baseName}"
    memory { 5.GB * task.attempt }
    time{ 6.hour * task.attempt }
    publishDir "${params.work_dir}/VCF_ANN/MERGED", overwrite: true, mode:'symlink'
    input:
        set val(chrm), file(vcf_file) from annotate_vcf_merge_1
    output:
        set val(chrm), file("${vcf_out}.gz") into annotate_dbsnp_merge
    script:
        vcf_out = "${file(vcf_file.baseName).baseName}_dbsnp.vcf"
        """
        SnpSift \
            annotate \
            ${params.dbsnp_vcf} \
            -c ${params.snpeff_config} \
            ${vcf_file} > ${vcf_out} -v
        bgzip -f ${vcf_out}
        bcftools index --tbi -f ${vcf_out}.gz
        """
}


'''
Step 1.3: Annotate VCF for Ancestral Allele (AA) using in-house python script
'''
annotate_dbsnp_merge.into { annotate_dbsnp_merge; annotate_dbsnp_merge_1 } // Duplicate channel so that it can be used multiple times
process add_ANC_to_VCF_merge {
    echo true
    tag "AA_${chrm}_${file(vcf_file.baseName).baseName}"
    memory { 4.GB * task.attempt }
    time{ 4.hour * task.attempt }
    publishDir "${params.work_dir}/VCF_ANN/MERGED", overwrite: true, mode:'symlink'
    input:
        set val(chrm), file(vcf_file) from annotate_dbsnp_merge_1
    output:
        set val(chrm), file("${vcf_out}.gz") into vcf_anc_merge
    script:
        vcf_out = "${file(vcf_file.baseName).baseName}_anc.vcf"
        """
        gunzip -c ${vcf_file} > ${vcf_file.baseName}
        ${params.python_27_env} ${params.homedir}/templates/add-ANC-to-vcf_new.py -g --in ${vcf_file.baseName} --out ${vcf_out} --genomedata ${params.genomedata_path}
        bgzip -f ${vcf_out}
        bcftools index --tbi -f ${vcf_out}.gz
        rm -f ${vcf_file.baseName}
        """
}


'''
Step 1.4.1: LiftOver from build37 to build37 using picard
'''
gwascatalog_cha = Channel.fromPath(params.gwascatalog)
process gwascatalog_b37tob38 {
    tag "b37tob38_${gwascatalog.baseName}"
    memory { 1.GB * task.attempt }
    time{ 1.hour * task.attempt }
    publishDir "${params.work_dir}/VCF_ANN/", overwrite: true, mode:'symlink'
    input:
        file(gwascatalog) from gwascatalog_cha
    output:
        file(file_out) into gwascatalog_b37tob38
    script:
        file_out = "${gwascatalog.baseName}_b37.tsv"
        bed_out = "${gwascatalog.baseName}_b37.bed"
        """
        ~/miniconda3/envs/ngs_py27/bin/python ${params.homedir}/templates/gwascatb83tob37.py \
            --toBed ${gwascatalog} \
            --outBed ${gwascatalog.baseName}.bed
        ~/miniconda3/envs/ngs_py27/bin/CrossMap.py bed \
            ${params.b37tob38_chain} \
            ${gwascatalog.baseName}.bed > \
            ${bed_out}
        ~/miniconda3/envs/ngs_py27/bin/python ${params.homedir}/templates/gwascatb83tob37.py \
            --mappedBed ${bed_out} \
            --gwascat ${gwascatalog} \
            --outb37 ${file_out}
        """
}


'''
Step 1.4.2: Annotate whole merge database with snpEff using gwas catalog
'''
vcf_anc_merge.into { vcf_anc_merge; vcf_anc_merge_2}
gwascatalog_b37tob38.into{ gwascatalog_b37tob38; gwascatalog_b37tob38_1}
vcf_anc_merge_3 = vcf_anc_merge_2.combine(gwascatalog_b37tob38_1)
process annotate_gwascat_merge {
    echo true
    tag "gwascat_${chrm}_${file(vcf_file.baseName).baseName}"
    memory { 5.GB * task.attempt }
    time{ 6.hour * task.attempt }
    publishDir "${params.work_dir}/VCF_ANN/MERGED", overwrite: true, mode:'symlink'
    input:
        set val(chrm), file(vcf_file), file(gwascat_b37) from vcf_anc_merge_3
    output:
        set val(chrm), file("${vcf_out}.gz") into annotate_gwascat_merge
    script:
        vcf_out = "${file(vcf_file.baseName).baseName}_gwascat.vcf"
        """
        ${params.snpSift} \
            gwasCat \
            -db ${gwascat_b37} \
            ${vcf_file} \
            -v -c ${params.snpeff_config} \
            > ${vcf_out}
        bgzip -f ${vcf_out}
        bcftools index --tbi -f ${vcf_out}.gz
        """
}


'''
Step 1.5: Annotate whole merge database with snpEff using clinvar
'''
annotate_gwascat_merge.into { annotate_gwascat_merge; annotate_gwascat_merge_2}
process annotate_clinvar_merge {
    echo true
    tag "clinvar_${chrm}_${file(vcf_file.baseName).baseName}"
    memory { 5.GB * task.attempt }
    time{ 6.hour * task.attempt }
    publishDir "${params.work_dir}/VCF_ANN/MERGED", overwrite: true, mode:'symlink'
    input:
        set val(chrm), file(vcf_file) from annotate_gwascat_merge_2
    output:
        set val(chrm), file("${vcf_out}.gz") into annotate_clinvar_merge
    script:
        vcf_out = "${file(vcf_file.baseName).baseName}_clinvar.vcf"
        """
        ${params.snpSift} \
            annotate \
            ${params.clinvar} \
            ${vcf_file} \
            -v -c ${params.snpeff_config} \
            > ${vcf_out}
        bgzip -f ${vcf_out}
        bcftools index --tbi -f ${vcf_out}.gz
        """
}


'''
Step 1.6: Annotate whole baylor database with snpEff using clinvar
'''
annotate_clinvar_merge.into { annotate_clinvar_merge; annotate_clinvar_merge_2}
process annotate_cosmic_merge {
    tag "cosmic_${chrm}_${file(vcf_file.baseName).baseName}"
    memory { 5.GB * task.attempt }
    time{ 6.hour * task.attempt }
    publishDir "${params.work_dir}/VCF_ANN/MERGED", overwrite: true, mode:'symlink'
    input:
        set val(chrm), file(vcf_file) from annotate_clinvar_merge_2
    output:
        set val(chrm), file("${vcf_out}.gz") into annotate_cosmic_merge
    script:
        vcf_out = "${file(vcf_file.baseName).baseName}_cosmic.vcf"
        """
        ${params.snpSift} \
            annotate \
            ${params.cosmic} \
            ${vcf_file} \
            -v -c ${params.snpeff_config} \
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
    CHRMS.each{ chrm ->
        mafs_annotations_data << [dataset.key, chrm, file(dataset.value)]
    }
}
mafs_annotations_data_cha = Channel.from(mafs_annotations_data)
process mafs_annot_chrs {
    tag "mafs_annot_${mafs_dataset}_${chrm}"
    memory { 5.GB * task.attempt }
    time { 3.hour * task.attempt }
    publishDir "${params.ref_dir}/pop_mafs/${mafs_dataset}", overwrite: true, mode:'symlink'
    input:
        set val(mafs_dataset), val(chrm), file(mafs_file) from mafs_annotations_data_cha
    output:
        set val(mafs_dataset), val(chrm), file(tsv_out) into mafs_annot_merge mode flatten
    script:
        tsv_out = "${mafs_dataset}_chr${chrm}_mafs.tsv"
        """
        ${params.homedir}/templates/annotateVCFwithTSV.py \
            --inTSV ${mafs_file} \
            --chrm ${chrm} \
            --outTSV ${tsv_out}
        """
}


'''
Step 1.7.2: Annotate whole merge database with mafs from AGVP, SAHGP, TRYPANOGEN, gnomAD, ExAC
'''
mafs_annot_merge.into { mafs_annot_merge; mafs_annot_merge_1 }
annotate_cosmic_merge.into { annotate_cosmic_merge; annotate_cosmic_merge_1}
mafs_annot_merge_list = mafs_annot_merge_1.toSortedList().val

def cosmic_annot_list_merge = {
    annotate_cosmic_merge = [:]
    chrm        = it[0]
    vcf_file    = it[1]
    annotate_cosmic_merge[chrm] = [chrm, file(vcf_file)]
    annot = []
    datasets = []
    mafs_annot_merge_list.each{ mafs_dataset, chr, mafs_file ->
        if (chrm == chr){
            annot << file(mafs_file)
            if( mafs_dataset.size() >= 6 ) {
                mafs_dataset = mafs_dataset[0..5]
            }
            datasets << mafs_dataset
        }
    }
    annotate_cosmic_merge[chrm] << annot.join(';')
    annotate_cosmic_merge[chrm] << datasets.toSorted().join('-')
    return annotate_cosmic_merge.values()
}
annotate_cosmic_merge_2_cha = annotate_cosmic_merge_1.
        flatMap{ it -> cosmic_annot_list_merge (it) }

process annotate_mafs_merge {
    tag "mafs_${file(vcf_file.baseName).baseName}"
    memory { 10.GB * task.attempt }
    time { 10.hour * task.attempt }
    maxRetries 5
    //maxForks 15
    publishDir "${params.work_dir}/VCF_ANN/MERGED", overwrite: true, mode:'symlink'
    input:
        set val(chrm), file(vcf_file), val(mafs_files), val(mafs_dataset) from annotate_cosmic_merge_2_cha
    output:
        set val(chrm), file("${vcf_out}.gz") into annotate_mafs_merge
    script:
        vcf_out = "${file(vcf_file.baseName).baseName}_mafs-${mafs_dataset}.vcf"
        """
        gunzip -c ${vcf_file} > ${vcf_file.baseName}
        ${params.homedir}/templates/annotateVCFwithTSV.py \
            --inVCF ${vcf_file.baseName} \
            --inTSV \'${mafs_files}\' \
            --outVCF ${vcf_out}
        ${params.venv}/bin/bgzip -f ${vcf_out}
        ${params.venv}/bin/bcftools index --tbi -f ${vcf_out}.gz
        rm -f ${vcf_file.baseName}
        
        """
}


'''
Step 1.8: Filter out related individuals and remove sites with high missingness
'''
annotate_mafs_merge.into { annotate_mafs_merge; annotate_mafs_merge_1}
process norelated_high_miss_merge {
    tag "norelated_${chrm}_${file(vcf_file.baseName).baseName}"
    memory { 10.GB * task.attempt }
    time { 5.hour * task.attempt }
    //maxForks 15
    publishDir "${params.work_dir}/VCF_FILTERED/MERGED", overwrite: true, mode:'symlink'
    input:
        set val(chrm), file(vcf_file) from annotate_mafs_merge_1
    output:
        set val(chrm), file(vcf_out) into norelated_merge
    script:
        vcf_out = "${file(vcf_file.baseName).baseName}_noRelated_noHighMiss.vcf.gz"
        """
        ${params.venv}/bin/vcftools \
            --gzvcf ${vcf_file} \
            --remove ${params.sample_to_remove_merge} \
            --exclude-positions ${params.high_missing_snp} \
            --recode --recode-INFO-all -c | \
        ${params.venv}/bin/bgzip -c > ${vcf_out}
        ${params.venv}/bin/bcftools index --tbi -f ${vcf_out}
        """
}



'''
Step 1.9: Filter out sites with ALT='.'
'''
norelated_merge.into { norelated_merge; norelated_merge_1}
process remove_noALT_merge {
    tag "norelated_${chrm}_${file(vcf_file.baseName).baseName}"
    memory { 10.GB * task.attempt }
    time { 5.hour * task.attempt }
    //maxForks 15
    publishDir "${params.work_dir}/VCF_FILTERED/MERGED", overwrite: true, mode:'symlink'
    input:
        set val(chrm), file(vcf_file) from norelated_merge_1
    output:
        set val(chrm), file(snp_out), file(vcf_out) into remove_noALT_merge
    script:
        vcf_out = "${file(vcf_file.baseName).baseName}_noALT.vcf.gz"
        snp_out = "${file(vcf_file.baseName).baseName}_noALT.snp"
        """
        ${params.venv}/bin/bcftools \
            view -i 'ALT="."' ${vcf_file} | \
        ${params.venv}/bin/bcftools \
            query -f '%CHROM  %POS  %REF  %ALT\n' \
            > ${snp_out}
        ${params.venv}/bin/bcftools \
            view -e 'ALT="."' ${vcf_file} \
            -Oz -o ${vcf_out}
        ${params.venv}/bin/bcftools index --tbi -f ${vcf_out}
        """
}


'''
Step 1.8: Sites only
'''
remove_noALT_merge.into { remove_noALT_merge; remove_noALT_merge_1}
process vcf_sites_only_merge {
    tag "Sites_only_${file(vcf_file.baseName).baseName}"
    memory { 10.GB * task.attempt }
    time { 5.hour * task.attempt }
    publishDir "${params.work_dir}/SITE_ONLY/MERGED", overwrite: true, mode:'copy'
    input:
        set val(chrm), file(snp_file), file(vcf_file) from remove_noALT_merge_1
    output:
        set val(chrm), file(vcf_out) into vcf_sites_only_merge
    script:
        vcf_out = "${file(vcf_file.baseName).baseName}_sitesOnly.vcf.gz"
        """
        bcftools \
            view ${vcf_file} \
            --drop-genotypes \
            -Oz -o ${vcf_out}
        bcftools index --tbi -f ${vcf_out}
        """
}


'''
Step 1.9: Concatenate chromosome VCFs into one
'''
vcf_sites_only_merge.into { vcf_sites_only_merge; vcf_sites_only_merge_1}
vcf_sites_only_merge_list = []
vcf_sites_only_merge_1.toSortedList().val.each { chrm, vcf_chrm_pop ->
    vcf_sites_only_merge_list << vcf_chrm_pop
}
vcf_sites_only_merge_cha = Channel.from(vcf_sites_only_merge_list.join(' '))
process merge_vcf_sites_only {
    tag { "merge_vcf_sites_only" }
    memory { 40.GB * task.attempt }
    time{ 2.hour * task.attempt }
    publishDir "${params.work_dir}/SITE_ONLY/MERGED", overwrite: true, mode:'symlink'
    input:
        val datas from vcf_sites_only_merge_cha
    output:
        file(vcf_out) into merge_vcf_sites_only
    script:
        vcf_out = "Eagle.merged_${CHRMS[0]}-${CHRMS[-1]}.vcf.gz"
        """
        bcftools concat ${datas} -Oz -o ${vcf_out}
        """
}


'''
Step 1.10: Annotate merge sites-only dataset using snpEff database
'''
merge_vcf_sites_only.into{ merge_vcf_sites_only; merge_vcf_sites_only_1}
process annotate_vcf_merge_sitesOnly {
    //maxForks 15
    tag "annot_sitesOnly_merge"
    memory { 25.GB * task.attempt }
    time{ 6.hour * task.attempt }
    publishDir "${params.work_dir}/SITE_ONLY/MERGED/", overwrite: true, mode:'symlink'
    input:
        file(vcf_file) from merge_vcf_sites_only_1
    output:
        file("${vcf_out}.gz") into annotate_vcf_merge_sitesOnly
    script:
        vcf_out = "${file(vcf_file.baseName).baseName}_snpeff.vcf"
        """
        ${params.snpEff} \
            ${params.snpEff_human_db} \
            -lof \
            -stats ${vcf_file.baseName}.html \
            -csvStats ${vcf_file.baseName}.csv \
            -dataDir ${params.snpEff_database} \
            -c ${params.snpeff_config} \
            ${vcf_file} > ${vcf_out} -v -nodownload
        bgzip -f ${vcf_out}
        bcftools index --tbi -f ${vcf_out}.gz
        """
}


// BAYLOR dataset
'''
Step 2.1: Annotate whole baylor using snpEff database
'''
vcf_data.into{ vcf_data; vcf_data_2}
process annotate_vcf {
    tag "annot_${chrm}_${file(vcf_file.baseName).baseName}"
    //maxForks 15
    memory { 5.GB * task.attempt }
    time{ 6.hour * task.attempt }
    publishDir "${params.work_dir}/VCF_ANN/BAYLOR/", overwrite: true, mode:'symlink'
    input:
        set val(chrm), file(vcf_file) from vcf_data_2
    output:
        set val(chrm), file("${vcf_out}.gz") into main
    script:
        vcf_out = "${file(vcf_file.baseName).baseName}_snpeff.vcf"
        """
        ${params.snpEff} \
            ${params.snpEff_human_db} \
            -lof \
            -stats ${vcf_file.baseName}.html \
            -csvStats ${vcf_file.baseName}.csv \
            -dataDir ${params.snpEff_database} \
            -c ${params.snpeff_config} \
            ${vcf_file} > ${vcf_out} -v -nodownload
        bgzip -f ${vcf_out}
        bcftools index --tbi -f ${vcf_out}.gz
        """
}
main.into {main; annotate_vcf_sub}
annotate_vcf_sub.subscribe{
    println "Finished ${it.join(', ')}"
}


'''
Step 2.2: Annotate whole baylor database with snpEff using dbSNP database
'''
main.into { main; annotate_vcf_1}
process annotate_dbsnp_baylor {
    echo true
    //maxForks 15
    tag "dbSNP_${chrm}_${file(vcf_file.baseName).baseName}"
    memory { 5.GB * task.attempt }
    time{ 6.hour * task.attempt }
    publishDir "${params.work_dir}/VCF_ANN/BAYLOR", overwrite: true, mode:'symlink'
    input:
        set val(chrm), file(vcf_file) from annotate_vcf_1
    output:
        set val(chrm), file("${vcf_out}.gz") into annotate_dbsnp_baylor
    script:
        vcf_out = "${file(vcf_file.baseName).baseName}_dbsnp.vcf"
        """
        SnpSift \
            annotate \
            ${params.dbsnp_vcf} \
            -c ${params.snpeff_config} \
            ${vcf_file} > ${vcf_out} -v
        bgzip -f ${vcf_out}
        bcftools index --tbi -f ${vcf_out}.gz
        """
}


'''
Step 2.3: Annotate VCF for Ancestral Allele (AA) using in-house python script
'''
annotate_dbsnp_baylor.into { annotate_dbsnp_baylor; annotate_dbsnp_baylor_1 } // Duplicate channel so that it can be used multiple times
process add_ANC_to_VCF_baylor {
    tag "AA_${chrm}_${file(vcf_file.baseName).baseName}"
    memory { 4.GB * task.attempt }
    time{ 8.hour * task.attempt }
    publishDir "${params.work_dir}/VCF_ANN/BAYLOR", overwrite: true, mode:'symlink'
    input:
        set val(chrm), file(vcf_file) from annotate_dbsnp_baylor_1
    output:
        set val(chrm), file("${vcf_out}.gz") into vcf_anc_baylor
    script:
        vcf_out = "${file(vcf_file.baseName).baseName}_anc.vcf"
        """
        gunzip -c ${vcf_file} > ${vcf_file.baseName}
        ${params.python_27_env} ${params.homedir}/templates/add-ANC-to-vcf_new.py -g --in ${vcf_file.baseName} --out ${vcf_out} --genomedata ${params.genomedata_path}
        bgzip -f ${vcf_out}
        bcftools index --tbi -f ${vcf_out}.gz
        rm -f ${vcf_file.baseName}
        """
}
vcf_anc_baylor.into {vcf_anc_baylor; vcf_anc_baylor_sub}
vcf_anc_baylor_sub.subscribe{
    println "Finished ${it.join(', ')}"
}


'''
Step 2.4.1: Annotate whole merge database with snpEff using gwas catalog
'''
vcf_anc_baylor.into { vcf_anc_baylor; vcf_anc_baylor_2}
gwascatalog_b37tob38.into{ gwascatalog_b37tob38; gwascatalog_b37tob38_2}
vcf_anc_baylor_3 = vcf_anc_baylor_2.combine(gwascatalog_b37tob38_2)
process annotate_gwascat_baylor {
    tag "gwascat_${chrm}_${file(vcf_file.baseName).baseName}"
    memory { 5.GB * task.attempt }
    time{ 6.hour * task.attempt }
    publishDir "${params.work_dir}/VCF_ANN/BAYLOR", overwrite: true, mode:'symlink'
    input:
        set val(chrm), file(vcf_file), file(gwascat_b37) from vcf_anc_baylor_3
    output:
        set val(chrm), file("${vcf_out}.gz") into annotate_gwascat_baylor
    script:
        vcf_out = "${file(vcf_file.baseName).baseName}_gwascat.vcf"
        """
        ${params.snpSift} \
            gwasCat \
            -db ${gwascat_b37} \
            ${vcf_file} \
            -v -c ${params.snpeff_config} \
            > ${vcf_out}
        bgzip -f ${vcf_out}
        bcftools index --tbi -f ${vcf_out}.gz
        """
}


'''
Step 2.5: Annotate whole baylor database with snpEff using clinvar
'''
annotate_gwascat_baylor.into { annotate_gwascat_baylor; annotate_gwascat_baylor_2}
process annotate_clinvar_baylor {
    tag "clinvar_${chrm}_${file(vcf_file.baseName).baseName}"
    memory { 5.GB * task.attempt }
    time{ 6.hour * task.attempt }
    publishDir "${params.work_dir}/VCF_ANN/BAYLOR", overwrite: true, mode:'symlink'
    input:
        set val(chrm), file(vcf_file) from annotate_gwascat_baylor_2
    output:
        set val(chrm), file("${vcf_out}.gz") into annotate_clinvar_baylor
    script:
        vcf_out = "${file(vcf_file.baseName).baseName}_clinvar.vcf"
        """
        ${params.snpSift} \
            annotate \
            ${params.clinvar} \
            ${vcf_file} \
            -v -c ${params.snpeff_config} \
            > ${vcf_out}
        bgzip -f ${vcf_out}
        bcftools index --tbi -f ${vcf_out}.gz
        """
}


'''
Step 2.6.2: Annotate whole baylor database with snpEff using cosmic
'''
annotate_clinvar_baylor.into { annotate_clinvar_baylor; annotate_clinvar_baylor_2}
process annotate_cosmic_baylor {
    tag "cosmic_${chrm}_${file(vcf_file.baseName).baseName}"
    memory { 5.GB * task.attempt }
    time { 6.hour * task.attempt }
    publishDir "${params.work_dir}/VCF_ANN/BAYLOR", overwrite: true, mode:'symlink'
    input:
        set val(chrm), file(vcf_file) from annotate_clinvar_baylor_2
    output:
        set val(chrm), file("${vcf_out}.gz") into annotate_cosmic_baylor
    script:
        vcf_out = "${file(vcf_file.baseName).baseName}_cosmic.vcf"
        """
        ${params.snpSift} \
            annotate \
            ${params.cosmic} \
            ${vcf_file} \
            -v -c ${params.snpeff_config} \
            > ${vcf_out}
        bgzip -f ${vcf_out}
        bcftools index --tbi -f ${vcf_out}.gz
        """
}


'''
Step 2.7.1: combine frequency annotation files AGVP, SAHGP, TRYPANOGEN, gnomAD, ExAC per chromosome
'''

'''
Step 2.7.2: Annotate whole baylor database with mafs from AGVP, SAHGP, TRYPANOGEN, gnomAD, ExAC
'''

annotate_cosmic_baylor.into { annotate_cosmic_baylor; annotate_cosmic_baylor_1}

def cosmic_annot_list = {
    annotate_cosmic_baylor_2 = [:]
    chrm        = it[0]
    vcf_file    = it[1]
    annotate_cosmic_baylor_2[chrm] = [chrm, file(vcf_file)]
    annot = []
    datasets = []
    mafs_annot_merge_list.each{ mafs_dataset, chr, mafs_file ->
        if (chrm == chr){
            annot << file(mafs_file)
            if( mafs_dataset.size() >= 6 ) {
                mafs_dataset = mafs_dataset[0..5]
            }
            datasets << mafs_dataset
        }
    }
    annotate_cosmic_baylor_2[chrm] << annot.join(';')
    annotate_cosmic_baylor_2[chrm] << datasets.toSorted().join('-')
    return annotate_cosmic_baylor_2.values()
}
annotate_cosmic_baylor_2_cha = annotate_cosmic_baylor_1.
        flatMap{ it -> cosmic_annot_list (it) }

process annotate_mafs_baylor {
    tag "mafs_${file(vcf_file.baseName).baseName}"
    memory { 10.GB * task.attempt }
    time { 10.hour * task.attempt }
    maxRetries 5
    //maxForks 15
    publishDir "${params.work_dir}/VCF_ANN/BAYLOR", overwrite: true, mode:'symlink'
    input:
        set val(chrm), file(vcf_file), val(mafs_files), val(mafs_dataset) from annotate_cosmic_baylor_2_cha
    output:
        set val(chrm), file("${vcf_out}.gz") into annotate_mafs_baylor
    script:
       vcf_out = "${file(vcf_file.baseName).baseName}_mafs-${mafs_dataset}.vcf"
        """
        gunzip -c ${vcf_file} > ${vcf_file.baseName}
        ${params.homedir}/templates/annotateVCFwithTSV.py \
            --inVCF ${vcf_file.baseName} \
            --inTSV \'${mafs_files}\' \
            --outVCF ${vcf_out}
        ${params.venv}/bin/bgzip -f ${vcf_out}
        ${params.venv}/bin/bcftools index --tbi -f ${vcf_out}.gz
        rm -f ${vcf_file.baseName}
        """
}


'''
Step 2.8: Filter out related individuals and remove sites with high missingness
'''
annotate_mafs_baylor.into { annotate_mafs_baylor; annotate_mafs_baylor_1}
process norelated_high_miss_baylor {
    tag "norelated_${chrm}_${file(vcf_file.baseName).baseName}"
    memory { 10.GB * task.attempt }
    time { 5.hour * task.attempt }
    //maxForks 15
    publishDir "${params.work_dir}/VCF_FILTERED/BAYLOR", overwrite: true, mode:'symlink'
    input:
        set val(chrm), file(vcf_file) from annotate_mafs_baylor_1
    output:
        set val(chrm), file(vcf_out), file("${vcf_out}.tbi") into norelated_baylor
    script:
        vcf_out = "${file(vcf_file.baseName).baseName}_noRelated_noHighMiss.vcf.gz"
        """
        vcftools \
            --gzvcf ${vcf_file} \
            --remove ${params.related_sample} \
            --exclude-positions ${params.high_missing_snp} \
            --recode --recode-INFO-all -c | \
        bgzip -c > ${vcf_out}
        bcftools index --tbi -f ${vcf_out}
        """
}


'''
Step 2.9: Filter VCF by coverage depth using bcftools
'''
norelated_baylor.into{ norelated_baylor; norelated_baylor_1}
process filter_vcf_by_depth_baylor {
    tag "dp6_${chrm}_${file(vcf_file.baseName).baseName}"
    memory { 10.GB * task.attempt }
    time{ 2.hour * task.attempt }
    //maxForks 15
    publishDir "${params.work_dir}/VCF_FILTERED/BAYLOR", overwrite: true, mode:'symlink'
    input:
        set val(chrm), file(vcf_file), file(vcf_file_tbi) from norelated_baylor_1
    output:
        set val(chrm), file(vcf_out), file("${vcf_out}.tbi") into vcf_6depth_baylor
    script:
        vcf_out = "${file(vcf_file.baseName).baseName}_dp6.vcf.gz"
        """
        bcftools view -i 'DP>6' ${vcf_file} | bgzip -c > ${vcf_out}
        bcftools index --tbi -f ${vcf_out}
        """
}


'''
Step 2.10: Filter sites with AA
'''
vcf_6depth_baylor.into { vcf_6depth_baylor; vcf_6depth_baylor_1}
process vcf_anc_only {
    echo true
    tag "AA_only_${chrm}_${file(vcf_file.baseName).baseName}"
    memory { 5.GB * task.attempt }
    time{ 2.hour * task.attempt }
    publishDir "${params.work_dir}/VCF_FILTERED/BAYLOR", overwrite: true, mode:'symlink'
    input:
        set val(chrm), file(vcf_file), file(vcf_file_tbi) from vcf_6depth_baylor_1
    output:
        set val(chrm), file(vcf_out), file("${vcf_out}.tbi") into vcf_anc_only
    script:
        vcf_out = "${file(vcf_file.baseName).baseName}_ancOnly.vcf.gz"
        """
        bcftools view -i 'AA!="." & AA!="-" & AA!="N"' ${vcf_file} | bgzip -c > ${vcf_out}
        bcftools index --tbi -f ${vcf_out}
        """
}


'''
Step 2.11: Filter out multi-allelic sites
'''
vcf_anc_only.into { vcf_anc_only; vcf_anc_only_1}
process biallOnly_baylor {
    tag "biallOnly_${chrm}_${file(vcf_file.baseName).baseName}"
    memory { 5.GB * task.attempt }
    time { 5.hour * task.attempt }
    //maxForks 15
    publishDir "${params.work_dir}/VCF_FILTERED/BAYLOR", overwrite: true, mode:'symlink'
    input:
        set val(chrm), file(vcf_file), file(vcf_file_tbi) from vcf_anc_only_1
    output:
        set val(chrm), file(vcf_out), file("${vcf_out}.tbi") into biallOnly_baylor
    script:
        vcf_out = "${file(vcf_file.baseName).baseName}_biall.vcf.gz"
        """
        vcftools \
            --gzvcf ${vcf_file} \
            --exclude-positions ${params.biall_only} \
            --recode --recode-INFO-all -c | \
        bgzip -c > ${vcf_out}
        bcftools index --tbi -f ${vcf_out}
        """
}


POPS_ALL1 = Channel.from(POPS_ALL).combine([file(params.sample_file)])
'''
Step 2.12: Generate sample file
'''
process split_POP_samples {
    echo true
    tag { "split_POP_samples_${POP}" }
    memory { 2.GB  * task.attempt }
    time{ 1.hour * task.attempt }
    publishDir "${params.work_dir}/samples/", overwrite: true, mode:'symlink'
    input:
        set val(POP), file(sample_file) from POPS_ALL1
    output:
        set val(POP), file(sample_out) into split_POP_samples_all
    script:
        sample_out = "${POP}.sample"
        """
        grep ${POP} ${sample_file} | cut -f1- > ${sample_out}
        """
}

'''
Step 2.12.1: Generate reduced sample file with 26 individuals each
'''
split_POP_samples_all.into { split_POP_samples_all; split_POP_samples_2 }
process reduce_POP_samples {
    echo true
    tag { "reduce_samples_${sample_size}_${POP}" }
    memory { 2.GB  * task.attempt }
    time{ 1.hour * task.attempt }
    publishDir "${params.work_dir}/${sample_size}/samples/", overwrite: true, mode:'copy'
    input:
        set val(POP), file(sample_file) from split_POP_samples_2
        each sample_size from params.reduce_sample_size
    output:
        set val(POP), file(sample_out), val(sample_size) into reduce_POP_samples_all
    script:
        sample_out = "${POP}_reduce${sample_size}.sample"
        """
        python2.7 ${params.homedir}/templates/random_sample.py \
            --sampleIN ${sample_file} \
            --sampleOUT ${sample_out} \
            --sample_size ${sample_size} \
            --dataset ${POP}
        """
}

'''
Step 2.13: Split vcf per population
'''
biallOnly_baylor.into { biallOnly_baylor; biallOnly_baylor_split_pop}
split_POP_samples_all.into {split_POP_samples_all; split_POP_samples_all__split_pop}
pop_sample_data = split_POP_samples_all__split_pop.combine(biallOnly_baylor_split_pop)
process split_vcf_per_pop {
    echo true
    tag { "split_vcf_${POP}_${chrm}_${file(vcf_file.baseName).baseName}" }
    memory { 5.GB * task.attempt }
    time{ 2.hour * task.attempt }
    publishDir "${params.work_dir}/VCF_POP/BAYLOR/VCF/${POP}", overwrite: true, mode:'symlink'
    input:
        set val(POP), file(POP_sample_file), val(chrm), file(vcf_file), file(vcf_file_tbi) from pop_sample_data
    output:
        set val(POP), val(chrm), file(vcf_out) into split_vcf_per_pop
    script:
        vcf_out = "${POP}_${vcf_file}"
        """
        vcftools \
            --gzvcf ${vcf_file} \
            --keep ${POP_sample_file} \
            --recode --recode-INFO-all -c | \
        bgzip -c > ${vcf_out}
        bcftools index --tbi -f ${vcf_out}
        """
}


////'''
////Step 2.13.1: Split vcf per reduced population
////'''
////biallOnly_baylor.into { biallOnly_baylor; biallOnly_baylor_split_pop_1}
////reduce_POP_samples_all.into {reduce_POP_samples_all; reduce_POP_samples_all__split_pop}
////pop_sample_data_1 = reduce_POP_samples_all__split_pop.combine(biallOnly_baylor_split_pop_1)
////process split_vcf_per_pop_reduced {
////    echo true
////    tag { "split_vcf_${POP}_${chrm}_${file(vcf_file.baseName).baseName}" }
////    memory { 5.GB * task.attempt }
////    time{ 2.hour * task.attempt }
////    publishDir "${params.work_dir}/VCF_POP/BAYLOR/VCF/${POP}", overwrite: true, mode:'symlink'
////    input:
////        set val(POP), file(POP_sample_file), val(chrm), file(vcf_file), file(vcf_file_tbi) from pop_sample_data_1
////    output:
////        set val(POP), val(chrm), file(vcf_out) into split_vcf_per_pop_reduced
////    script:
////        vcf_out = "${POP}_${vcf_file}_reduced"
////        """
////        vcftools \
////            --gzvcf ${vcf_file} \
////            --keep ${POP_sample_file} \
////            --recode --recode-INFO-all -c | \
////        bgzip -c > ${vcf_out}
////        bcftools index --tbi -f ${vcf_out}
////        """
////}
//
'''
Step 2.14: Concatenate chromosome VCFs into one
'''
split_vcf_per_pop.into { split_vcf_per_pop; split_vcf_per_pop__merge_vcf_pop}
merge_vcf_pop_list = [:]
split_vcf_per_pop__merge_vcf_pop_list = split_vcf_per_pop__merge_vcf_pop.toSortedList().val
split_vcf_per_pop__merge_vcf_pop_list.each { POP, chrm, vcf_chrm_pop ->
    if ( !(POP in merge_vcf_pop_list.keySet()) ) {
        merge_vcf_pop_list[POP] = [POP]
    }
    merge_vcf_pop_list[POP] << vcf_chrm_pop
}
merge_vcf_pop_cha = Channel.from(merge_vcf_pop_list.values())
merge_vcf_pop_cha.into { merge_vcf_pop_cha; merge_vcf_pop__merge}
process merge_vcf_pop {
    tag { "merge_vcf_pop_${POP}" }
    memory { 4.GB * task.attempt }
    time{ 2.hour * task.attempt }
    publishDir "${params.work_dir}/VCF_POP/BAYLOR/VCF/${POP}", overwrite: true, mode:'symlink'
    input:
        val datas from merge_vcf_pop__merge
    output:
        set val(POP), file(vcf_out) into merge_vcf_pop
    script:
        POP = datas[0]
        vcf_data = datas[1..-1].join(' ')
        vcf_f = file(file(datas[1]).baseName).baseName.split('_')
        vcf_out = "${POP}_${vcf_f[1].split('\\.')[0..1].join('.')}_${vcf_f[2..-1].join('_')}.vcf.gz"
        """
        bcftools concat ${vcf_data} -Oz -o ${vcf_out}
        bcftools index --tbi -f ${vcf_out}
        """
}


//'''
//Step 2.14.1: Concatenate chromosome VCFs into one (Reduced)
//'''
//split_vcf_per_pop_reduced.into { split_vcf_per_pop_reduced; split_vcf_per_pop_reduced__merge_vcf_pop}
//merge_vcf_pop_list_reduced = [:]
//split_vcf_per_pop_reduced__merge_vcf_pop_list = split_vcf_per_pop__merge_vcf_pop.toSortedList().val
//split_vcf_per_pop_reduced__merge_vcf_pop_list.each { POP, chrm, vcf_chrm_pop ->
//    if ( !(POP in merge_vcf_pop_list_reduced.keySet()) ) {
//        merge_vcf_pop_list_reduced[POP] = [POP]
//    }
//    merge_vcf_pop_list_reduced[POP] << vcf_chrm_pop
//}
//merge_vcf_pop_reduced_cha = Channel.from(merge_vcf_pop_list_reduced.values())
//merge_vcf_pop_reduced_cha.into { merge_vcf_pop_reduced_cha; merge_vcf_pop_reduced__merge}
//process merge_vcf_pop_reduced {
//    tag { "merge_vcf_pop_${POP}" }
//    memory { 4.GB * task.attempt }
//    time{ 2.hour * task.attempt }
//    publishDir "${params.work_dir}/VCF_POP/BAYLOR/VCF/${POP}", overwrite: true, mode:'symlink'
//    input:
//        val datas from merge_vcf_pop_reduced__merge
//    output:
//        set val(POP), file(vcf_out) into merge_vcf_pop_reduced
//    script:
//        POP = datas[0]
//        vcf_data = datas[1..-1].join(' ')
//        vcf_f = file(file(datas[1]).baseName).baseName.split('_')
//        vcf_out = "${POP}_${vcf_f[1].split('\\.')[0..1].join('.')}_${vcf_f[2..-1].join('_')}.vcf.gz"
//        """
//        bcftools concat ${vcf_data} -Oz -o ${vcf_out}
//        bcftools index --tbi -f ${vcf_out}
//        """
//}
//
'''
Step 2.15.1: DAF by chromosomes
'''
split_vcf_per_pop.into { split_vcf_per_pop; split_vcf_per_pop_2}
process daf_by_pop_chrm {
    tag { "daf_${chrm}_${POP}" }
    memory { 5.GB * task.attempt }
    time{ 2.hour * task.attempt }
    publishDir "${params.work_dir}/VCF_POP/BAYLOR/DAF/${POP}", overwrite: true, mode:'symlink'
    input:
        set val(POP), val(chrm), file(vcf_pop) from split_vcf_per_pop_2
    output:
        set val(POP), val(chrm), file("${vcf_out}.frq.count") into daf_by_pop_chrm
    script:
        vcf_out = "${file(vcf_pop.baseName).baseName}.daf"
        """
        vcftools \
            --gzvcf ${vcf_pop} \
            --counts --derived \
            --out ${vcf_out}
        vcftools \
            --gzvcf ${vcf_pop} \
            --counts --derived \
            --out ${vcf_out}
        """
}

//
//'''
//Step 2.15
//'''
//merge_vcf_pop.into { merge_vcf_pop; merge_vcf_pop__daf}
//process daf_by_pop {
//    tag { "daf_${POP}" }
//    memory { 5.GB * task.attempt }
//    time{ 2.hour * task.attempt }
//    publishDir "${params.work_dir}/VCF_POP/BAYLOR/DAF/${POP}", overwrite: true, mode:'symlink'
//    input:
//        set val(POP), file(vcf_pop) from merge_vcf_pop__daf
//    output:
//        set val(POP), file("${vcf_out}.frq.count") into daf_by_pop
//    script:
//        vcf_out = "${file(vcf_pop.baseName).baseName}.daf"
//        """
//        vcftools \
//            --gzvcf ${vcf_pop} \
//            --counts --derived \
//            --out ${vcf_out}
//        """
//}

//'''
//Step 2.16:
//'''
//process sfs_by_pop {
//    tag { "sfs_${POP}" }
//    memory { 10.GB * task.attempt }
//    time{ 2.hour * task.attempt }
//    publishDir "${params.work_dir}/VCF_POP/BAYLOR/SFS/${POP}", overwrite: true, mode:'symlink'
////    publishDir "${params.group_dir}/VCF_POP/SFS", overwrite: true, mode:'symlink'
//    input:
//        set val(POP), file(vcf_frq) from daf_by_pop
//    output:
//        set val(POP), file(sfs_plot) into sfs_by_pop
//    script:
//        popID = POP
//        frq_input = vcf_frq
//        sfs_plot = "${vcf_frq.baseName}.sfs.pdf"
//        template "CreateSFSandDerivedCounts.R"
//}


//
//'''
//Extract number of singletons per population
//'''
//merge_vcf_pop.into { merge_vcf_pop; merge_vcf_pop__singl}
//process singl_by_pop {
//    echo true
//    tag { "singl_${POP}" }
//    memory { 4.GB * task.attempt }
//    time{ 2.hour * task.attempt }
//    publishDir "${params.work_dir}/VCF_POP/SINGL", overwrite: true, mode:'symlink'
////    publishDir "${params.group_dir}/VCF_POP/SINGL", overwrite: true, mode:'symlink'
//    input:
//        set val(POP), file(vcf_pop) from merge_vcf_pop__singl
//    output:
//        set val(POP), file("${POP}.singletons.per.sample") into singl_by_pop
//    script:
//        vcf_out = "${POP}_dp6_anc_f_dbsnp_snpeff_derived"
//        """
//        vcftools \
//            --gzvcf ${vcf_pop} \
//            --singletons --derived \
//            --out ${vcf_out}
//        grep ${POP} ${params.sample_file} | cut -f1 > ${POP}.sample
//        echo "singletons" > ${POP}.singleton.count
//        while read SAMPLE; do
//            grep \$SAMPLE ${vcf_out}.singletons | wc -l;
//        done < ${POP}.sample >> ${POP}.singleton.count
//        echo "sample" > ${POP}_.sample | cat ${POP}.sample >> ${POP}_.sample
//        paste ${POP}_.sample ${POP}.singleton.count > ${POP}.singletons.per.sample
//        """
//}

//
//'''
//Step 12.6: Plot average singletons per population
//'''
//popFrq_4_all.into{ popFrq_4_all; popFrq_4_1 }
//process plot_pgxAverageSingletonsPerPopulation{
//    tag "plot_lof_${dataset}"
//    cpus { 2 * task.attempt }
//    memory { 2.GB * task.cpus }
//    publishDir "${params.work_dir}/data/PGX_ONLY/FRQ/${dataset}", overwrite: true, mode:'copy'
//    input:
//        set val(dataset), file(dataset_singletons_per_sample) from popFrq_4_1
//    output:
//        set val(dataset), file(dataset_pgxAverageSingletonsPerPopulation) into plot_pgxAverageSingletonsPerPopulation
//    script:
//        sample_file = params.sample_file
//        dataset_pgxAverageSingletonsPerPopulation = "${dataset}_pgxAverageSingletonsPerPopulation.tiff"
//        template "step12_6_average_singletons_per_population.R"
//}

//
//'''
//Step 1.8: Filter out related individuals and remove sites with high missingness
//'''
//biallOnly_baylor.into { biallOnly_baylor; biallOnly_baylor_1}
//process vcf_sites_only_baylor {
//    tag "Sites_only_${file(vcf_file.baseName).baseName}"
//    memory { 10.GB * task.attempt }
//    time { 5.hour * task.attempt }
//    publishDir "${params.work_dir}/VCF_FILTERED/SITE_ONLY/BAYLOR", overwrite: true, mode:'copy'
//    input:
//        set val(chrm), file(vcf_file), file(vcf_file_tbi) from biallOnly_baylor_1
//    output:
//        set val(chrm), file(vcf_out), file("${vcf_out}.tbi") into vcf_sites_only_baylor
//    script:
//        vcf_out = "${file(vcf_file.baseName).baseName}_sitesOnly.vcf.gz"
//        """
//        bcftools \
//            view ${vcf_file} \
//            --drop-genotypes \
//            -Oz -o ${vcf_out}
//        bcftools index --tbi -f ${vcf_out}
//        """
//}

//
//'''
//Step 2.7: Annotate whole baylor database with snpEff using dbNSFP
//'''
//annotate_cosmic_baylor.into { annotate_cosmic_baylor; annotate_cosmic_baylor_1}
//process annotate_dbnsfp_baylor {
//    tag "cosmic_${chrm}_${file(vcf_file.baseName).baseName}"
//    memory { 5.GB * task.attempt }
//    time{ 6.hour * task.attempt }
//    publishDir "${params.work_dir}/VCF_ANN/BAYLOR", overwrite: true, mode:'symlink'
//    input:
//        set val(chrm), file(vcf_file) from annotate_cosmic_baylor_1
//    output:
//        set val(chrm), file("${vcf_out}.gz") into annotate_dbnsfp_baylor
//    script:
//        vcf_out = "${file(vcf_file.baseName).baseName}_dbnsfp.vcf"
//        """
//        SnpSift \
//            annotate \
//            ${params.dbnsfp} \
//            ${vcf_file} \
//            -v -c ${params.snpeff_config} \
//            > ${vcf_out}
//        bgzip -f ${vcf_out}
//        bcftools index --tbi -f ${vcf_out}.gz
//        """
//}


//'''
//Step : NOT USED .... Download Genome VCF Files from http://gnomad.broadinstitute.org/downloads
//'''
//vcf_data_merge.into{ vcf_data_merge; vcf_data_merge_1}
//process download_gnomad_genome_vcfs {
//    tag "gnomad_g${chrm}"
//    maxForks 15
//    memory { 5.GB * task.attempt }
//    time{ 6.hour * task.attempt }
//    publishDir "${params.ref_dir}/gnomAD", overwrite: true, mode:'symlink'
//    input:
//        val chrm from CHRMS
//    output:
//        set val(chrm), file(gFile) into download_gnomad_genomes
//    script:
//        inFile = sprintf(params.gnomad_genome, chrm)
//        gFile = sprintf("gnomad.genomes.r2.0.1.sites.%s.vcf.gz", chrm)
//        """
//        wget ${inFile} -O ${gFile}
//        """
//}

//
//'''
//Step 2: Annotate whole baylor database with snpEff using dbSNP database
//'''
//main.into { main; annotate_vcf_1}
//process annotate_dbsnp_snpeff_all {
//    echo true
//    tag "dbSNP_${chrm}"
//    memory { 4.GB * task.attempt }
//    time{ 6.hour * task.attempt }
//    publishDir "${params.work_dir}/VCF_ANN/", overwrite: true, mode:'symlink'
//    input:
//        set val(chrm), file(vcf_file) from annotate_vcf_1
//    output:
//        set val(chrm), file("${vcf_out}.gz") into annotate_dbsnp_snpeff_all
//    script:
//        vcf_out = "${params.prefix_new}${chrm}_dp6_anc_f_dbsnp.vcf"
//        """
//        SnpSift \
//            annotate \
//            ${params.dbsnp_vcf} \
//            -c ${params.snpeff_config} \
//            ${vcf_dp6_anc_f} > ${vcf_out} -v
//        bgzip -f ${vcf_out}
//        """
//}

//
//'''
//Step 1: Filter VCF by coverage depth using bcftools
//'''
//vcf_data.into{ vcf_data; vcf_data_1}
//process filter_vcf_by_depth {
//    tag "dp6_${chrm}"
//    memory { 4.GB * task.attempt }
//    time{ 2.hour * task.attempt }
//    publishDir "${params.work_dir}/VCF_FILTERED/", overwrite: true, mode:'symlink'
////    publishDir "${params.group_dir}/VCF_FILTERED/", overwrite: true, mode:'symlink'
//    input:
//        set val(chrm), file(vcf_file) from vcf_data_1
//    output:
//        set val(chrm), file("${params.prefix_new}${chrm}_dp6.vcf.gz"), file("${params.prefix_new}${chrm}_dp6.vcf.gz.tbi") into vcf_6depth
//    script:
//        """
//        bcftools view -i 'DP>6' ${vcf_file} | \
//            bgzip -c > ${params.prefix_new}${chrm}_dp6.vcf.gz
//        bcftools index --tbi -f ${params.prefix_new}${chrm}_dp6.vcf.gz
//        """
//}
//
//
//'''
//Step 2: Annotate VCF for Ancestral Allele (AA) using in-house python script
//'''
//vcf_6depth.into { vcf_6depth; vcf_6depth__anc } // Duplicate channel so that it can be used multiple times
//process add_ANC_to_VCF {
//    echo true
//    tag "AA_${chrm}"
//    memory { 4.GB * task.attempt }
//    time{ 4.hour * task.attempt }
//    publishDir "${params.work_dir}/VCF_ANC/", overwrite: true, mode:'symlink'
////    publishDir "${params.group_dir}/VCF_ANC/", overwrite: true, mode:'symlink'
//    input:
//        set val(chrm), file(vcf_file), file(vcf_file_tbi) from vcf_6depth__anc
//    output:
//        set val(chrm), file("${vcf_file_out}.gz"), file("${vcf_file_out}.gz.tbi") into vcf_anc
//    script:
//        vcf_file_out = "${params.prefix_new}${chrm}_dp6_anc.vcf"
//        """
//        gunzip -c ${vcf_file} > ${vcf_file.baseName}
//        ${params.python_27_env} ${params.homedir}/templates/add-ANC-to-vcf_new.py -g --in ${vcf_file.baseName} --out ${vcf_file_out} --genomedata ${params.genomedata_path}
//        bgzip -f ${vcf_file_out}
//        bcftools index --tbi -f ${vcf_file_out}.gz
//        rm -f ${vcf_file.baseName}
//        """
//}
//
//
//'''
//Step 3: Filter sites with AA
//'''
//vcf_anc.into { vcf_anc; vcf_anc__only}
//process add_ANC_to_VCF_only {
//    echo true
//    tag "AA_only_${chrm}"
//    memory { 1.GB * task.attempt }
//    time{ 2.hour * task.attempt }
//    publishDir "${params.work_dir}/VCF_ANC/", overwrite: true, mode:'symlink'
////    publishDir "${params.group_dir}/VCF_ANC/", overwrite: true, mode:'symlink'
//    input:
//        set val(chrm), file(vcf_dp6_anc), file(vcf_dp6_anc_tbi) from vcf_anc__only
//    output:
//        set val(chrm), file("${params.prefix_new}${chrm}_dp6_anc_f.vcf.gz"), file("${params.prefix_new}${chrm}_dp6_anc_f.vcf.gz.tbi") into vcf_anc_f
//    script:
//        vcf_dp6_anc_out = "${params.prefix_new}${chrm}_dp6_anc_f.vcf.gz"
//        """
//        bcftools view -i 'AA!="." & AA!="-" & AA!="N"' ${vcf_dp6_anc} | bgzip -c > ${vcf_dp6_anc_out}
//        bcftools index --tbi -f ${vcf_dp6_anc_out}
//        """
//}
//
//// TODO Download dbSNP database if not exists
//vcf_anc_f.into { vcf_anc_f; vcf_anc_f__dbsnp}
//process download_snpeff_db {
//    echo true
//    tag "down_db"
//    memory { 8.GB * task.attempt }
//    time{ 6.hour * task.attempt }
//    publishDir "${params.work_dir}/", overwrite: true, mode:'symlink'
//    output:
//        set val(chrm), file("${vcf_out}.gz") into annotate_dbsnp_snpeff_all
//    script:
//        """
//        wget http://downloads.sourceforge.net/project/snpeff/databases/v4_3/snpEff_v4_3_GRCh37.75.zip
//        touch
//        """
//}


//'''
//Step 4:
//'''
//vcf_anc_f.into { vcf_anc_f; vcf_anc_f__dbsnp}
//process annotate_dbsnp_snpeff {
//    echo true
//    tag "dbSNP_${chrm}"
//    memory { 4.GB * task.attempt }
//    time{ 6.hour * task.attempt }
//    publishDir "${params.work_dir}/VCF_ANN/", overwrite: true, mode:'symlink'
////    publishDir "${params.group_dir}/VCF_ANN/", overwrite: true, mode:'symlink'
//    input:
//        set val(chrm), file(vcf_dp6_anc_f) from vcf_anc_f__dbsnp
//    output:
//        set val(chrm), file("${vcf_out}.gz") into annotate_dbsnp_snpeff_all
//    script:
//        vcf_out = "${params.prefix_new}${chrm}_dp6_anc_f_dbsnp.vcf"
//        """
//        SnpSift \
//            annotate \
//            ${params.dbsnp_vcf} \
//            -c ${params.snpeff_config} \
//            ${vcf_dp6_anc_f} > ${vcf_out} -v
//        bgzip -f ${vcf_out}
//        """
//}
//
//
//'''
//Step 5
//'''
//process annotate_snpeff {
//    echo true
//    tag { "snpEff_${chrm}" }
//    memory { 4.GB * task.attempt }
//    time{ 6.hour * task.attempt }
//    publishDir "${params.work_dir}/VCF_ANN/", overwrite: true, mode:'symlink'
////    publishDir "${params.group_dir}/VCF_ANN/", overwrite: true, mode:'symlink'
//    input:
//        set val(chrm), file(vcf_dp6_anc_f_dbsnp) from annotate_dbsnp_snpeff_all
//    output:
//        set val(chrm), file("${vcf_out}.gz") into annotate_snpeff_all
//    script:
//        vcf_in = "${params.prefix_new}${chrm}_dp6_anc_f_dbsnp_snpeff"
//        vcf_out = "${params.prefix_new}${chrm}_dp6_anc_f_dbsnp_snpeff.vcf"
//        """
//        ${params.snpEff} \
//            ${params.snpEff_human_db} \
//            -stats ${vcf_in}.html \
//            -csvStats ${vcf_in}.csv \
//            -dataDir ${params.snpEff_database} \
//            -c ${params.snpeff_config} \
//            ${vcf_dp6_anc_f_dbsnp} > ${vcf_out} -v -nodownload
//        bgzip -f ${vcf_out}
//        """
//}
//
//

// TODO Not working
//'''
//Step: Copy files to group shared folder
//'''
//daf_by_chrm_pop.into { daf_by_chrm_pop; daf_by_chrm_pop_5 }
//copyToshared_ = daf_by_chrm_pop_5.collect()
//process copyToshared {
//    tag { "copyToshared" }
//    memory { 2.GB * task.attempt }
//    time{ 1.hour * task.attempt }
//    validExitStatus 0,1,2
//    beforeScript "source activate ngs"
//    input:
//        val(dat) from copyToshared_
//    output:
//        file(out) into copyToshared
//    script:
//        source = "${params.work_dir}/samples ${params.work_dir}/VCF_ANN ${params.work_dir}/VCF_FILTERED ${params.work_dir}/VCF_POP "
//        dest = "${params.group_dir}"
//        Date date = new Date()
//        String datePart = date.format("dd_MM_yyyy")
//        String timePart = date.format("HH_mm_ss")
//        out = "copy_${datePart}_${timePart}.txt"
//        """
//        mkdir -p ${params.group_dir}
//        rsync -avzr --progress ${source} ${dest}
//        echo \' rsync -avzr --progress ${source} > ${dest} \' > ${out}
//        """
//}
//copyToshared.subscribe{ println "|-- Results copied to ${it}" }

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