nextflow.enable.dsl=2

params.input = ""
params.rfile = ""
process smoove {
  publishDir "PipeRes", mode: 'copy'
	input:
	      tuple(val(sample_name), path(bam), path(index))
        path(fasta)
        path(fai)
	output:
	      tuple(file("${sample_name}_smoove.vcf.gz"),file("${sample_name}_smoove.vcf.gz.tbi") , val(sample_name))	
	script:
	"""
  echo ${sample_name} > sample.txt
	smoove call  --outdir . \
	 	--name ${sample_name} --fasta ${fasta}  -p ${task.cpus} ${bam} 
	bcftools view -O u -o ${sample_name}.R.bcf ${sample_name}-smoove.vcf.gz
  bcftools sort --temp-dir /scratch/dysgu  -m 2G -O z -o ${sample_name}-smoove.vcf.gz  ${sample_name}.R.bcf
  bcftools reheader -s sample.txt ${sample_name}-smoove.vcf.gz > ${sample_name}_smoove.vcf.gz
  bcftools index --tbi ${sample_name}_smoove.vcf.gz

  """

	}


process manta {
    errorStrategy 'terminate' // TODO: change after debugging is done
    // container = 'docker://brentp/rare-disease-sv:v0.1.2'
    publishDir "PipeRes", mode: 'copy'
    shell = ['/bin/bash', '-euo', 'pipefail']
    input:
         tuple(val(sample_name), path(bam), path(index))
         path(fasta)
         path(fai)
         path(rfile)
         path(rindex)
    output:
         tuple(file("${sample_name}_diploidSV.vcf.gz"),file("${sample_name}_diploidSV.vcf.gz.tbi"), val(sample_name))
   script:
         """
         echo ${sample_name} > sample.txt
         configManta.py --bam ${bam} --referenceFasta $fasta --runDir . --callRegions $rfile
         python2 ./runWorkflow.py -j ${task.cpus} -m local 
         bcftools reheader -s sample.txt results/variants/diploidSV.vcf.gz >${sample_name}_diploidSV.vcf.gz
         bcftools index -f --tbi ${sample_name}_diploidSV.vcf.gz
         """
}
process expansion_hunter {
    errorStrategy 'terminate' // TODO: change after debugging is done

    publishDir "PipeRes", mode: 'copy'
    shell = ['/bin/bash', '-euo', 'pipefail']

    input:
        tuple(val(sample_name), path(bam), path(index))
        path(fasta)
        path(fai)
        path(catalog)

    output:
        tuple(file("${sample_name}_str.vcf"), val(sample_name))

script:
    output_file = "${sample_name}_str.vcf"
    Expansion_Dir="/home/moabd/projects/def-wyeth/moabd/ExpansionHunter-v5.0.0-linux_x86_64/bin"
    """
    ${Expansion_Dir}/ExpansionHunter \
    --output-prefix ${sample_name} --reference $fasta  \
    --reads ${bam} \
    --variant-catalog ${catalog} -n ${task.cpus} 
    mv ${sample_name}.vcf ${output_file}
    """

}

process melt {
// merge all files into one file
// add gene annotation file as a parameter
 publishDir "PipeRes", mode: 'copy'
      shell = ['/bin/bash', '-euo', 'pipefail']
      module = "java/1.8.0_192:bowtie2:bcftools:tabix"   
      input:
            tuple(val(sample_name), path(bam), path(index))
            path(fasta)
            path(fai)
            path(tfile)
            path(annotation_file)
      output:
            tuple(file("${sample_name}_mei.vcf.gz"),file("${sample_name}_mei.vcf.gz.tbi"), val(sample_name))
      script:
        MELT_DIR="/home/moabd/projects/def-wyeth/moabd/MELTv2.2.2"
        output_file="${sample_name}_mei.vcf.gz"
         """
            mkdir -p  ${sample_name}
            java -Xmx8G -jar ${MELT_DIR}/MELT.jar Single \
            -b hs37d5 \
            -t ${tfile}  \
            -h $fasta \
            -bamfile $bam \
            -w ${sample_name} \
            -n ${annotation_file}
            # fix issues with  MELT vcfs
            for name in {ALU,SVA,LINE1};do bcftools annotate -x FMT/GL ${sample_name}/\${name}.final_comp.vcf > \${name}.vcf;done
            for name in {ALU,SVA,LINE1};do bcftools view -O u -o \${name}.bcf \${name}.vcf;done
            for name in {ALU,SVA,LINE1};do bcftools sort  -m 2G -O z -o \${name}.vcf.gz \${name}.bcf;done
            for name in {ALU,SVA,LINE1};do bcftools index --tbi \${name}.vcf.gz;done
            bcftools concat -a -Oz  -o ${output_file} *vcf.gz
            bcftools index --tbi ${output_file}

       """

}

process concat_by_sample {
    errorStrategy 'terminate' // TODO: change after debugging is done

    publishDir "SVRes", mode: 'copy'
    shell = ['/bin/bash', '-euo', 'pipefail']

    input:
        tuple(path(vcfs), path(indexes), val(sample_name))
    output: file("${output_file}")

    script:
    output_file = "${sample_name}.concat-svs.vcf"
    """
    bcftools concat -a -O v -o ${output_file} *.vcf.gz
    """
}
process jasmine {
    container = 'docker://brentp/rare-disease-sv:v0.1.2'
    publishDir "SVRes/", mode: 'copy'
    shell = ['/bin/bash', '-euo', 'pipefail']

    input:
        val(sample_vcfs)
  	path(fasta)
        path(fai)
    output: tuple(file("${output_file}"), file("${output_file}.tbi"))

    script:
    output_file = "jasmine.merged.vcf.gz"
    file("$workDir/vcfs.list").withWriter { fh ->
            sample_vcfs.each { vcf ->
                // write .vcf even though it's really .vcf.gz since jasmine doesn't accept gz
                // and we change the file below.
                fh.write(vcf.toString()); fh.write("\n")
            }
    }
    // we don't merge if we only have a single sample.
    if(sample_vcfs.size() > 1) {
        """
        # jasmine can't do gzip.
        jasmine  -Xmx6g --dup_to_ins  --threads ${task.cpus} --output_genotypes --allow_intrasample --file_list=${workDir}/vcfs.list out_file=${output_file}  genome_file=$fasta
        # NOTE: removing BNDs and setting min start to > 150 as paragraph fails if start < readlength
        tiwih setsvalt --drop-bnds --inv-2-ins -o ${output_file}.tmp.vcf.gz $fasta $output_file
        bcftools sort --temp-dir /scratch/tmp/jasmine  -m 2G -O z -o ${output_file} ${output_file}.tmp.vcf.gz
        bcftools index --tbi $output_file
        """
    } else {
        """
        tiwih setsvalt --drop-bnds -o ${output_file} $fasta ${sample_vcfs[0]}
        tabix $output_file
        """
    }

}

process paragraph_duphold {
  errorStrategy 'terminate' // TODO: change after debugging is done
  shell = ['/bin/bash', '-euo', 'pipefail']
  container = 'docker://brentp/rare-disease-sv:v0.1.2'

  publishDir "SVRes/", mode: 'copy'

  input:
     tuple(path(site_vcf), path(site_vcf_index))
     tuple(val(sample), path(bam), path(index))
     path(fasta)
     path(fai)

  output: tuple(path("${output_file}"), path("${output_file}.csi"))

  script:
  output_file = "${sample}.paragraph.vcf.gz"

  """
    dp=\$(tiwih meandepth $bam)
    tsample=\$(tiwih samplename $bam)
    echo "id\tpath\tdepth\tread length" > sample.manifest
    echo "\$tsample\t$bam\t\$dp\t150" >> sample.manifest
    M=\$((dp * 5))
    cat sample.manifest
    
    # this is the main paragraph entrypoint
    multigrmpy.py -i $site_vcf \
        -m sample.manifest \
        -r $fasta \
        -o t \
        -t ${task.cpus} \
        -M \$M


    # duphold adds depth annotations looking at coverage fold-change around Svs
    duphold -d -v t/genotypes.vcf.gz -b $bam -f $fasta -t 4 -o $output_file
    bcftools index --threads 3 $output_file
  """

}

params.help = false
if (params.help) {
    log.info("""
-------------------------------------------------------------------------
rare-disease-wf sv
Required Arguments:
   --input           A csv file with columns for sample, bam, and bam index
   --fasta           Path to reference fasta
   --tfile           Path to transposon file for MELT
   --rfile            Path to bgziped region bed file
Optional Arguments:
-------------------
   --cohort_name     optional name for the cohort (default: "rare-disease")
   --output_dir      optional name for where to place results (default: "results-rare-disease")
   """)

}

// process bcftools {
//     publishDir "NewRes", mode: 'copy'
//     input:
//         file vcf from vcf_ch
//     output:
//         file "*_stats.txt" into stats_ch
//     script:
//     """
//     bcftools stats -s - $vcf > ${vcf}_stats.txt
//     """
//     }

// process multiqc {
//     publishDir "NewRes", mode: 'copy'
// 
//     input:
//     file('Res/*')   from stats_ch.collect()
// 
//     output:
//     file "multiqc_report.html" 
//     file "multiqc_data"
// 
//     script:
//     """
//     multiqc .
//     """
// 
//     }

workflow {
    //  Input Data //

    input = Channel.fromPath(params.input, checkIfExists: true)
            .splitCsv(header:true)
            .map{ row-> tuple(row.sample, file(row.bam), file(row.index))  }
    fasta = file(params.fasta, checkIfExists: true)
    fai = file(params.fasta + ".fai", checkIfExists: true)
    tfile = file(params.tfile)
    rfile = file(params.rfile)
    rindex = file(params.rfile + ".tbi")
    annotation_file = file(params.annotation_file)
    catalog = file(params.catalog)
    
    // Actual Workflow //

    sm = smoove(input, fasta, fai).view()
    mr = manta(input, fasta, fai,rfile,rindex).view()
    melt(input, fasta, fai,tfile, annotation_file)
    expansion_hunter(input, fasta, fai, catalog)
    sv_groups = mr.concat(sm) | groupTuple(by: 2)
    sv_groups.view()
    svs = concat_by_sample(sv_groups) | collect
    svs.view()
    sv_merged = jasmine(svs, fasta, fasta + ".fai")
    genotyped = paragraph_duphold(sv_merged, input, fasta, fasta + ".fai")
}
