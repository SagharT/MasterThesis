'''
This snake fike automatically picks up new .raw files from the specified directory for processing.
It includes Parsing, Feature extraction, Quality control, peptide identification and creating reports
Check the config file and edit the nessesary paths
'''
configfile: "config.yaml"

import os

# Generate the list of samples to analyse
samples_to_analyse = [sample.replace('.raw', '') for sample in os.listdir(config["path_to_samples"]) if sample.endswith(".raw")]

rule all:
    input:
        expand(config["path_to_csv_output"] + "{sample}.features.csv", sample=samples_to_analyse),
        expand(config["qc_summary_folder_path"]+"{sample}.mzmlsummary.txt", sample=samples_to_analyse),
        expand(config["qc_summary_folder_path"]+"{sample}.featuresummary.txt", sample=samples_to_analyse), #feature detection (Dinosaur)
        expand(config["path_to_report"]+"{sample}.pdf", sample=samples_to_analyse),
        expand(config["path_to_report"]+"{sample}.csv", sample=samples_to_analyse),
        config["path_to_report"] + "Comparison.pdf",
        config["path_to_report"]+ ".cleanup_complete"

        
# Raw file conversion 
rule convert_raw_file:
    input:
        raw = config["path_to_samples"] + "{sample}.raw"
    params:
        rawfile = "{sample}.raw"
    output:
        mzML = config["path_to_mzml"] + "{sample}.mzML"
        
    shell:
        """
        set -e
        if docker run --rm -v ./{config[path_to_samples]}:/indata \
            -v ./{config[path_to_mzml]}:/outdata chambm/pwiz-skyline-i-agree-to-the-vendor-licenses wine msconvert \
            -o /outdata -v --filter="peakPicking true 1-" --filter="demultiplex" /indata/{params.rawfile}; then
            echo "Success with demultiplex filter for {wildcards.sample}"
        else
            echo "Failed with demultiplex, trying peakPicking true 1- for {wildcards.sample}"
            docker run --rm -v ./{config[path_to_samples]}:/indata \
            -v ./{config[path_to_mzml]}:/outdata chambm/pwiz-skyline-i-agree-to-the-vendor-licenses wine msconvert \
            -o /outdata -v --filter="peakPicking true 1-" /indata/{params.rawfile} && \
            echo "Success with peakPicking true 1- for {wildcards.sample}" 
        fi
        sudo chmod -R ou+w {config[path_to_mzml]}
        """

# mzML file parsing
rule parse_mzML:
    input:
        mzML = config["path_to_mzml"] + "{sample}.mzML"
    output:
        csv = config["path_to_csv_output"] + "{sample}.features.csv"
    shell:
        "python3 mzML-Parser.py {input.mzML} {output.csv}"

# Feature extraction, Dinosaur for both DIA and DDA
rule feat_extraction:
    input:
        config["path_to_mzml"]+"{sample}.mzML"
    params:
        qc = config["dinosaur_output_path"]+"{sample}.qc.zip"
    output:
        config["dinosaur_output_path"]+"{sample}.features.tsv"
    shell:
        """
        java {config[max_heap_size]} -Djava.awt.headless=true \
          -jar {config[dinosaur_jar_path]} \
          --zipQcFolder=true \
          {config[dinosaur_params]} \
          --outDir={config[dinosaur_output_path]} \
          {input}
        """
# Quality control summaries for both DIA and DDA (mzmlsummary) and (featuresummary)
rule qc_mzML_summary:
    input:
        config["path_to_mzml"]+"{sample}.mzML"
    output:
        config["qc_summary_folder_path"]+"{sample}.mzmlsummary.txt"
    shell:
        """
        java -cp {config[qc_lcms_jar_path]} \
          se.nbis.omics.MzMLSummary \
          {input} \
          >{output}
        """

rule qc_feature_summary:
    input:
        config["dinosaur_output_path"]+"{sample}.features.tsv"
    output:
        config["qc_summary_folder_path"]+"{sample}.featuresummary.txt"
    shell:
        """
        java -cp {config[qc_lcms_jar_path]} \
          se.nbis.omics.FeaturesQCSummary \
          {input} \
          >{output}
        """

############---------- DDA ----------############
# Peptide identification msgfplus
rule peptide_identification:
    input:
        mzml = config["path_to_mzml"]+"{sample}.mzML",
        db = config["fasta"],
        featurescsv = config["path_to_csv_output"] + "{sample}.features.csv",
    params:
        sample = "{sample}",
    output:
       mzid = config["peptide_id_output"]+"{sample}.mzid"
    shell:
        """
        set -e
        if python3 DIA-or-DDA.py {input.featurescsv}; then       
            echo dia;
            touch {output};
        else
        java {config[max_heap_size]} -Djava.awt.headless=true \
         -jar {config[msgfplus_jar_path]} \
        {config[msgfplus_extra]} \
         -mod {config[msgfplus_mod_path]} \
         -d {input.db} \
         -o {output.mzid} \
         -s {input.mzml}
        fi
       """
rule qc_peptideid_mzid_summary:
    input:
        mzid = config["peptide_id_output"]+"{sample}.mzid",
        featurescsv = config["path_to_csv_output"] + "{sample}.features.csv",
    output:
        config["qc_summary_folder_path"]+"{sample}.mzidsummary.txt"
    shell:
        """
        set -e
         if python3 DIA-or-DDA.py {input.featurescsv}; then
            echo dia;
            touch {output};
        else
          java -cp {config[qc_lcms_jar_path]} \
          se.nbis.omics.MzIdentMLSummary \
          {input.mzid} \
          {config[qc_fdr]} \
          >{output}
         java -cp {config[qc_lcms_jar_path]} se.nbis.omics.QCMerge \
          {config[qc_summary_folder_path]}
         fi
        """

# Run only if summary files present
rule qc_report:
    input:
        mzidsummary = config["qc_summary_folder_path"]+"{s}.mzidsummary.txt"
    shell:
        """
        set -e
        if python3 DIA-or-DDA.py {input.featurescsv}; then
            echo dia;
            touch {output};
        else
        java -cp {config[qc_lcms_jar_path]} se.nbis.omics.QCMerge \
          {config[qc_summary_folder_path]}
        fi
       """


############---------- DIA ----------############
# DIA-NN
rule Diann:
    input:
        fasta = config["fasta"],
        mzML = config["path_to_mzml"] + "{sample}.mzML",
        featurescsv = config["path_to_csv_output"] + "{sample}.features.csv"
    output:
        result = config["path_to_result_output"] + "{sample}_result.tsv",
        tsv = config["path_to_tsv_output"] + "{sample}.tsv",
        statsfile = config["path_to_result_output"] + "{sample}.stats.tsv"
    shell:
    #library/50ngHeLas_lib.tsv.speclib
    #library/report-lib.predicted.speclib
        """
                set -e
        if python3 DIA-or-DDA.py {input.featurescsv}; then
            /usr/bin/time -v diann-1.8.1 --f {input.mzML}  --lib library/report-lib.predicted.speclib --threads 4 --verbose 1 --out {output.tsv} \
            --qvalue 0.01 --matrices  --out-lib {output.result} --gen-spec-lib --fasta {input.fasta} \
            --met-excision --cut K*,R* --relaxed-prot-inf --smart-profiling --peak-center --no-ifs-removal 

        else
            echo dda;
            touch {output};
        fi
        """

#Report
rule Report:
    input:
        mzmlcsv = config["path_to_csv_output"] + "{sample}.features.csv",
        dianntsv = config["path_to_result_output"] + "{sample}.stats.tsv",
        dinosaurtsv = config["dinosaur_output_path"]+"{sample}.features.tsv",
        diannresulttsv = config["path_to_result_output"]+"{sample}_result.tsv",
        diannbigtsv= config["path_to_tsv_output"]+"{sample}.tsv",
        xml=config["peptide_id_output"]+"{sample}.mzid",
    output:
        pdf = config["path_to_report"]+"{sample}.pdf",
        table = config["path_to_report"]+"{sample}.csv"
    shell:
        """
        set -e
        if python3 DIA-or-DDA.py {input.mzmlcsv}; then
            python3 DIA-Report.py {input.mzmlcsv} {input.dianntsv} {input.dinosaurtsv} {input.diannresulttsv} {input.diannbigtsv} {output.pdf} {output.table}

        else
            python3 DDA-Report.py {input.mzmlcsv} {input.dinosaurtsv} {input.xml} {output.pdf} {output.table} 
        fi
        """
        
#Comparison Report
def get_dianntsv_inputs(wildcards):
    return glob.glob(config["path_to_result_output"] + "*.stats.tsv")

def get_mzmlcsv_inputs(wildcards):
    return glob.glob(config["path_to_csv_output"] + "*.features.csv")

rule ComparisonReport:
    input:
        csv = expand(config["path_to_csv_output"] + "{sample}.features.csv", sample= samples_to_analyse),
        tsv = expand(config["path_to_result_output"] + "{sample}.stats.tsv", sample= samples_to_analyse),
        mzidsummary = expand(config["qc_summary_folder_path"]+"{sample}.mzidsummary.txt", sample= samples_to_analyse),
    output:
        pdf = config["path_to_report"]+ "Comparison.pdf"
    shell:
        "python3 Compare.py {input.csv} {input.tsv} {input.mzidsummary} {output.pdf}"

rule cleanup_empty_outputs:
    input:
        expand(config["peptide_id_output"]+"{sample}.mzid", sample= samples_to_analyse),
        expand(config["qc_summary_folder_path"]+"{sample}.mzidsummary.txt", sample= samples_to_analyse),
        expand(config["path_to_result_output"] + "{sample}_result.tsv", sample= samples_to_analyse),
        expand(config["path_to_tsv_output"] + "{sample}.tsv", sample= samples_to_analyse),
        expand(config["path_to_result_output"] + "{sample}.stats.tsv", sample= samples_to_analyse),
        pdf1 = config["path_to_report"]+ "Comparison.pdf",
        pdf = expand(config["path_to_report"]+"{sample}.pdf", sample= samples_to_analyse),
    output:
        clean = config["path_to_report"]+ ".cleanup_complete"
    shell:
        """
        python3 cleanup_empty_files.py {config[path_to_tsv_output]} {config[path_to_result_output]} {config[peptide_id_output]} {config[qc_summary_folder_path]}
        touch {output.clean}
        """

