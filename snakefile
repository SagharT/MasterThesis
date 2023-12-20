'''
This snake fike automatically picks up new .raw files from the specified directory for processing.
Converts .raw files to .mzML format using the msconvert tool.
Parses the .mzML files to extract specific features and saves this data into CSV files.
'''
configfile: "config.yaml"

import os

# Generate the list of samples to analyse
samples_to_analyse = [sample.replace('.raw', '') for sample in os.listdir(config["path_to_samples"]) if sample.endswith(".raw")]

rule all:
    input:
        expand(config["path_to_csv_output"] + "{sample}_features.csv", sample=samples_to_analyse),
        expand(config["qc_summary_folder_path"]+"{s}.mzmlsummary.txt", s=samples_to_analyse),
        expand(config["qc_summary_folder_path"]+"{s}.mzidsummary.txt", s=samples_to_analyse), #peptide identification (MSGFPlus)
        expand(config["qc_summary_folder_path"]+"{s}.featuresummary.txt", s=samples_to_analyse), #feature detection (Dinosaur)
        expand(config["path_to_tsv_output"] + "{sample}_result.tsv", sample=samples_to_analyse),
        expand(config["path_to_tsv_output"] + "{sample}.tsv", sample=samples_to_analyse)
        

# Raw file conversion
rule convert_raw_file:
    input:
        raw = config["path_to_samples"] + "{sample}.raw"
    params:
        rawfile = "{sample}.raw"
    output:
        mzML = temp(config["path_to_mzml"] + "{sample}.mzML")
    shell:
        "docker run --rm -v {config[path_to_samples]}:/indata "
        "-v {config[path_to_mzml]}:/outdata chambm/pwiz-skyline-i-agree-to-the-vendor-licenses wine msconvert "
        "-o /outdata -v --filter=\"peakPicking true 2-\" /indata/{params.rawfile};"
        "sudo chmod -R ou+w {config[path_to_mzml]}"

# Peptide identification
rule peptide_identification:
    input:
        mzml = config["path_to_mzml"]+"{s}.mzML",
        db = config["fasta"]
    params:
        sample = "{s}",
    output:
        config["peptide_id_output"]+"{s}.mzid"
    shell:
       """
       java {config[max_heap_size]} -Djava.awt.headless=true \
         -jar {config[msgfplus_jar_path]} \
         {config[msgfplus_extra]} \
         -mod {config[msgfplus_mod_path]} \
         -d {input.db} \
         -o {output} \
         -s {input.mzml}
       """
# Feature extraction, Dinosaur
rule feat_extraction:
    input:
        config["path_to_mzml"]+"{s}.mzML"
    params:
        qc = config["dinosaur_output_path"]+"{s}.qc.zip"
    output:
        config["dinosaur_output_path"]+"{s}.features.tsv"
    shell:
        """
        java {config[max_heap_size]} -Djava.awt.headless=true \
          -jar {config[dinosaur_jar_path]} \
          --zipQcFolder=true \
          {config[dinosaur_params]} \
          --outDir={config[dinosaur_output_path]} \
          {input}
        """
# Quality control summaries

rule qc_mzML_summary:
    input:
        config["path_to_mzml"]+"{s}.mzML"
    output:
        config["qc_summary_folder_path"]+"{s}.mzmlsummary.txt"
    shell:
        """
        java -cp {config[qc_lcms_jar_path]} \
          se.nbis.omics.MzMLSummary \
          {input} \
          >{output}
        """

rule qc_feature_summary:
    input:
        config["dinosaur_output_path"]+"{s}.features.tsv"
    output:
        config["qc_summary_folder_path"]+"{s}.featuresummary.txt"
    shell:
        """
        java -cp {config[qc_lcms_jar_path]} \
          se.nbis.omics.FeaturesQCSummary \
          {input} \
          >{output}
        """

rule qc_peptideid_mzid_summary:
    input:
        config["peptide_id_output"]+"{s}.mzid"
    output:
        config["qc_summary_folder_path"]+"{s}.mzidsummary.txt"
    shell:
        """
        java -cp {config[qc_lcms_jar_path]} \
          se.nbis.omics.MzIdentMLSummary \
          {input} \
          {config[qc_fdr]} \
          >{output}

        java -cp {config[qc_lcms_jar_path]} se.nbis.omics.QCMerge \
          {config[qc_summary_folder_path]}

        """

# Run only if summary files present

rule qc_report:
    input:
        config["qc_summary_folder_path"]+"{s}.mzidsummary.txt"
    shell:
        """
        java -cp {config[qc_lcms_jar_path]} se.nbis.omics.QCMerge \
          {config[qc_summary_folder_path]}
        """

# mzML file parsing
rule parse_mzML:
    input:
        mzML = config["path_to_mzml"] + "{sample}.mzML"
    output:
        csv = config["path_to_csv_output"] + "{sample}_features.csv"
    shell:
        "python3 mzML-Parser.py {input.mzML} {output.csv}"

# DIA-NN
rule Diann:
    input:
        fasta = config["fasta"],
        mzML = config["path_to_mzml"] + "{sample}.mzML"
    output:
        result = config["path_to_result_output"] + "{sample}_result.tsv",
        tsv = config["path_to_tsv_output"] + "{sample}.tsv"
    shell:
#        "diann-1.8.1 --f {input.mzML}  --lib \"\" --threads 1 --verbose 1 --out {output.tsv}"
#       " --qvalue 0.01 --matrices  --out-lib {output.result} --gen-spec-lib --predictor --fasta {input.fasta}"
#      " --fasta-search --min-fr-mz 200 --max-fr-mz 1800 --met-excision --cut K*,R* --missed-cleavages 1"
#        " --min-pep-len 7 --max-pep-len 30 --min-pr-mz 300 --max-pr-mz 1800 --min-pr-charge 1 --max-pr-charge 4"
#        " --unimod4 --var-mods 1 --var-mod UniMod:21,79.966331,STY --monitor-mod UniMod:21 --reanalyse"
#        " --relaxed-prot-inf --smart-profiling --peak-center"
        "diann-1.8.1 --f {input.mzML}  --lib library/50ngHeLas_lib.tsv.speclib --threads 4 --verbose 1 --out {output.tsv}"
        " --qvalue 0.01 --matrices  --out-lib {output.result} --gen-spec-lib --fasta {input.fasta}"
        " --met-excision --cut K*,R* --relaxed-prot-inf --smart-profiling --peak-center --no-ifs-removal "

