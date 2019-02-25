SEGMENTS = ["HA"]

rule all:
    input:
        expand("results/tree_{segment}.nwk", segment=SEGMENTS)
        #auspice_tree = expand("auspice/H9N2_{segment}_tree.json", segment=SEGMENTS),
        #auspice_meta = expand("auspice/H9N2_{segment}_meta.json", segment=SEGMENTS)


rule files:
    params:
        raw_fasta = "raw_data/FilteredH9HAaligned_FWv1.fasta",
        geo_synonyms = "config/geo_synonyms.tsv",
        dropped_strains = "",
        reference = "config/GCF_000851145.1.gb"
        #colors = "config/colors.tsv",
        #auspice_config = "config/auspice_config.json"

files = rules.files.params


rule create_meta:
    message:
        """
        take raw H9N2 sequences after Ivano's pipeline
        will fix headers and extract required metadata
        sequences not meeting requirements will be filtered
        a geo_synonyms file is required to check country
        """
    input:
        sequences = files.raw_fasta,
        geo_synonyms = files.geo_synonyms
    output:
        sequences = "data/H9N2_{segment}.fasta",
        meta = "data/H9N2_{segment}.meta.tsv"
    log:
        "logs/{segment}.create_meta.log"
    shell:
        """
        ./scripts/fix_h9n2_headers.py \
            --fasta_file {input.sequences} \
            --geo {input.geo_synonyms} \
            --fasta_output {output.sequences} \
            --meta_output {output.meta} > {log}
        """

rule filter:
    message:
        """
        Filtering to sequences
        """
    input:
        sequences = rules.create_meta.output.sequences,
        metadata = rules.create_meta.output.meta
    output:
        sequences = "results/filtered_{segment}.fasta"
    shell:
        """
        augur filter \
            --sequences {input.sequences} \
            --metadata {input.metadata} \
            --output {output.sequences} \
        """

rule align:
    message:
        """
        Aligning sequences to {input.reference}
          - filling gaps with N
        """
    input:
        sequences = rules.filter.output.sequences,
        reference = files.reference
    output:
        alignment = "results/aligned_{segment}.fasta"
    threads:
        8
    shell:
        """
        augur align \
            --sequences {input.sequences} \
            --reference-sequence {input.reference} \
            --output {output.alignment} \
            --nthreads {threads} \
            --fill-gaps
        """

rule tree:
    message: "Building tree"
    input:
        alignment = rules.align.output.alignment
    output:
        "results/tree_raw_{segment}.nwk"
    params:
        method = "iqtree"
    threads:
        8
    shell:
        """
        augur tree \
            --alignment {input.alignment} \
            --output {output} \
            --method {params.method} \
            --nthreads {threads}
        """

rule refine:
    message:
        """
        Refining tree
          - estimate timetree
          - use {params.coalescent} coalescent timescale
          - estimate {params.date_inference} node dates
        """
    input:
        tree = rules.tree.output,
        alignment = rules.align.output,
        metadata = rules.create_meta.output.meta
    output:
        tree = "results/tree_{segment}.nwk",
        node_data = "results/branch_lengths_{segment}.json"
    params:
        coalescent = "opt",
        date_inference = "marginal"
    shell:
        """
        augur refine \
            --tree {input.tree} \
            --alignment {input.alignment} \
            --metadata {input.metadata} \
            --output-tree {output.tree} \
            --output-node-data {output.node_data} \
            --timetree \
            --coalescent {params.coalescent} \
            --date-confidence \
            --date-inference {params.date_inference}
        """
