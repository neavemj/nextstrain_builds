SEGMENTS = ["HA"]
FILTERING = ["ALL", "10-per-year"]

rule all:
    input:
        auspice_tree = expand("auspice/H9N2_{filter}_tree.json", filter=FILTERING),
        auspice_meta = expand("auspice/H9N2_{filter}_meta.json", filter=FILTERING)


rule files:
    params:
        raw_fasta = "raw_data/FilteredH9HAaligned_FWv1.fasta",
        geo_synonyms = "config/geo_synonyms.tsv",
        dropped_strains = "",
        reference = "config/GCF_000851145.1.gb",
        #colors = "config/colors.tsv",
        auspice_config = "config/auspice_config.json"

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
        sequences = "data/H9N2_{filter}.fasta",
        meta = "data/H9N2_{filter}.meta.tsv"
    log:
        "logs/{filter}.create_meta.log"
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
        sequences = "results/filtered_{filter}.fasta"
    run:
        if wildcards.filter == "ALL":
            shell(
            """
            augur filter \
                --sequences {input.sequences} \
                --metadata {input.metadata} \
                --output {output.sequences}
            """)
        elif wildcards.filter == "10-per-year":
            shell(
            """
            augur filter \
                --sequences {input.sequences} \
                --metadata {input.metadata} \
                --output {output.sequences} \
                --group-by country date \
                --sequences-per-group 10
            """)

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
        alignment = "results/aligned_{filter}.fasta"
    threads:
        32
    shell:
        """
        augur align \
            --sequences {input.sequences} \
            --reference-sequence {input.reference} \
            --output {output.alignment} \
            --nthreads {threads} \
            --remove-reference \
            --fill-gaps
        """

rule tree:
    message: "Building tree"
    input:
        alignment = rules.align.output.alignment
    output:
        "results/tree_raw_{filter}.nwk"
    params:
        method = "iqtree"
    threads:
        32
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
        tree = "results/tree_{filter}.nwk",
        node_data = "results/branch_lengths_{filter}.json"
    params:
        coalescent = "const",
        date_inference = "marginal",
        clock_filter_iqd = 4
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
            --date-inference {params.date_inference} \
            --clock-filter-iqd {params.clock_filter_iqd}
        """

rule ancestral:
    message: "Reconstructing ancestral sequences and mutations"
    input:
        tree = rules.refine.output.tree,
        alignment = rules.align.output
    output:
        node_data = "results/nt_muts_{filter}.json"
    params:
        inference = "joint"
    shell:
        """
        augur ancestral \
            --tree {input.tree} \
            --alignment {input.alignment} \
            --output {output.node_data} \
            --inference {params.inference}
        """

rule translate:
    message: "Translating amino acid sequences"
    input:
        tree = rules.refine.output.tree,
        node_data = rules.ancestral.output.node_data,
        reference = files.reference
    output:
        node_data = "results/aa_muts_{filter}.json"
    shell:
        """
        augur translate \
            --tree {input.tree} \
            --ancestral-sequences {input.node_data} \
            --reference-sequence {input.reference} \
            --output {output.node_data}
        """

rule traits:
    message: "Inferring ancestral traits for {params.columns!s}"
    input:
        tree = rules.refine.output.tree,
        metadata = rules.create_meta.output.meta
    output:
        node_data = "results/traits_{filter}.json"
    params:
        columns = "location country"
    shell:
        """
        augur traits \
            --tree {input.tree} \
            --metadata {input.metadata} \
            --output {output.node_data} \
            --columns {params.columns} \
            --confidence
        """

rule export:
    message: "Exporting data files for for auspice"
    input:
        tree = rules.refine.output.tree,
        metadata = rules.create_meta.output.meta,
        branch_lengths = rules.refine.output.node_data,
        traits = rules.traits.output.node_data,
        nt_muts = rules.ancestral.output.node_data,
        aa_muts = rules.translate.output.node_data,
        lat_longs = "config/lat_longs.tsv",
        auspice_config = files.auspice_config
    output:
        auspice_tree = "auspice/H9N2_{filter}_tree.json",
        auspice_meta = "auspice/H9N2_{filter}_meta.json"
    shell:
        """
        augur export \
            --tree {input.tree} \
            --metadata {input.metadata} \
            --node-data {input.branch_lengths} {input.traits} {input.nt_muts} {input.aa_muts} \
            --auspice-config {input.auspice_config} \
            --lat-longs {input.lat_longs} \
            --output-tree {output.auspice_tree} \
            --output-meta {output.auspice_meta}
        """
