
ANALYSES = ["fillgaps"]
ORFS = ["ORF125"]

rule all:
    input:
        auspice_tree = expand("auspice/wssv_{analysis}_{orf}_tree.json", analysis=ANALYSES, orf=ORFS),
        auspice_meta = expand("auspice/wssv_{analysis}_{orf}_meta.json", analysis=ANALYSES, orf=ORFS)

rule files:
    params:
        input_alignment = "results/aligned_{analysis}.fasta",
        input_meta = "data/wssv_{analysis}.meta.tsv",
        dropped_strains = "",
        reference = "config/whitespot_AF369029.2.gb",
        #colors = "config/colors.tsv",
        auspice_config = "config/auspice_config.json"

files = rules.files.params

rule chop_alignment:
    message: "Chopping the alignment at the selected ORF boundaries"
    input:
        alignment = files.input_alignment
    output:
        chopped_alignment = "results/aligned_{analysis}_{orf}.fasta"
    params:
        orf: "ORF125"
    shell:
        """
        chop_align.py \
            -i {input.alignment}
            -o {output.chopped_alignment}
        """

rule tree:
    message: "Building tree"
    input:
        alignment = rules.chop_alignment.output
    output:
        tree = "results/tree_raw_{analysis}_{orf}.nwk"
    params:
        method = "iqtree"
    shell:
        """
        augur tree \
            --alignment {input.alignment} \
            --output {output.tree} \
            --method {params.method}
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
        tree = rules.tree.output.tree,
        alignment = files.input_alignment,
        metadata = files.input_meta
    output:
        tree = "results/tree_{analysis}_{orf}.nwk",
        node_data = "results/branch_lengths_{analysis}_{orf}.json"
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
            --date-format %d/%m/%Y \
            --coalescent {params.coalescent} \
            --date-confidence \
            --date-inference {params.date_inference}
        """

rule ancestral:
    message: "Reconstructing ancestral sequences and mutations"
    input:
        tree = rules.refine.output.tree,
        alignment = files.input_alignment
    output:
        node_data = "results/nt_muts_{analysis}_{orf}.json"
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
        node_data = "results/aa_muts_{analysis}_{orf}.json"
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
        metadata = files.input_meta
    output:
        node_data = "results/traits_{analysis}_{orf}.json"
    params:
        columns = "country"
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
        metadata = files.input_meta,
        branch_lengths = rules.refine.output.node_data,
        traits = rules.traits.output.node_data,
        nt_muts = rules.ancestral.output.node_data,
        aa_muts = rules.translate.output.node_data,
        auspice_config = files.auspice_config,
        lat_longs = "data/lat_longs_aus.tsv"
    output:
        auspice_tree = "auspice/wssv_{analysis}_{orf}_tree.json",
        auspice_meta = "auspice/wssv_{analysis}_{orf}_meta.json"
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
