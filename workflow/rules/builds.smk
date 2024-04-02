"""Rules for generating simulated HA sequences for validation of forecasting models.
"""
BUILD_PATH = "results/builds/{type}/{sample}/"
BUILD_LOG_STEM = "{type}_{sample}"
BUILD_TIMEPOINT_PATH = BUILD_PATH + "timepoints/{timepoint}/"
BUILD_SEGMENT_LOG_STEM = "{type}_{sample}_{timepoint}"


rule get_strains_by_timepoint:
    input:
        metadata = _get_metadata_by_wildcards
    output:
        strains = BUILD_TIMEPOINT_PATH + "strains.txt"
    params:
        date_field = _get_date_field_to_partition_strains_by,
        years_back = _get_years_back_to_build_trees,
        reference_strains = _get_required_strains_argument
    conda: "../envs/anaconda.python3.yaml"
    shell:
        """
        python3 workflow/scripts/partition_strains_by_timepoint.py \
            {input.metadata} \
            {wildcards.timepoint} \
            {output} \
            --date-field {params.date_field} \
            --years-back {params.years_back} \
            {params.reference_strains}
        """


rule extract:
    input:
        sequences = _get_sequences_by_wildcards,
        strains = rules.get_strains_by_timepoint.output.strains,
    output:
        sequences = BUILD_TIMEPOINT_PATH + "filtered_sequences.fasta"
    conda: "../envs/anaconda.python3.yaml"
    shell:
        """
        python3 workflow/scripts/extract_sequences.py \
            --sequences {input.sequences} \
            --samples {input.strains} \
            --output {output}
        """


rule align:
    input:
        sequences = rules.extract.output.sequences,
        reference = _get_reference
    output:
        alignment = BUILD_TIMEPOINT_PATH + "aligned.fasta"
    conda: "../envs/anaconda.python3.yaml"
    benchmark: "benchmarks/align_" + BUILD_SEGMENT_LOG_STEM + ".txt"
    threads: 4
    resources:
        mem_mb=16000,
        runtime="0:30:00",
        partition="campus-new",
        qos="campus-new",
    shell:
        """
        augur align \
            --sequences {input.sequences} \
            --reference-sequence {input.reference} \
            --output {output.alignment} \
            --remove-reference \
            --fill-gaps \
            --nthreads {threads}
        """


def _get_alignment_by_wildcards(wildcards):
    """Return filtered sequences for simulated builds, since these are already
    aligned. Otherwise, return the name of the multiple sequence alignment
    output.

    """
    if wildcards.type == "simulated":
        return BUILD_TIMEPOINT_PATH + "filtered_sequences.fasta"
    else:
        return BUILD_TIMEPOINT_PATH + "aligned.fasta"


rule tree:
    message: "Building tree ({wildcards})"
    input:
        alignment = _get_alignment_by_wildcards,
    output:
        tree = BUILD_TIMEPOINT_PATH + "tree_raw.nwk"
    conda: "../envs/anaconda.python3.yaml"
    shadow: "minimal"
    benchmark: "benchmarks/tree_" + BUILD_SEGMENT_LOG_STEM + ".txt"
    log: "logs/tree_" + BUILD_SEGMENT_LOG_STEM + ".log"
    threads: 4
    resources:
        mem_mb=16000,
        runtime="2:00:00",
        partition="campus-new",
        qos="campus-new",
    shell:
        """
        augur tree \
            --alignment {input.alignment} \
            --output {output.tree} \
            --method iqtree \
            --nthreads {threads} &> {log}
        """


rule refine:
    input:
        tree = rules.tree.output.tree,
        alignment = _get_alignment_by_wildcards,
        metadata = _get_metadata_by_wildcards
    output:
        tree = BUILD_TIMEPOINT_PATH + "tree.nwk",
        node_data = BUILD_TIMEPOINT_PATH + "branch_lengths.json"
    params:
        coalescent = "const",
        date_inference = "marginal",
        clock_rate = _get_clock_rate_argument,
        clock_std_dev = _get_clock_std_dev_argument
    conda: "../envs/anaconda.python3.yaml"
    benchmark: "benchmarks/refine_" + BUILD_SEGMENT_LOG_STEM + ".txt"
    log: "logs/refine_" + BUILD_SEGMENT_LOG_STEM + ".log"
    resources:
        mem_mb=16000,
        runtime="4:00:00",
        partition="campus-new",
        qos="campus-new",
    shell:
        """
        augur refine \
            --tree {input.tree} \
            --alignment {input.alignment} \
            --metadata {input.metadata} \
            --output-tree {output.tree} \
            --output-node-data {output.node_data} \
            --timetree \
            --no-covariance \
            --use-fft \
            --stochastic-resolve \
            {params.clock_rate} \
            {params.clock_std_dev} \
            --coalescent {params.coalescent} \
            --date-confidence \
            --date-inference {params.date_inference} &> {log}
        """


rule tip_frequencies:
    message:
        """
        Estimating tip frequencies for {input.tree}
          - narrow bandwidth: {params.narrow_bandwidth}
          - wide bandwidth: {params.wide_bandwidth}
          - proportion wide: {params.proportion_wide}
        """
    input:
        tree = rules.refine.output.tree,
        metadata = _get_metadata_by_wildcards
    output:
        frequencies = "results/auspice/flu_" + BUILD_SEGMENT_LOG_STEM + "_original-tip-frequencies.json"
    params:
        narrow_bandwidth=config["frequencies"]["narrow_bandwidth"],
        wide_bandwidth=config["frequencies"]["wide_bandwidth"],
        proportion_wide=config["frequencies"]["proportion_wide"],
        pivot_frequency=_get_pivot_interval,
        min_date=_get_min_date_for_augur_frequencies_by_wildcards,
        max_date=_get_max_date_for_augur_frequencies_by_wildcards
    conda: "../envs/anaconda.python3.yaml"
    benchmark: "benchmarks/tip_frequencies_" + BUILD_SEGMENT_LOG_STEM + ".txt"
    log: "logs/tip_frequencies_" + BUILD_SEGMENT_LOG_STEM + ".log"
    shell:
        """
        augur frequencies \
            --method kde \
            --tree {input.tree} \
            --metadata {input.metadata} \
            --narrow-bandwidth {params.narrow_bandwidth} \
            --wide-bandwidth {params.wide_bandwidth} \
            --proportion-wide {params.proportion_wide} \
            --min-date {params.min_date} \
            --max-date {params.max_date} \
            --pivot-interval {params.pivot_frequency} \
            --output {output}
        """

rule convert_frequencies_to_table:
    input:
        tree = rules.refine.output.tree,
        frequencies = rules.tip_frequencies.output.frequencies
    output:
        table = BUILD_TIMEPOINT_PATH + "frequencies.tsv"
    params:
        method = "kde"
    conda: "../envs/anaconda.python3.yaml"
    shell:
        """
        python3 workflow/scripts/frequencies_to_table.py \
            --tree {input.tree} \
            --frequencies {input.frequencies} \
            --method {params.method} \
            --output {output} \
            --annotations timepoint={wildcards.timepoint}
        """


rule ancestral:
    message: "Reconstructing ancestral sequences and mutations for {wildcards}"
    input:
        tree = rules.refine.output.tree,
        alignment = _get_alignment_by_wildcards,
    output:
        node_data = BUILD_TIMEPOINT_PATH + "nt_muts.json"
    params:
        inference = "joint"
    conda: "../envs/anaconda.python3.yaml"
    benchmark: "benchmarks/ancestral_" + BUILD_SEGMENT_LOG_STEM + ".txt"
    log: "logs/ancestral_" + BUILD_SEGMENT_LOG_STEM + ".txt"
    resources:
        mem_mb=4000,
        runtime="0:15:00",
    shell:
        """
        augur ancestral \
            --tree {input.tree} \
            --alignment {input.alignment} \
            --output-node-data {output.node_data} \
            --inference {params.inference} &> {log}
        """


rule translate:
    message: "Translating amino acid sequences"
    input:
        tree = rules.refine.output.tree,
        node_data = rules.ancestral.output.node_data,
        reference = _get_reference
    output:
        node_data = BUILD_TIMEPOINT_PATH + "aa_muts.json"
    conda: "../envs/anaconda.python3.yaml"
    benchmark: "benchmarks/translate_" + BUILD_SEGMENT_LOG_STEM + ".txt"
    log: "logs/translate_" + BUILD_SEGMENT_LOG_STEM + ".txt"
    shell:
        """
        augur translate \
            --tree {input.tree} \
            --ancestral-sequences {input.node_data} \
            --reference-sequence {input.reference} \
            --output-node-data {output.node_data} &> {log}
        """


rule reconstruct_translations:
    message: "Reconstructing translations for {wildcards.gene}"
    input:
        tree = rules.refine.output.tree,
        node_data = rules.translate.output.node_data
    output:
        aa_alignment = BUILD_TIMEPOINT_PATH + "aa-seq_{gene}.fasta"
    conda: "../envs/anaconda.python3.yaml"
    benchmark: "benchmarks/reconstruct_translations_{gene}_" + BUILD_SEGMENT_LOG_STEM + ".txt"
    log: "logs/reconstruct_translations_{gene}_" + BUILD_SEGMENT_LOG_STEM + ".txt"
    resources:
        mem_mb=4000,
    shell:
        """
        augur reconstruct-sequences \
            --tree {input.tree} \
            --mutations {input.node_data} \
            --gene {wildcards.gene} \
            --output {output.aa_alignment} \
            --internal-nodes &> {log}
        """


rule convert_translations_to_json:
    input:
        tree = rules.refine.output.tree,
        translations = translations(segment="ha", path=BUILD_TIMEPOINT_PATH)
    output:
        translations = BUILD_TIMEPOINT_PATH + "aa_seq.json"
    params:
        gene_names = gene_names(segment="ha")
    conda: "../envs/anaconda.python3.yaml"
    shell:
        """
        python3 workflow/scripts/convert_translations_to_json.py \
            --tree {input.tree} \
            --alignment {input.translations} \
            --gene-names {params.gene_names} \
            --output {output.translations}
        """


rule traits:
    message: "Inferring ancestral traits for {params.columns!s}"
    input:
        tree = rules.refine.output.tree,
        metadata = _get_metadata_by_wildcards
    output:
        node_data = BUILD_TIMEPOINT_PATH + "traits.json",
    params:
        columns = "region country"
    conda: "../envs/anaconda.python3.yaml"
    benchmark: "benchmarks/traits_" + BUILD_SEGMENT_LOG_STEM + ".txt"
    shell:
        """
        augur traits \
            --tree {input.tree} \
            --metadata {input.metadata} \
            --output {output.node_data} \
            --columns {params.columns} \
            --confidence
        """


rule distances:
    input:
        tree = rules.refine.output.tree,
        alignments = translations(segment="ha", path=BUILD_TIMEPOINT_PATH),
        # TODO: define distance maps in build configs
        distance_maps = _get_distance_maps_for_simulations,
        date_annotations = rules.refine.output.node_data
    params:
        genes = gene_names(segment="ha"),
        comparisons = _get_distance_comparisons_for_simulations,
        attribute_names = _get_distance_attributes_for_simulations,
        earliest_date = _get_distance_earliest_date_by_wildcards,
        latest_date = _get_distance_latest_date_by_wildcards
    output:
        distances = BUILD_TIMEPOINT_PATH + "distances.json",
    conda: "../envs/anaconda.python3.yaml"
    resources:
        mem_mb=8000,
        runtime="00:30:00",
    shell:
        """
        augur distance \
            --tree {input.tree} \
            --alignment {input.alignments} \
            --gene-names {params.genes} \
            --compare-to {params.comparisons} \
            --attribute-name {params.attribute_names} \
            --map {input.distance_maps} \
            --date-annotations {input.date_annotations} \
            --earliest-date {params.earliest_date} \
            --latest-date {params.latest_date} \
            --output {output}
        """


rule lbi:
    message: "Calculating LBI"
    input:
        tree = rules.refine.output.tree,
        branch_lengths = rules.refine.output.node_data
    params:
        tau = config["lbi"]["tau"],
        window = config["lbi"]["window"],
        names = "lbi"
    output:
        lbi = BUILD_TIMEPOINT_PATH + "lbi.json"
    conda: "../envs/anaconda.python3.yaml"
    shell:
        """
        augur lbi \
            --tree {input.tree} \
            --branch-lengths {input.branch_lengths} \
            --output {output} \
            --attribute-names {params.names} \
            --tau {params.tau} \
            --window {params.window}
        """


rule unnormalized_lbi:
    input:
        tree = rules.refine.output.tree,
        branch_lengths = rules.refine.output.node_data
    params:
        tau = config["lbi"]["tau"],
        window = config["lbi"]["window"],
        names = "unnormalized_lbi"
    output:
        lbi = BUILD_TIMEPOINT_PATH + "unnormalized_lbi.json"
    conda: "../envs/anaconda.python3.yaml"
    shell:
        """
        augur lbi \
            --tree {input.tree} \
            --branch-lengths {input.branch_lengths} \
            --output {output} \
            --attribute-names {params.names} \
            --tau {params.tau} \
            --window {params.window} \
            --no-normalization
        """


rule normalize_fitness:
    input:
        metadata = _get_metadata_by_wildcards,
        frequencies = rules.convert_frequencies_to_table.output.table
    output:
        fitness = BUILD_TIMEPOINT_PATH + "normalized_fitness.json"
    params:
        preferred_frequency_method = config["frequencies"]["preferred_method"]
    conda: "../envs/anaconda.python3.yaml"
    shell:
        """
        python3 workflow/scripts/normalize_fitness.py \
            --metadata {input.metadata} \
            --frequencies-table {input.frequencies} \
            --frequency-method {params.preferred_frequency_method} \
            --output {output.fitness}
        """


def _get_node_data_for_export(wildcards):
    """Return a list of node data files to include for a given build's wildcards.
    """
    # Define inputs shared by specific builds.
    inputs = [
        rules.refine.output.node_data,
        rules.ancestral.output.node_data,
        rules.translate.output.node_data,
        rules.convert_translations_to_json.output.translations,
        rules.distances.output.distances,
        rules.lbi.output.lbi,
    ]

    # Define node data that only make sense for natural populations
    # such as titer models.
    if wildcards.type == "natural":
        inputs.extend([
            rules.traits.output.node_data,
        ])
    elif wildcards.type == "simulated":
        inputs.extend([
            rules.normalize_fitness.output.fitness
        ])

    # Convert input files from wildcard strings to real file names.
    inputs = [input_file.format(**wildcards) for input_file in inputs]
    return inputs


rule convert_node_data_to_table:
    input:
        tree = rules.refine.output.tree,
        metadata = _get_metadata_by_wildcards,
        node_data = _get_node_data_for_export
    output:
        table = BUILD_TIMEPOINT_PATH + "node_data.tsv"
    params:
        excluded_fields_arg = _get_excluded_fields_arg,
        lineage = _get_lineage,
        segment = _get_segment
    conda: "../envs/anaconda.python3.yaml"
    shell:
        """
        python3 workflow/scripts/node_data_to_table.py \
            --tree {input.tree} \
            --metadata {input.metadata} \
            --jsons {input.node_data} \
            --output {output} \
            {params.excluded_fields_arg} \
            --annotations timepoint={wildcards.timepoint} \
                          lineage={params.lineage} \
                          segment={params.segment}
        """


rule merge_node_data_and_frequencies:
    input:
        node_data = rules.convert_node_data_to_table.output.table,
        kde_frequencies = rules.convert_frequencies_to_table.output.table,
    output:
        table = BUILD_TIMEPOINT_PATH + "tip_attributes.tsv"
    params:
        preferred_frequency_method = config["frequencies"]["preferred_method"]
    conda: "../envs/anaconda.python3.yaml"
    shell:
        """
        python3 workflow/scripts/merge_node_data_and_frequencies.py \
            --node-data {input.node_data} \
            --kde-frequencies {input.kde_frequencies} \
            --preferred-frequency-method {params.preferred_frequency_method} \
            --output {output.table}
        """


rule collect_tip_attributes:
    input:
        _get_tip_attributes_by_wildcards
    output:
        attributes = BUILD_PATH + "tip_attributes.tsv"
    conda: "../envs/anaconda.python3.yaml"
    shell:
        """
        python3 workflow/scripts/collect_tables.py \
            --tables {input} \
            --output {output.attributes}
        """


rule annotate_sample_and_delay_type_to_tip_attributes:
    input:
        attributes = BUILD_PATH + "tip_attributes.tsv"
    output:
        attributes = BUILD_PATH + "tip_attributes_with_sample_and_delay_type.tsv"
    conda: "../envs/csv.yaml"
    params:
        columns="timepoint,strain,date,frequency,aa_sequence",
        delay_type=lambda wildcards: config["builds"].get(wildcards.type, {}).get(wildcards.sample, {}).get("delay_type"),
    shell:
        """
        csvtk cut -t -f {params.columns} {input.attributes} \
            | csvtk mutate2 -t -n sample -e "'{wildcards.sample}'" \
            | csvtk mutate2 -t -n delay_type -e "'{params.delay_type}'" > {output}
        """


rule annotate_naive_tip_attribute:
    input:
        attributes = rules.collect_tip_attributes.output.attributes
    output:
        attributes = BUILD_PATH + "tip_attributes_with_naive_predictor.tsv"
    conda: "../envs/anaconda.python3.yaml"
    shell:
        """
        python3 workflow/scripts/annotate_naive_tip_attribute.py \
            --tip-attributes {input.attributes} \
            --output {output.attributes}
        """


rule target_distances:
    input:
        attributes = rules.annotate_naive_tip_attribute.output.attributes
    output:
        distances = BUILD_PATH + "target_distances.tsv",
    params:
        delta_months = config["fitness_model"]["delta_months_to_fit"],
        sequence_attribute_name = "aa_sequence"
    benchmark: "benchmarks/target_distances_" + BUILD_LOG_STEM + ".txt"
    conda: "../envs/anaconda.python3.yaml"
    shell:
        """
        python3 workflow/scripts/calculate_target_distances.py \
            --tip-attributes {input.attributes} \
            --delta-months {params.delta_months} \
            --sequence-attribute-name {params.sequence_attribute_name} \
            --output {output}
        """


rule annotate_weighted_distances_for_tip_attributes:
    input:
        attributes = rules.annotate_naive_tip_attribute.output.attributes,
        distances = rules.target_distances.output.distances
    output:
        attributes = BUILD_PATH + "tip_attributes_with_weighted_distances.tsv"
    params:
        delta_months = config["fitness_model"]["delta_months_to_fit"]
    conda: "../envs/anaconda.python3.yaml"
    shell:
        """
        python3 src/weighted_distances.py \
            --tip-attributes {input.attributes} \
            --distances {input.distances} \
            --delta-months {params.delta_months} \
            --output {output}
        """


rule fit_models_by_distances:
    input:
        attributes = rules.annotate_weighted_distances_for_tip_attributes.output.attributes,
        distances = rules.target_distances.output.distances
    output:
        model = BUILD_PATH + "models_by_distances/{predictors}.json",
        errors = BUILD_PATH + "models_by_distances_errors/{predictors}.tsv",
        coefficients = BUILD_PATH + "models_by_distances_coefficients/{predictors}.tsv"
    params:
        predictors = _get_predictor_list,
        delta_months = config["fitness_model"]["delta_months_to_fit"],
        training_window = _get_fitness_model_training_window,
        cost_function = config["fitness_model"]["distance_cost_function"],
        l1_lambda = config["fitness_model"]["l1_lambda"]
    conda: "../envs/anaconda.python3.yaml"
    benchmark: "benchmarks/fitness_model_distances_" + BUILD_LOG_STEM + "_{predictors}.txt"
    log: "logs/fitness_model_distances_" + BUILD_LOG_STEM + "_{predictors}.txt"
    resources:
        mem_mb=20000,
        runtime="06:00:00",
        partition="campus-new",
        qos="campus-new",
    shell:
        """
        python3 src/fit_model.py \
            --tip-attributes {input.attributes} \
            --training-window {params.training_window} \
            --delta-months {params.delta_months} \
            --predictors {params.predictors} \
            --cost-function {params.cost_function} \
            --l1-lambda {params.l1_lambda} \
            --target distances \
            --distances {input.distances} \
            --errors-by-timepoint {output.errors} \
            --coefficients-by-timepoint {output.coefficients} \
            --include-scores \
            --output {output.model} &> {log}
        """


rule extract_minimal_models_by_distances:
    input:
        model = rules.fit_models_by_distances.output.model
    output:
        model = BUILD_PATH + "minimal_models_by_distances/{predictors}.json",
    conda: "../envs/anaconda.python3.yaml"
    shell:
        """
        python3 workflow/scripts/extract_minimal_models_by_distances.py \
            --model {input.model} \
            --output {output.model}
        """


rule annotate_distance_models:
    input:
        attributes = rules.annotate_weighted_distances_for_tip_attributes.output.attributes,
        model = rules.fit_models_by_distances.output.model,
        errors = rules.fit_models_by_distances.output.errors,
        coefficients = rules.fit_models_by_distances.output.coefficients
    output:
        errors = BUILD_PATH + "annotated_models_by_distances_errors/{predictors}.tsv",
        coefficients = BUILD_PATH + "annotated_models_by_distances_coefficients/{predictors}.tsv"
    params:
        delta_months = config["fitness_model"]["delta_months_to_fit"],
        error_type = "validation"
    conda: "../envs/anaconda.python3.yaml"
    shell:
        """
        python3 workflow/scripts/annotate_model_tables.py \
            --tip-attributes {input.attributes} \
            --model {input.model} \
            --errors-by-timepoint {input.errors} \
            --coefficients-by-timepoint {input.coefficients} \
            --annotated-errors-by-timepoint {output.errors} \
            --annotated-coefficients-by-timepoint {output.coefficients} \
            --delta-months {params.delta_months} \
            --annotations type="{wildcards.type}" sample="{wildcards.sample}" error_type="{params.error_type}"
        """


rule plot_tree:
    input:
        tree = rules.refine.output.tree
    output:
        tree = BUILD_TIMEPOINT_PATH + "tree.pdf"
    conda: "../envs/anaconda.python3.yaml"
    benchmark: "benchmarks/plot_tree_" + BUILD_SEGMENT_LOG_STEM + ".txt"
    log: "logs/plot_tree_" + BUILD_SEGMENT_LOG_STEM + ".log"
    shell:
        """
        python3 workflow/scripts/plot_tree.py {input} {output} &> {log}
        """


rule aggregate_tree_plots:
    input: _get_tree_plots_by_wildcards
    output:
        trees="results/figures/trees_" + BUILD_LOG_STEM + ".pdf"
    conda: "../envs/anaconda.python3.yaml"
    shell: "gs -dBATCH -dNOPAUSE -q -sDEVICE=pdfwrite -sOutputFile={output} {input}"


rule target_distances_by_timepoint:
    input:
        attributes = rules.merge_node_data_and_frequencies.output.table
    output:
        distances = BUILD_TIMEPOINT_PATH + "target_distances.tsv",
    params:
        delta_months = _get_delta_months_to_forecast
    conda: "../envs/anaconda.python3.yaml"
    shell:
        """
        python3 workflow/scripts/calculate_target_distances.py \
            --tip-attributes {input.attributes} \
            --delta-months {params.delta_months} \
            --output {output}
        """


rule forecast_tips:
    input:
        attributes = rules.merge_node_data_and_frequencies.output.table,
        frequencies = rules.tip_frequencies.output.frequencies,
        model = lambda wildcards: config["builds"][wildcards.type][wildcards.sample]["best_predictor"]
    output:
        node_data = BUILD_TIMEPOINT_PATH + "forecasts.json",
        frequencies = "results/auspice/flu_" + BUILD_SEGMENT_LOG_STEM + "_tip-frequencies.json"
    params:
        delta_months = _get_delta_months_to_forecast
    conda: "../envs/popcast.yaml"
    shell:
        """
        popcast forecast \
            --tip-attributes {input.attributes} \
            --frequencies {input.frequencies} \
            --model {input.model} \
            --delta-months {params.delta_months} \
            --output-node-data {output.node_data} \
            --output-frequencies {output.frequencies}
        """


rule forecast_all_tips:
    input:
        attributes = rules.annotate_weighted_distances_for_tip_attributes.output.attributes,
        model = _get_model_from_validation
    output:
        table = BUILD_PATH + "forecasts_{predictors}.tsv",
    params:
        delta_months = config["fitness_model"]["delta_months_to_fit"]
    conda: "../envs/popcast.yaml"
    shell:
        """
        popcast forecast \
            --tip-attributes {input.attributes} \
            --model {input.model} \
            --delta-months {params.delta_months} \
            --output-table {output.table}
        """


def _get_target_tip_attributes_by_wildcards(wildcards):
    validation_build = _get_validation_sample_by_wildcards(wildcards)
    wildcards_dict = dict(wildcards)
    wildcards_dict["sample"] = validation_build
    return (BUILD_PATH + "tip_attributes_with_naive_predictor.tsv").format(
        **wildcards_dict,
    )


rule test_distance_models:
    input:
        attributes = rules.annotate_naive_tip_attribute.output.attributes,
        target_attributes = _get_target_tip_attributes_by_wildcards,
        model = _get_model_to_test_by_wildcards,
    output:
        model = BUILD_PATH + "test_models_by_distances/{delta_month}/{predictors}.json",
        errors = BUILD_PATH + "test_models_by_distances_errors/{delta_month}/{predictors}.tsv",
        coefficients = BUILD_PATH + "test_models_by_distances_coefficients/{delta_month}/{predictors}.tsv"
    conda: "../envs/popcast.yaml"
    benchmark: "benchmarks/test_fitness_model_distances_" + BUILD_LOG_STEM + "_{delta_month}_{predictors}.txt"
    log: "logs/test_fitness_model_distances_" + BUILD_LOG_STEM + "_{delta_month}_{predictors}.txt"
    resources:
        mem_mb=20000,
        runtime="0:20:00",
        partition="campus-new",
        qos="campus-new",
    shell:
        """
        popcast fit \
            --tip-attributes {input.attributes} \
            --target-tip-attributes {input.target_attributes} \
            --target distances \
            --fixed-model {input.model} \
            --delta-months {wildcards.delta_month} \
            --prefer-user-delta-months \
            --errors-by-timepoint {output.errors} \
            --coefficients-by-timepoint {output.coefficients} \
            --include-scores \
            --output {output.model} &> {log}
        """


rule annotate_test_distance_models:
    input:
        attributes = rules.annotate_naive_tip_attribute.output.attributes,
        model = rules.test_distance_models.output.model,
        errors = rules.test_distance_models.output.errors,
        coefficients = rules.test_distance_models.output.coefficients
    output:
        errors = BUILD_PATH + "annotated_test_models_by_distances_errors_by_horizon/{delta_month}/{predictors}.tsv",
        coefficients = BUILD_PATH + "annotated_test_models_by_distances_coefficients_by_horizon/{delta_month}/{predictors}.tsv",
        frequencies = BUILD_PATH + "annotated_test_models_by_distances_frequencies_by_horizon/{delta_month}/{predictors}.tsv",
    conda: "../envs/anaconda.python3.yaml"
    params:
        delay_type=lambda wildcards: config["builds"].get(wildcards.type, {}).get(wildcards.sample, {}).get("delay_type"),
    shell:
        """
        python3 workflow/scripts/annotate_model_tables.py \
            --tip-attributes {input.attributes} \
            --model {input.model} \
            --errors-by-timepoint {input.errors} \
            --coefficients-by-timepoint {input.coefficients} \
            --annotated-errors-by-timepoint {output.errors} \
            --annotated-coefficients-by-timepoint {output.coefficients} \
            --annotated-frequencies-by-timepoint {output.frequencies} \
            --delta-months {wildcards.delta_month} \
            --annotations type="{wildcards.type}" sample="{wildcards.sample}" delay_type="{params.delay_type}" delta_month="{wildcards.delta_month}"
        """


rule aggregate_annotated_test_distance_models_errors:
    input:
        errors = expand(BUILD_PATH.replace("{", "{{").replace("}", "}}") + "annotated_test_models_by_distances_errors_by_horizon/{delta_month}/{{predictors}}.tsv", delta_month=config["fitness_model"]["delta_months"]),
    output:
        errors = BUILD_PATH + "annotated_test_models_by_distances_errors/{predictors}.tsv",
    conda: "../envs/anaconda.python3.yaml"
    shell:
        """
        python3 workflow/scripts/concatenate_tables.py \
            --tables {input.errors} \
            --output {output.errors}
        """


rule aggregate_annotated_test_distance_models_frequencies:
    input:
        frequencies = expand(BUILD_PATH.replace("{", "{{").replace("}", "}}") + "annotated_test_models_by_distances_frequencies_by_horizon/{delta_month}/{{predictors}}.tsv", delta_month=config["fitness_model"]["delta_months"]),
    output:
        frequencies = BUILD_PATH + "annotated_test_models_by_distances_frequencies/{predictors}.tsv",
    conda: "../envs/anaconda.python3.yaml"
    shell:
        """
        python3 workflow/scripts/concatenate_tables.py \
            --tables {input.frequencies} \
            --output {output.frequencies}
        """


rule export_for_clade_assignment:
    input:
        tree = rules.refine.output.tree,
        metadata = _get_metadata_by_wildcards,
        auspice_config = "config/auspice_config.json",
        node_data = [
            rules.refine.output.node_data,
            rules.ancestral.output.node_data,
            rules.translate.output.node_data,
        ],
    output:
        auspice_tree = BUILD_TIMEPOINT_PATH + "auspice.json",
    conda: "../envs/anaconda.python3.yaml"
    shell:
        """
        augur export v2 \
            --tree {input.tree} \
            --metadata {input.metadata} \
            --node-data {input.node_data} \
            --auspice-config {input.auspice_config} \
            --output {output.auspice_tree} \
            --minify-json
        """


rule assign_clades:
    input:
        auspice_tree=BUILD_TIMEPOINT_PATH + "auspice.json",
        weights="config/weights_per_site_for_clades.json",
    output:
        clades=BUILD_TIMEPOINT_PATH + "clades.json",
    conda: "../envs/anaconda.python3.yaml"
    params:
        lineage="h3n2",
        segment="ha",
        clade_attribute="clade_membership",
    shell:
        """
        python3 workflow/scripts/add_new_clades.py \
            --input {input.auspice_tree} \
            --weights {input.weights} \
            --lineage {params.lineage} \
            --segment {params.segment} \
            --new-key {params.clade_attribute} \
            --output {output.clades}
        """


rule extract_clade_per_tip:
    input:
        clades_tree=BUILD_TIMEPOINT_PATH + "clades.json",
    output:
        tips_to_clades=BUILD_TIMEPOINT_PATH + "tips_to_clades.tsv",
    conda: "../envs/anaconda.python3.yaml"
    params:
        clade_attribute="clade_membership",
    shell:
        """
        python3 workflow/scripts/auspice_tree_to_table.py \
            --tree {input.clades_tree} \
            --attributes {params.clade_attribute:q} \
            --output-metadata {output.tips_to_clades}
        """


rule collect_annotated_tip_clade_tables:
    input:
        table=_get_tip_clades_by_wildcards,
    output:
        tip_clade_table = BUILD_PATH + "tips_to_clades.tsv",
    conda: "../envs/csv.yaml"
    shell:
        """
        csvtk rename -t -f name -n strain {input.table} > {output.tip_clade_table}
        """


def _get_tips_to_clades_for_full_tree_by_wildcards(wildcards):
    full_tree_sample = _get_full_tree_sample_by_wildcards(wildcards)
    tip_clades_path = BUILD_PATH + "tips_to_clades.tsv"
    return tip_clades_path.format(type=wildcards.type, sample=full_tree_sample)


rule merge_metadata_and_clades:
    input:
        metadata=_get_metadata_by_wildcards,
        clades=_get_tips_to_clades_for_full_tree_by_wildcards,
    output:
        metadata=BUILD_PATH + "metadata_with_clades.tsv",
    conda: "../envs/csv.yaml"
    shell:
        """
        csvtk join -t \
            --fields strain \
            {input.metadata} \
            {input.clades} > {output}
        """


rule export:
    input:
        tree = rules.refine.output.tree,
        metadata = BUILD_PATH + "metadata_with_clades.tsv",
        auspice_config = "config/auspice_config.json",
        node_data = _get_node_data_for_export,
        forecasts = rules.forecast_tips.output.node_data,
        colors = "config/colors.tsv"
    output:
        auspice_tree = "results/auspice/flu_" + BUILD_SEGMENT_LOG_STEM + ".json",
    conda: "../envs/anaconda.python3.yaml"
    shell:
        """
        augur export v2 \
            --tree {input.tree} \
            --metadata {input.metadata} \
            --node-data {input.node_data} {input.forecasts} \
            --colors {input.colors} \
            --auspice-config {input.auspice_config} \
            --output {output.auspice_tree} \
            --minify-json
        """


rule plot_validation_figure:
    input:
        attributes = rules.annotate_weighted_distances_for_tip_attributes.output.attributes,
        tips_to_clades = _get_tips_to_clades_for_full_tree_by_wildcards,
        forecasts = rules.forecast_all_tips.output.table,
        model_errors = _get_model_validation_errors
    output:
        figure = "manuscript/figures/validation_figure_{type}-{sample}-{predictors}.pdf",
        clades = "manuscript/figures/validation_figure_clades_{type}-{sample}-{predictors}.tsv",
        ranks = "manuscript/figures/validation_figure_ranks_{type}-{sample}-{predictors}.tsv"
    conda: "../envs/anaconda.python3.yaml"
    shell:
        """
        python3 workflow/scripts/plot_validation_figure_by_population.py \
            --tip-attributes {input.attributes} \
            --tips-to-clades {input.tips_to_clades} \
            --forecasts {input.forecasts} \
            --model-errors {input.model_errors} \
            --population {wildcards.type} \
            --sample {wildcards.sample} \
            --predictors {wildcards.predictors} \
            --output {output.figure} \
            --output-clades-table {output.clades} \
            --output-ranks-table {output.ranks}
        """
