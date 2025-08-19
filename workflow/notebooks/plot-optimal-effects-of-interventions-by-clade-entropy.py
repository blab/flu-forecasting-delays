import marimo

__generated_with = "0.14.17"
app = marimo.App(width="medium")


@app.cell
def _():
    import marimo as mo

    import argparse

    import altair as alt
    import pandas as pd
    from scipy.stats import entropy, linregress
    return alt, argparse, entropy, linregress, pd


@app.cell
def _(argparse):
    parser = argparse.ArgumentParser()

    parser.add_argument("--clade-frequencies", help="path to CSV of clade frequencies per timepoint, delay type, and forecast horizon")
    parser.add_argument("--effects", help="path to CSV with effects of realistic interventions on optimal distances to the future")
    parser.add_argument("--output", help="plot showing intervention effects by future clade entropy")

    args = parser.parse_args()
    return (args,)


@app.cell
def _(args):
    clade_frequencies_path = args.clade_frequencies if args.clade_frequencies else "results/clade_frequencies_for_h3n2.csv"
    return (clade_frequencies_path,)


@app.cell
def _(args):
    effects_of_interventions_path = args.effects if args.effects else "manuscript/tables/h3n2_optimal_effects_of_realistic_interventions_on_distances_to_the_future.csv"
    return (effects_of_interventions_path,)


@app.cell
def _(args):
    output_path = args.output if args.output else "manuscript/figures/h3n2_optimal_effects_of_realistic_interventions_on_distances_to_the_future_by_future_clade_entropy.pdf"
    return (output_path,)


@app.cell
def _(clade_frequencies_path, pd):
    clade_frequencies = pd.read_csv(
        clade_frequencies_path,
    )
    return (clade_frequencies,)


@app.cell
def _(clade_frequencies):
    clade_frequencies
    return


@app.cell
def _(clade_frequencies):
    clade_frequencies["timepoint"].drop_duplicates()
    return


@app.cell
def _(effects_of_interventions_path, pd):
    effects_of_interventions = pd.read_csv(
        effects_of_interventions_path,
    )
    return (effects_of_interventions,)


@app.cell
def _(effects_of_interventions):
    effects_of_interventions
    return


@app.cell
def _(effects_of_interventions):
    future_timepoints = effects_of_interventions["future_timepoint"].drop_duplicates().values
    return (future_timepoints,)


@app.cell
def _(future_timepoints):
    future_timepoints
    return


@app.cell
def _(clade_frequencies, future_timepoints):
    future_clade_frequencies = clade_frequencies.loc[
        (
            (clade_frequencies["sample"] == "h3n2_no_delay") &
            (clade_frequencies["delta_month"] == 12) &
            (clade_frequencies["future_timepoint"].isin(future_timepoints))
        ),
        [
            "clade_membership",
            "observed_frequency",
            "future_timepoint",
        ]
    ].copy()
    return (future_clade_frequencies,)


@app.cell
def _(future_clade_frequencies):
    future_clade_frequencies
    return


@app.cell
def _(entropy, future_clade_frequencies):
    entropy_per_future_timepoint = future_clade_frequencies.groupby(["future_timepoint"]).aggregate(
        entropy=("observed_frequency", entropy),
    ).reset_index()
    return (entropy_per_future_timepoint,)


@app.cell
def _(entropy_per_future_timepoint):
    entropy_per_future_timepoint
    return


@app.cell
def _(effects_of_interventions, entropy_per_future_timepoint):
    effects_of_interventions_with_entropy = effects_of_interventions.merge(
        entropy_per_future_timepoint,
        on="future_timepoint",
    )
    return (effects_of_interventions_with_entropy,)


@app.cell
def _(effects_of_interventions_with_entropy):
    effects_of_interventions_with_entropy
    return


@app.cell
def _(effects_of_interventions_with_entropy, linregress):
    r_by_intervention = {}
    for intervention, intervention_df in effects_of_interventions_with_entropy.groupby("intervention_name"):
        slope, intercept, rval, pval, stderr = linregress(
            intervention_df["entropy"],
            intervention_df["difference_in_optimal_distance"],
        )
        r_by_intervention[intervention] = rval
    return (r_by_intervention,)


@app.cell
def _(r_by_intervention):
    r_by_intervention
    return


@app.cell
def _(effects_of_interventions_with_entropy, r_by_intervention):
    effects_of_interventions_with_entropy_and_r = effects_of_interventions_with_entropy.assign(
        intervention_name_with_r=effects_of_interventions_with_entropy["intervention_name"].apply(
            lambda name: name + f" (r={r_by_intervention[name]:.2f})"
        )
    )
    return (effects_of_interventions_with_entropy_and_r,)


@app.cell
def _(alt, effects_of_interventions_with_entropy_and_r, output_path):
    chart = alt.Chart(effects_of_interventions_with_entropy_and_r).mark_point().encode(
        x=alt.X("entropy:Q", title="future clade entropy").scale(zero=False),
        y=alt.Y("difference_in_optimal_distance:Q", title="difference in optimal distance to the future"),
    ).facet(
        alt.Column("intervention_name_with_r:N", title=None),
    )

    chart.save(output_path, ppi=300)
    chart
    return


@app.cell
def _():
    return


if __name__ == "__main__":
    app.run()
