import marimo

__generated_with = "0.14.17"
app = marimo.App(width="medium")


@app.cell
def _():
    import marimo as mo

    import argparse

    import matplotlib.pyplot as plt
    import numpy as np
    import pandas as pd
    import seaborn as sns
    return argparse, pd, plt, sns


@app.cell
def _(argparse):
    parser = argparse.ArgumentParser()

    parser.add_argument("--effects", help="path to CSV with effects of realistic interventions on optimal distances to the future")
    parser.add_argument("--output", help="plot showing intervention effects by future clade entropy")

    args = parser.parse_args()
    return (args,)


@app.cell
def _(args):
    effects_of_interventions_path = args.effects if args.effects else "manuscript/tables/h3n2_optimal_effects_of_realistic_interventions_on_distances_to_the_future.csv"
    return (effects_of_interventions_path,)


@app.cell
def _(args):
    output_path = args.output if args.output else "manuscript/figures/h3n2_optimal_effects_of_realistic_interventions_on_distances_to_the_future_by_hemisphere.pdf"
    return (output_path,)


@app.cell
def _(effects_of_interventions_path, pd):
    effects_of_interventions = pd.read_csv(
        effects_of_interventions_path,
    )

    effects_of_interventions["hemisphere"] = effects_of_interventions["future_timepoint"].apply(
        lambda date: "northern" if date.split("-")[1] in ["10", "01"] else "southern"
    )
    return (effects_of_interventions,)


@app.cell
def _(effects_of_interventions):
    effects_of_interventions
    return


@app.cell
def _():
    intervention_order = [
        "improved vaccine",
        "improved surveillance",
        "improved vaccine and surveillance",
    ]
    return (intervention_order,)


@app.cell
def _(effects_of_interventions, intervention_order, output_path, plt, sns):
    fig, ax = plt.subplots(1, 1, figsize=(10, 5), dpi=300)

    sns.violinplot(
        x="intervention_name",
        y="difference_in_optimal_distance",
        hue="hemisphere",
        hue_order=["northern", "southern"],
        data=effects_of_interventions,
        order=intervention_order,
        #color="#FFFFFF",
        fill=False,
        cut=0,
        inner="quartile",
        ax=ax,
    )
    sns.stripplot(
        x="intervention_name",
        y="difference_in_optimal_distance",
        hue="hemisphere",
        hue_order=["northern", "southern"],
        data=effects_of_interventions,
        order=intervention_order,
        alpha=0.35,
        ax=ax,
        dodge=True,
        legend=False,
    )

    ax.axhline(
        y=0,
        color="#000000",
        zorder=-10,
        linewidth=1,
    )

    ax.set_xlabel("Intervention")
    ax.set_ylabel("Difference in optimal distance to future per timepoint\n(status quo - intervention)")

    ax.set_ylim(
        bottom=effects_of_interventions["difference_in_optimal_distance"].min() - 0.75,
    )

    ax.text(
        0.5,
        0.95,
        "more accurate",
        horizontalalignment='center',
        verticalalignment='center',
        transform=ax.transAxes,
    )

    ax.text(
        0.5,
        0.05,
        "less accurate",
        horizontalalignment='center',
        verticalalignment='center',
        transform=ax.transAxes,
    )

    sns.despine()
    plt.tight_layout()
    plt.savefig(output_path)
    plt.show()
    return


@app.cell
def _():
    return


if __name__ == "__main__":
    app.run()
