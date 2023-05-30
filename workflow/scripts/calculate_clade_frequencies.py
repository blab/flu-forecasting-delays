import argparse
import pandas as pd


if __name__ == "__main__":
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("--tips", required=True, help="tip attributes and forecasts merged by future timepoint across all samples")
    parser.add_argument("--tip-clades", required=True, help="clade to tip mappings from a single 'full tree' build without delays, such that clade labels are the same for tips across all timepoints")
    parser.add_argument("--output", required=True, help="clade frequencies per sample/delay type, forecast horizon, and timepoint")

    args = parser.parse_args()

    tips = pd.read_csv(args.tips, sep="\t")
    tip_clades = pd.read_csv(
        args.tip_clades,
        sep="\t",
        usecols=[
            "strain",
            "clade_membership",
            "depth",
        ]
    )

    all_clade_frequencies = []
    for (sample, initial_timepoint, future_timepoint, delta_month, delay_type) in tips.loc[:, ["sample", "timepoint", "future_timepoint", "delta_month", "delay_type"]].dropna().drop_duplicates().values:
        # Select future tips, keeping only distinct strains
        # and their frequencies at that timepoint.
        future_tips = tips.loc[
            (
                (tips["delay_type"] == "none") &
                (tips["timepoint"] == future_timepoint)
            ),
            ("strain", "frequency")
        ].drop_duplicates()

        # Find all clades associated with future tips.
        future_tips_with_clades = future_tips.merge(
            tip_clades,
            on=["strain"],
        )

        # Select current tips by delay type, horizon, and timepoint.
        current_tips = tips[
            (tips["delay_type"] == delay_type) &
            (tips["timepoint"] == initial_timepoint) &
            (tips["future_timepoint"] == future_timepoint)
        ]

        # Calculate initial and estimated clade frequencies using clades present
        # at the current timepoint.
        current_clade_frequencies = current_tips.groupby([
            "clade_membership"
        ]).aggregate({
            "frequency": "sum",
            "projected_frequency": "sum",
        }).reset_index()

        # Filter future tips to current clades and select the most derived clade
        # label for each future tip so we get one clade per future tip.
        current_clade_names = set(current_clade_frequencies["clade_membership"].values)
        future_tips_with_derived_clades = future_tips_with_clades[
            future_tips_with_clades["clade_membership"].isin(current_clade_names)
        ].sort_values([
            "strain",
            "depth"
        ]).groupby([
            "strain",
            "frequency"
        ])["clade_membership"].first().reset_index()

        # Group filtered future tips by clade to get observed future
        # clade frequencies without delay.
        future_clade_frequencies = future_tips_with_derived_clades.groupby([
            "clade_membership"
        ]).aggregate({
            "frequency": "sum",
        }).reset_index()

        # Collect clade frequencies from current and future timepoints.
        clade_frequencies = current_clade_frequencies.merge(
            future_clade_frequencies,
            on="clade_membership",
            how="left",
            suffixes=["", "_observed"],
        ).fillna(0).rename(
            columns={"frequency_observed": "observed_frequency"},
        )

        # Reannotate attributes for this set of records.
        clade_frequencies["sample"] = sample
        clade_frequencies["timepoint"] = initial_timepoint
        clade_frequencies["future_timepoint"] = future_timepoint
        clade_frequencies["delta_month"] = delta_month
        clade_frequencies["delay_type"] = delay_type

        all_clade_frequencies.append(clade_frequencies)

    # Collect clade frequencies across all samples, timepoints, etc.
    all_clade_frequencies_df = pd.concat(all_clade_frequencies, ignore_index=True)
    all_clade_frequencies_df.to_csv(
        args.output,
        index=False,
        sep="\t",
    )
