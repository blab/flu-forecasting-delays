import argparse
import pandas as pd


if __name__ == "__main__":
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("--tips", required=True, help="tip attributes and forecasts merged by future timepoint across all samples")
    parser.add_argument("--tip-clades", required=True, help="clade to tip mappings from builds without delays")
    parser.add_argument("--output", required=True, help="clade frequencies per sample/delay type, forecast horizon, and timepoint")

    args = parser.parse_args()

    tips = pd.read_csv(args.tips, sep="\t")
    tip_clades = pd.read_csv(args.tip_clades, sep="\t")

    all_clade_frequencies = []
    for (sample, initial_timepoint, future_timepoint, delta_month, delay_type) in tips.loc[:, ["sample", "timepoint", "future_timepoint", "delta_month", "delay_type"]].dropna().drop_duplicates().values:
        # Select future tip clades.
        future_tip_clades = tip_clades[tip_clades["timepoint"] == future_timepoint]

        # Select future tips, keeping only distinct strains
        # and their frequencies at that timepoint.
        future_tips = tips.loc[
            (
                (tips["delay_type"] == "none") &
                (tips["timepoint"] == future_timepoint)
            ),
            ("timepoint", "strain", "frequency")
        ].drop_duplicates()

        # Find all clades associated with future tips.
        future_tips_with_clades = future_tips.merge(
            future_tip_clades,
            on=["timepoint", "strain"],
        )

        # Find distinct clades associated with future tips.
        distinct_future_clades = set(future_tips_with_clades["clade_membership"].drop_duplicates().values)

        # Select current tips by delay type, horizon, and timepoint.
        current_tips = tips[
            (tips["delay_type"] == delay_type) &
            (tips["timepoint"] == initial_timepoint) &
            (tips["future_timepoint"] == future_timepoint)
        ]

        # Find all future clades associate with current tips.
        current_tips_with_clades = current_tips.merge(
            future_tip_clades,
            left_on=["future_timepoint", "strain"],
            right_on=["timepoint", "strain",],
            suffixes=["", "_future"],
        )

        # Filter to future clades that are present for future tips
        # and current tips.
        current_tips_with_clades_in_future = current_tips_with_clades[current_tips_with_clades["clade_membership_future"].isin(distinct_future_clades)]

        # Sort strains by clade depth and take the first value
        # by timepoint, strain, and frequencies to get the most
        # derived clade label that is present in both current and
        # future timepoints.
        current_tips_with_derived_clades = current_tips_with_clades_in_future.sort_values([
            "strain",
            "depth",
        ]).groupby([
            "timepoint",
            "strain",
            "frequency",
            "projected_frequency"
        ])["clade_membership_future"].first().reset_index(
            name="clade_membership_future"
        )

        # Calculate initial and estimated clade frequencies.
        current_clade_frequencies = current_tips_with_derived_clades.groupby([
        "clade_membership_future"
        ]).aggregate({
            "frequency": "sum",
            "projected_frequency": "sum",
        }).reset_index().rename(
            columns={"clade_membership_future": "clade_membership"},
        )

        # Find distinct clades used for current tips.
        distinct_current_clades = set(current_clade_frequencies["clade_membership"].drop_duplicates().values)

        # Filter future tips by clade to distinct clades for current and future tips.
        future_tips_with_clades_in_current = future_tips_with_clades[future_tips_with_clades["clade_membership"].isin(distinct_current_clades)]

        # Sort future tips by clade depth (ascending) and take first
        # record per tip.
        future_tips_with_derived_clades = future_tips_with_clades_in_current.sort_values([
            "strain",
            "depth",
        ]).groupby([
            "timepoint",
            "strain",
            "frequency",
        ])["clade_membership"].first().reset_index(
            name="clade_membership"
        )

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
        ).fillna(0)

        # Reannotate attributes for this set of records.
        clade_frequencies["sample"] = sample
        clade_frequencies["timepoint"] = initial_timepoint
        clade_frequencies["future_timepoint"] = future_timepoint
        clade_frequencies["delta_month"] = delta_month
        clade_frequencies["delay_type"] = delay_type

        all_clade_frequencies.append(clade_frequencies)

    # Collect clade frequencies across all samples, timepoints, etc.
    all_clade_frequencies_df = pd.concat(all_clade_frequencies, ignore_index=True).rename(
        columns={"frequency_observed": "observed_frequency"},
    )
    all_clade_frequencies_df.to_csv(
        args.output,
        index=False,
        sep="\t",
    )
