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

        # Select current tips by delay type, horizon, and timepoint.
        current_tips = tips[
            (tips["delay_type"] == delay_type) &
            (tips["timepoint"] == initial_timepoint) &
            (tips["future_timepoint"] == future_timepoint)
        ]

        # Filter tip-clade mappings to clades associated with current and future strains.
        current_strains = set(current_tips["strain"].values)
        future_strains = set(future_tips["strain"].values)
        current_clade_names = set(tip_clades.loc[tip_clades["strain"].isin(current_strains), "clade_membership"].drop_duplicates().values)
        future_clade_names = set(tip_clades.loc[tip_clades["strain"].isin(future_strains), "clade_membership"].drop_duplicates().values)

        # Keep tip-clade mappings for clades that appear for both the current and future tips.
        current_and_future_clade_names = current_clade_names & future_clade_names
        current_and_future_tip_clades = tip_clades[tip_clades["clade_membership"].isin(current_and_future_clade_names)]

        # Find all clades associated with future tips.
        future_tips_with_clades = future_tips.merge(
            current_and_future_tip_clades,
            on=["strain"],
        )

        # Find all clades associate with current tips.
        current_tips_with_clades = current_tips.merge(
            current_and_future_tip_clades,
            on=["strain"],
            suffixes=["", "_future"],
        )

        # Sort strains by clade depth and take the first value by strain and
        # frequencies to get the most derived clade label that is present in
        # both current and future timepoints.
        current_tips_with_derived_clades = current_tips_with_clades.sort_values([
            "strain",
            "depth",
        ]).groupby([
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

        # Sort future tips by the status of whether they need to use a clade
        # from the future timepoint (instead of a matching clade from the
        # current timepoint) and then by clade depth (ascending) and take first
        # record per tip. This approach prioritizes clades that exist in the
        # current timepoint, but it also allows us to track future tips without
        # any matches in the most derived clades from the current timepoint.
        future_tips_with_derived_clades = future_tips_with_clades.sort_values([
            "strain",
            "depth",
        ]).groupby([
            "strain",
            "frequency",
        ])["clade_membership"].first().reset_index(
            name="clade_membership"
        )
        assert future_tips["strain"].values.shape[0] == future_tips_with_derived_clades["strain"].values.shape[0]

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
            how="outer",
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
