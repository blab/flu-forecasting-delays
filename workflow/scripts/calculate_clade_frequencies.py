import argparse
import ipdb
import pandas as pd


if __name__ == "__main__":
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("--tips", required=True, help="tip attributes and forecasts merged by future timepoint across all samples")
    parser.add_argument("--tip-clades", required=True, help="clade to tip mappings from a single 'full tree' build without delays, such that clade labels are the same for tips across all timepoints")
    parser.add_argument("--min-clade-frequency", type=float, default=0.1, help="minimum clade frequency for an initial timepoint clade for that clade to be available to link to the future timepoint")
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
    ).sort_values([
        "strain",
        "depth",
    ])

    # Find all combinations of t - h, t, and h.
    initial_and_future_timepoints = tips.loc[
        :,
        [
            "timepoint",
            "future_timepoint",
            "delta_month",
        ]
    ].dropna().drop_duplicates().values

    # Find all delay values s.
    delays_and_samples = tips.loc[
        :,
        [
            "delay_type",
            "sample",
        ]
    ].dropna().drop_duplicates().values

    # For each combination of t, h, and t - h, find clades circulating at both t and t - h.
    all_clade_frequencies = []
    for (initial_timepoint, future_timepoint, delta_month) in initial_and_future_timepoints:
        print(f"Initial time: {initial_timepoint}, future time: {future_timepoint}")
        # Find all tips at time t without delay.
        future_tips = tips.loc[
            (
                (tips["delay_type"] == "none") &
                (tips["timepoint"] == future_timepoint)
            ),
            ("strain", "frequency")
        ].drop_duplicates()
        future_tip_names = future_tips["strain"].values

        # Find all tips at time t - h without delay.
        initial_tips = tips.loc[
            (
                (tips["delay_type"] == "none") &
                (tips["timepoint"] == initial_timepoint)
            ),
            ("strain", "frequency")
        ].drop_duplicates()
        initial_tip_names = initial_tips["strain"].values

        # Find all clades for tips at time t - h.
        initial_tip_clades = tip_clades.loc[
            tip_clades["strain"].isin(initial_tip_names)
        ].drop_duplicates()

        # Annotate tips from time t - h with clades.
        initial_tips_with_clades = initial_tip_clades.merge(
            initial_tips,
            on="strain",
        )

        # Find clades associated with each tip. This mapping of tips to clades
        # with tip frequencies will be used to link timepoints t - h and t.
        initial_tips_with_derived_clades = initial_tips_with_clades.query("depth == 0").copy()

        # Filter the global tip-to-clade mapping to just derived clades so we
        # can only map tips to extant clades.
        derived_clades = set(initial_tips_with_derived_clades["clade_membership"].values)
        initial_tips_with_clades_in_derived_clades = initial_tips_with_clades[
            initial_tips_with_clades["clade_membership"].isin(derived_clades)
        ].copy()

        # Find clade frequencies at the initial timepoint based on tips alone.
        initial_derived_clade_frequencies = initial_tips_with_derived_clades.groupby(
            "clade_membership"
        )["frequency"].sum().reset_index()

        # Find all clades with a frequency less than the minimum threshold.
        low_frequency_derived_clades = initial_derived_clade_frequencies[
            initial_derived_clade_frequencies["frequency"] < args.min_clade_frequency
        ].sort_values("frequency").reset_index(drop=True)

        # Find tips from low-frequency clades.
        low_frequency_tips = set(
            initial_tips_with_derived_clades.loc[
                initial_tips_with_derived_clades["clade_membership"].isin(
                    low_frequency_derived_clades["clade_membership"].values
                ),
                "strain"
            ]
        )

        # While there are tips from low frequency clades to process, try to
        # collapse each tip to an ancestral clade in the set of extant
        # clades. If no ancestral clade exists for a given tip, remove it from
        # the set to process. If an ancestral clade exists, collapse the tip
        # into that clade, recalculate clade frequencies, remove tips now in
        # high frequency clades, and continue to try collapsing the remaining
        # tips.
        print(f"Found {len(low_frequency_tips)} low frequency tips")
        while low_frequency_tips:
            # Select a tip from a low frequency clade.
            tip = low_frequency_tips.pop()
            print(f"Processing tip {tip}")

            # Find the parent clade for the tip from the list of extant clades.
            # The parent clade will be the first clade to have a depth greater
            # than the current clade's depth.
            current_clade = initial_tips_with_derived_clades.loc[
                initial_tips_with_derived_clades["strain"] == tip,
                "clade_membership"
            ].values[0]
            current_depth = initial_tips_with_derived_clades.loc[
                initial_tips_with_derived_clades["strain"] == tip,
                "depth"
            ].values[0]
            parent_clades = initial_tips_with_clades_in_derived_clades[
                (
                    (initial_tips_with_clades_in_derived_clades["strain"] == tip) &
                    (initial_tips_with_clades_in_derived_clades["depth"] > current_depth)
                )
            ]

            # If no parent clade exists, continue to the next tip.
            if len(parent_clades) == 0:
                print(f"No parent clade for tip {tip}, removing from low frequency tips as clade {current_clade}")
                continue

            # If a parent clade exists, assign the tip to that clade in the
            # tip-to-clade mapping.
            parent_clade = parent_clades.head(1)["clade_membership"].values[0]
            parent_depth = parent_clades.head(1)["depth"].values[0]
            initial_tips_with_derived_clades.loc[
                initial_tips_with_derived_clades["strain"] == tip,
                "clade_membership"
            ] = parent_clade

            initial_tips_with_derived_clades.loc[
                initial_tips_with_derived_clades["strain"] == tip,
                "depth"
            ] = parent_depth

            # Recalculate clade frequencies from the tip-to-clade mapping.
            initial_derived_clade_frequencies = initial_tips_with_derived_clades.groupby(
                "clade_membership"
            )["frequency"].sum().reset_index()

            # If the parent clade of the tip is above the minimum threshold,
            # remove all tips for that clade from the low frequency clades
            # list. Otherwise, add the current tip back to the list of low
            # frequency tips to process.
            if (initial_derived_clade_frequencies.loc[
                    initial_derived_clade_frequencies["clade_membership"] == parent_clade,
                    "frequency"
                ].values[0] >= args.min_clade_frequency):
                print(f"Removing tip {tip} from low frequency tips after collapsing it from {current_clade} into clade {parent_clade}")
                parent_clade_tips = set(
                    initial_tips_with_derived_clades.loc[
                        initial_tips_with_derived_clades["clade_membership"] == parent_clade,
                        "strain"
                    ].values
                )
                low_frequency_tips.difference_update(parent_clade_tips)
            else:
                print(f"Tip {tip} remains in low frequency tips after collapsing from {current_clade} to {parent_clade}")
                low_frequency_tips.add(tip)

            print(f"{len(low_frequency_tips)} low frequency tips remain")

        # Find the names of clades at the initial timepoint after collapsing
        # tips into ancestral clades.
        initial_tip_clade_names = initial_tips_with_derived_clades[
            "clade_membership"
        ].drop_duplicates().values

        # Find all clades for tips at time t that are also present at time t - h.
        all_future_tip_clades = tip_clades[
            (tip_clades["strain"].isin(future_tip_names)) &
            (tip_clades["clade_membership"].isin(initial_tip_clade_names))
        ]

        # Assign each tip at time t to its most derived clade from time t - h.
        future_tip_clades = all_future_tip_clades.sort_values([
            "strain",
            "depth",
        ]).groupby(
            "strain",
            sort=False,
        ).first()

        assert future_tip_clades.shape[0] == future_tips.shape[0]

        future_tips_with_clades = future_tips.merge(
            future_tip_clades,
            on="strain",
        )

        # Create the set of all clades present at both times t and t - h without
        # delay, C.
        future_tip_clade_names = future_tips_with_clades["clade_membership"].drop_duplicates().values

        # Sum tip frequencies by clade at time t to get observed future clade
        # frequencies.
        observed_future_clade_frequencies = future_tips_with_clades.groupby(
            "clade_membership"
        )["frequency"].sum().reset_index().rename(
            columns={"frequency": "observed_frequency"}
        )

        if initial_timepoint == "2011-04-01" and future_timepoint == "2011-07-01":
            ipdb.set_trace()

        # For each delay value s, calculate predicted clade frequencies for t
        # from t - h under s.
        for (delay, sample) in delays_and_samples:
            # Find all tips at time t - h under delay of s.
            initial_delayed_tips = tips.loc[
                (
                    (tips["delay_type"] == delay) &
                    (tips["timepoint"] == initial_timepoint) &
                    (tips["future_timepoint"] == future_timepoint)
                ),
                (
                    "strain",
                    "frequency",
                    "projected_frequency",
                )
            ]
            initial_delayed_tip_names = initial_delayed_tips["strain"].drop_duplicates().values

            # Assign each tip at time t - h to its most derived clade from the
            # set of C.
            all_initial_delayed_tip_clades = tip_clades[
                (tip_clades["strain"].isin(initial_delayed_tip_names)) &
                (tip_clades["clade_membership"].isin(future_tip_clade_names))
            ]

            initial_tip_clades = all_initial_delayed_tip_clades.sort_values([
                "strain",
                "depth",
            ]).groupby(
                "strain",
                sort=False,
            ).first()

            initial_delayed_tips_with_clades = initial_delayed_tips.merge(
                initial_tip_clades,
                on="strain",
            )

            # Sum predicted tip frequencies by clade at time t - h to get
            # predicted future clade frequencies.
            initial_delayed_clade_frequencies = initial_delayed_tips_with_clades.groupby(
                "clade_membership"
            ).aggregate({
                "frequency": "sum",
                "projected_frequency": "sum",
            }).reset_index()

            # Calculate forecast error from the difference between observed and
            # predicted future clade frequencies. It is possible that a clade
            # that existed at both times t and t - h without delay does not
            # exist at time t - h under a delay, so we left join from the
            # observed clade frequencies to the initial delay clade frequencies
            # and fill missing clade frequencies with zeros.
            clade_frequencies = observed_future_clade_frequencies.merge(
                initial_delayed_clade_frequencies,
                on="clade_membership",
                how="left",
            ).fillna(0.0)
            clade_frequencies["forecast_error"] = clade_frequencies["observed_frequency"] - clade_frequencies["projected_frequency"]
            clade_frequencies["absolute_forecast_error"] = clade_frequencies["forecast_error"].abs()

            # Record the values of t, t - h, h, s, c, observed future clade
            # frequency, predicted future clade frequency, forecast error, and
            # absolute forecast error. Reannotate attributes for this set of
            # records.
            clade_frequencies["sample"] = sample
            clade_frequencies["timepoint"] = initial_timepoint
            clade_frequencies["future_timepoint"] = future_timepoint
            clade_frequencies["delta_month"] = delta_month
            clade_frequencies["delay_type"] = delay

            all_clade_frequencies.append(clade_frequencies)

    # Collect clade frequencies across all samples, timepoints, etc.
    all_clade_frequencies_df = pd.concat(all_clade_frequencies, ignore_index=True)
    all_clade_frequencies_df.to_csv(
        args.output,
        index=False,
        sep="\t",
    )
