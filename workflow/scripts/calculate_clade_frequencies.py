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
        initial_tip_clade_names = initial_tips_with_clades["clade_membership"].drop_duplicates().values

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
        future_tip_clade_names = set(future_tips_with_clades["clade_membership"].drop_duplicates().values)

        # Sum tip frequencies by clade at time t to get observed future clade
        # frequencies.
        observed_future_clade_frequencies = future_tips_with_clades.groupby(
            "clade_membership"
        )["frequency"].sum().reset_index().rename(
            columns={"frequency": "observed_frequency"}
        )

        # Update C to include the most derived clades present at time t - h
        # (i.e., clades from tips) that did not have a descendant tip in time t.
        # First, find all initial tips that have a clade present in C.
        initial_tip_names_with_clades_in_future = initial_tips_with_clades.loc[
            initial_tips_with_clades["clade_membership"].isin(future_tip_clade_names),
            "strain"
        ].drop_duplicates().values

        # Find remaining initial tips that don't have clades in the future.
        initial_tips_without_clades_in_future = initial_tips_with_clades[
            ~initial_tips_with_clades["strain"].isin(initial_tip_names_with_clades_in_future)
        ]

        # Find most derived clade for each initial tip that does not have a
        # clade at time t.
        initial_tip_clades_not_in_future = initial_tips_without_clades_in_future.sort_values([
            "strain",
            "depth",
        ]).groupby(
            "strain",
            sort=False,
        ).first()

        initial_tip_clade_names_not_in_future = set(
            initial_tip_clades_not_in_future["clade_membership"].drop_duplicates().values
        )

        if len(initial_tip_clade_names_not_in_future) > 0:
            # We should not have any overlap between future clade names and initial
            # clade names that weren't found in the future.
            assert len(future_tip_clade_names & initial_tip_clade_names_not_in_future) == 0

            # Update C to include the most derived clades at time t - h missing from t.
            future_tip_clade_names.update(initial_tip_clades_not_in_future)

            # Add frequencies of zero for initial clades that didn't survive to time
            # t, so we can compare forecasts from t - h to the observed outcomes for
            # those clades.
            observed_future_clade_frequencies_for_missing_clades = pd.DataFrame({
                "clade_membership": list(initial_tip_clade_names_not_in_future),
                "observed_frequency": 0.0,
            })

            observed_future_clade_frequencies = pd.concat(
                [
                    observed_future_clade_frequencies,
                    observed_future_clade_frequencies_for_missing_clades,
                ],
                ignore_index=True,
            )

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
            # and fill missing clade frequencies with zeros. The observed future
            # clades include initial clades that did not survive to time t with
            # a future frequency of zero, so the left join should include those
            # initial clades from the delay time data.
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
