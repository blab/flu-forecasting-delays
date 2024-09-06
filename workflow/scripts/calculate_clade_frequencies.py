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
        print(f"Initial time: {initial_timepoint}, future time: {future_timepoint}")
        # Find all tips at time t without delay.
        future_tips = tips.loc[
            (
                (tips["delay_type"] == "none") &
                (tips["timepoint"] == future_timepoint)
            ),
            ("strain", "frequency")
        ].drop_duplicates()

        # Find the clade for each tip at time t.
        future_tips = future_tips.merge(
            tip_clades,
            on="strain",
        )

        # Calculate observed clade frequencies at time t.
        observed_future_clade_frequencies = future_tips.groupby(
            "clade_membership"
        ).aggregate({
            "frequency": "sum",
            "strain": "count",
        }).reset_index()

        future_clade_names = set(observed_future_clade_frequencies["clade_membership"].values)

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

            # Find the clade for each tip at time t - h with delay of s.
            initial_delayed_tips = initial_delayed_tips.merge(
                tip_clades,
                on="strain",
            )

            # Sum initial and predicted tip frequencies by clade at time t - h
            # to get initial and predicted future clade frequencies.
            initial_delayed_clade_frequencies = initial_delayed_tips.groupby(
                "clade_membership"
            ).aggregate({
                "frequency": "sum",
                "projected_frequency": "sum",
                "strain": "count",
            }).reset_index().rename(
                columns={
                    "strain": "n_strains",
                }
            )

            # Map clades from time t to clades at t - h with delay s.
            initial_clade_names = set(initial_delayed_clade_frequencies["clade_membership"].values)
            future_to_initial_clade_map = {}
            for future_clade in future_clade_names:
                mapped_future_clade = future_clade

                # Look for the future clade name in the list of initial clade
                # names, progressively stripping the suffix from the future
                # clade until we find a match in the initial clades or we run
                # out of suffixes to remove. For example, if the future clade is
                # A.1.1.3 and there is an initial clade named A.1, we remove
                # suffixes from the future clade (A.1.1.3 -> A.1.1 -> A.1) until
                # we find the match. If the future clade is A.2 and no initial
                # clade matches, we stop after removing the only suffix and use
                # "A".
                while mapped_future_clade not in initial_clade_names and "." in mapped_future_clade:
                    mapped_future_clade = ".".join(mapped_future_clade.split(".")[:-1])

                future_to_initial_clade_map[future_clade] = mapped_future_clade
                if future_clade != mapped_future_clade:
                    print(f"Mapping clade '{future_clade}' at future timepoint {future_timepoint} to clade '{mapped_future_clade}' from initial timepoint {initial_timepoint} under delay {delay}")

            # Add initial clade labels to the observed future clade frequencies.
            observed_future_clade_frequencies_with_initial_names = observed_future_clade_frequencies.copy()
            observed_future_clade_frequencies_with_initial_names["initial_clade_membership"] = observed_future_clade_frequencies_with_initial_names["clade_membership"].map(
                future_to_initial_clade_map
            )

            # Recalculate observed future frequencies for the initial clades.
            observed_future_clade_frequencies_with_initial_names = observed_future_clade_frequencies_with_initial_names.groupby(
                "initial_clade_membership"
            ).aggregate({
                "frequency": "sum",
                "strain": "sum",
            }).reset_index().rename(
                columns={
                    "initial_clade_membership": "clade_membership",
                    "frequency": "observed_frequency",
                    "strain": "future_n_strains",
                }
            )

            # Calculate forecast error from the difference between observed and
            # predicted future clade frequencies. All future clades should now
            # be mapped to a corresponding equal or ancestral clade at the
            # initial timepoint. However, it is possible that a clade that
            # existed at time t - h does not exist at time t because it died
            # out, so we left join from the initial clade frequencies to the
            # future clade frequencies and fill missing clade frequencies with
            # zeros to account.
            clade_frequencies = initial_delayed_clade_frequencies.merge(
                observed_future_clade_frequencies_with_initial_names,
                on="clade_membership",
                how="left",
            ).fillna(0.0)

            # Recalculate nested clade frequencies for each clade based on the
            # its frequencies assigned above and those of its children. This
            # step forces the observed future frequencies for each clade to be
            # the same at the same future timepoint even when the mapping of
            # clades between timepoints above varies with the delay.
            nested_clade_frequencies = []
            for clade in clade_frequencies["clade_membership"].drop_duplicates().values:
                # The following approach relies on the fact that child clades
                # always have a prefix that matches their parent clades, so we
                # can search on the current clade name plus a "." to find the
                # children. This approach would not work if we implemented
                # aliasing of clade names.
                clade_dataframe = clade_frequencies.loc[
                    (
                        (clade_frequencies["clade_membership"] == clade) |
                        (clade_frequencies["clade_membership"].str.startswith(f"{clade}."))
                    ),
                ].copy()
                clade_dataframe["clade_membership"] = clade
                clade_dataframe = clade_dataframe.groupby("clade_membership").sum().reset_index()
                nested_clade_frequencies.append(clade_dataframe)

            clade_frequencies = pd.concat(nested_clade_frequencies, ignore_index=True)

            # Calculate forecast error for the nested clade frequencies.
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
            clade_frequencies["future_n_strains"] = clade_frequencies["future_n_strains"].astype(int)

            all_clade_frequencies.append(clade_frequencies)

    # Collect clade frequencies across all samples, timepoints, etc.
    all_clade_frequencies_df = pd.concat(all_clade_frequencies, ignore_index=True)
    all_clade_frequencies_df.to_csv(
        args.output,
        index=False,
    )
