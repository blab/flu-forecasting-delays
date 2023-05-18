import argparse
from augur.dates import numeric_date
import pandas as pd


if __name__ == "__main__":
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("--input", required=True, help="distances to the future to filter and prepare")
    parser.add_argument("--samples", required=True, nargs="+", help="samples to include in the output")
    parser.add_argument("--output", required=True, help="distances to the future prepared for modeling by population")

    args = parser.parse_args()

    df = pd.read_csv(
        args.input,
        sep="\t",
        usecols=["delta_month", "sample", "validation_error", "validation_timepoint", "future_timepoint"],
    )

    month_delay_by_sample = {
        "simulated_no_delay": 0.0,
        "simulated_ideal_delay": 1.0,
        "simulated_realistic_delay": 3.0,
        "simulated_no_delay_with_bias": 0.0,
        "simulated_ideal_delay_with_bias": 1.0,
        "simulated_realistic_delay_with_bias": 3.0,
        "h3n2_no_delay": 0.0,
        "h3n2_ideal_delay": 1.0,
        "h3n2_observed_delay": 3.0,
    }

    df["delay"] = df["sample"].map(month_delay_by_sample)

    df = df[df["sample"].isin(args.samples)].copy()
    df = df.rename(
        columns={
            "validation_error": "distance",
            "delta_month": "horizon",
            "validation_timepoint": "initial_timepoint",
        }
    ).drop(columns=["sample"])

    index_by_timepoint = {timepoint: i + 1 for i, timepoint in enumerate(df["future_timepoint"].drop_duplicates().values)}
    df["t"] = df["future_timepoint"].map(index_by_timepoint)

    df["future_timepoint_numeric"] = df["future_timepoint"].apply(lambda date: round(numeric_date(date), 2))

    df = df.sort_values(["horizon", "delay", "t"])

    df.to_csv(args.output, sep="\t", index=False)
