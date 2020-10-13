import argparse
import numpy as np
import pandas as pd
from scipy.stats import gamma


if __name__ == "__main__":
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("--metadata", required=True, help="tab-delimited metadata with a date field")
    parser.add_argument("--date-field", required=True, help="name of the field containing sample collection date")
    parser.add_argument("--submission-field", required=True, nargs="+", help="name of output field(s) for resulting random submission dates")
    parser.add_argument("--shape", required=True, type=float, nargs="+", help="shape value(s) for gamma distribution(s)")
    parser.add_argument("--location", required=True, type=float, nargs="+", help="location value(s) for gamma distribution(s)")
    parser.add_argument("--scale", required=True, type=float, nargs="+", help="scale value(s) for gamma distribution(s)")
    parser.add_argument("--output", required=True, help="metadata annotated with random submission dates")

    args = parser.parse_args()

    metadata = pd.read_csv(args.metadata, sep="\t")

    # Find indices of records with unambiguous dates.
    has_unambiguous_dates = (~metadata["date"].str.contains("X"))
    number_of_unambiguous_dates = has_unambiguous_dates.sum()
    print(number_of_unambiguous_dates)

    # Generate samples from gamma distributions matching the given parameters
    # and store them in the given field name(s).
    for submission_field, shape, location, scale in zip(args.submission_field, args.shape, args.location, args.scale):
        random_samples = gamma.rvs(
            a=shape,
            loc=location,
            scale=scale,
            size=number_of_unambiguous_dates
        )

        random_offsets = np.array([
            pd.DateOffset(days=int(days))
            for days in random_samples
        ])

        metadata.loc[
            has_unambiguous_dates,
            submission_field
        ] = pd.to_datetime(metadata.loc[has_unambiguous_dates, args.date_field]) + random_offsets

    metadata.to_csv(args.output, sep="\t", index=False, na_rep="N/A")
