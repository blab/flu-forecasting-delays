import argparse
import numpy as np
import pandas as pd
from scipy.stats import gamma


if __name__ == "__main__":
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("--metadata", required=True, help="tab-delimited metadata with a date field")
    parser.add_argument("--date-field", required=True, help="name of the field containing sample collection date")
    parser.add_argument("--observed-submission-date-field", required=True, help="name of the field containing the observed submission date")
    parser.add_argument("--submission-field", required=True, help="name of output field for resulting random submission dates")
    parser.add_argument("--group-by", required=True, nargs="+", help="group metadata records by these fields and assign the same randomly selected submission date to all records from the same group")
    parser.add_argument("--shape", required=True, type=float, help="shape value for gamma distribution")
    parser.add_argument("--location", required=True, type=float, help="location value for gamma distribution")
    parser.add_argument("--scale", required=True, type=float, help="scale value for gamma distribution")
    parser.add_argument("--random-seed", type=int, default=0, help="random seed for sampling from gamma distribution")
    parser.add_argument("--output", required=True, help="metadata annotated with random submission dates")

    args = parser.parse_args()

    # Set the random seed.
    np.random.seed(args.random_seed)

    metadata = pd.read_csv(
        args.metadata,
        sep="\t",
        parse_dates=[
            args.date_field,
            args.observed_submission_date_field,
        ]
    )

    # Find all distinct groups for the requested fields, so we can assign the
    # same delayed submission date to all records from the same group.
    groups = metadata.loc[:, args.group_by].drop_duplicates()

    # Generate samples from gamma distributions matching the given parameters
    # and store them in the given field name(s).
    random_samples = gamma.rvs(
        a=args.shape,
        loc=args.location,
        scale=args.scale,
        size=groups.shape[0]
    )

    random_offsets = np.array([
        pd.DateOffset(days=int(days))
        for days in random_samples
    ])

    # Annotate random offsets to groups.
    groups["_offset"] = random_offsets

    # Annotate random offsets to original records through groups.
    metadata = metadata.merge(
        groups,
        on=args.group_by,
    )

    # Calculate submission date per record based on actual collection date and
    # random offset per group.
    metadata[args.submission_field] = metadata[args.date_field] + metadata["_offset"]

    # Find fields where the randomly assigned ideal submission date is greater
    # than the observed submission date.
    ideal_later_than_observed = (metadata[args.submission_field] > metadata[args.observed_submission_date_field])
    print(ideal_later_than_observed.sum())

    # Replace the ideal submission date with the observed when the ideal is
    # later than the observed.
    metadata.loc[ideal_later_than_observed, args.submission_field] = metadata.loc[
        ideal_later_than_observed,
        args.observed_submission_date_field
    ]

    metadata.drop(columns=["_offset"]).to_csv(args.output, sep="\t", index=False, na_rep="N/A")
