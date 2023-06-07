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
    parser.add_argument("--random-seed", type=int, default=0, help="random seed for sampling from gamma distribution(s)")
    parser.add_argument(
        "--bias-delay-by-fitness",
        action="store_true",
        help="""bias submission delays toward higher values for samples with the maximum fitness in their generation.
        Assumes that the metadata are from simulated populations with both 'generation' and 'fitness' columns."""
    )
    parser.add_argument("--output", required=True, help="metadata annotated with random submission dates")

    args = parser.parse_args()

    # Set the random seed.
    np.random.seed(args.random_seed)

    metadata = pd.read_csv(args.metadata, sep="\t", parse_dates=[args.date_field])

    if args.bias_delay_by_fitness:
        # Normalize fitness values per generation to simplify the logic downstream.
        max_fitness_by_generation = metadata.groupby("generation")["fitness"].max().reset_index()
        metadata = metadata.merge(max_fitness_by_generation, on="generation", how="left", suffixes=["", "_max"])
        metadata["normalized_fitness"] = metadata["fitness"] / metadata["fitness_max"]

    # Generate samples from gamma distributions matching the given parameters
    # and store them in the given field name(s).
    for submission_field, shape, location, scale in zip(args.submission_field, args.shape, args.location, args.scale):
        random_samples = gamma.rvs(
            a=shape,
            loc=location,
            scale=scale,
            size=metadata.shape[0]
        )

        if args.bias_delay_by_fitness:
            # When biasing the delay by normalized fitness, create an alternate
            # sample of values and select the maximum value between the original
            # and alternate when a sample has fitness == 1. This approach
            # increases the delay for the highest fitness samples on average
            # while keeping the average delay for other samples nearly the same.
            alternate_random_samples = gamma.rvs(
                a=shape,
                loc=location,
                scale=scale,
                size=metadata.shape[0]
            )
            random_samples = [
                random_samples[i] if fitness < 1 else max(random_samples[i], alternate_random_samples[i])
                for i, fitness in enumerate(metadata["normalized_fitness"].values)
            ]

        random_offsets = np.array([
            pd.DateOffset(days=int(days))
            for days in random_samples
        ])

        metadata[submission_field] = metadata[args.date_field] + random_offsets

    # Find fields where the randomly assigned ideal submission date is greater
    # than the realistic submission date.
    ideal_later_than_realistic = (metadata["ideal_submission_date"] > metadata["realistic_submission_date"])
    print(ideal_later_than_realistic.sum())

    # Swap the ideal and realistic submission dates when the ideal is
    # later than the realistic.
    early_realistic = metadata.loc[ideal_later_than_realistic, "realistic_submission_date"].values
    late_ideal = metadata.loc[ideal_later_than_realistic, "ideal_submission_date"].values

    metadata.loc[ideal_later_than_realistic, "ideal_submission_date"] = early_realistic
    metadata.loc[ideal_later_than_realistic, "realistic_submission_date"] = late_ideal

    metadata.to_csv(args.output, sep="\t", index=False, na_rep="N/A")
