import argparse
import pandas as pd


if __name__ == "__main__":
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("--metadata", required=True, help="tab-delimited metadata with a date field")
    parser.add_argument("--submission-date-field", required=True, help="name of metadata field with the observed submission date for each strain")
    parser.add_argument("--ideal-submission-date-field", required=True, help="name of metadata field with the ideal submission date for each strain")
    parser.add_argument("--output", required=True, help="metadata with ideal submission date less than or equal to the observed submission date")

    args = parser.parse_args()

    metadata = pd.read_csv(args.metadata, sep="\t", parse_dates=[args.submission_date_field, args.ideal_submission_date_field])

    # Find fields where the randomly assigned ideal submission date is greater than the observed submission date.
    ideal_later_than_observed = (metadata[args.ideal_submission_date_field] > metadata[args.submission_date_field])

    # Replace the ideal submission date with the observed when the ideal is later than the observed.
    metadata.loc[ideal_later_than_observed, args.ideal_submission_date_field] = metadata.loc[ideal_later_than_observed, args.submission_date_field]

    metadata.to_csv(args.output, sep="\t", index=False, na_rep="N/A")
