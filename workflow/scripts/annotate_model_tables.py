import argparse
import json
import pandas as pd


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--tip-attributes", required=True, help="tab-delimited file describing tip attributes at all timepoints with standardized predictors and weighted distances to the future")
    parser.add_argument("--model", required=True, help="JSON representing the model fit with training and cross-validation results, beta coefficients for predictors, and summary statistics")
    parser.add_argument("--errors-by-timepoint", help="data frame of cross-validation errors by validation timepoint")
    parser.add_argument("--coefficients-by-timepoint", help="data frame of coefficients by validation timepoint")
    parser.add_argument("--annotated-errors-by-timepoint", help="annotated model errors by timepoint")
    parser.add_argument("--annotated-coefficients-by-timepoint", help="annotated model coefficients by timepoint")
    parser.add_argument("--delta-months", type=int, help="number of months to project clade frequencies into the future")
    parser.add_argument("--annotations", nargs="+", help="additional annotations to add to the output table in the format of 'key=value' pairs")

    args = parser.parse_args()

    # Load tip attributes to calculate within-timepoint diversity.
    tips = pd.read_csv(args.tip_attributes, sep="\t", parse_dates=["timepoint"])

    # Load the model JSON to get access to projected frequencies for tips.
    with open(args.model, "r") as fh:
        model = json.load(fh)

    # Collect all projected frequencies and weighted distances, to enable
    # calculation of weighted average distances within and between seasons.
    df = pd.concat([
        pd.DataFrame(scores["validation_data"]["y_hat"])
        for scores in model["scores"]
    ])

    # Load the original model table output for validation/test errors.
    errors = pd.read_csv(args.errors_by_timepoint, sep="\t", parse_dates=["validation_timepoint"])
    errors["future_timepoint"] = errors["validation_timepoint"] + pd.DateOffset(months=args.delta_months)

    # Load coefficients to which annotations will be added.
    coefficients = pd.read_csv(args.coefficients_by_timepoint, sep="\t")

    # Add any additional annotations requested by the user in the format of
    # "key=value" pairs where each key becomes a new column with the given
    # value.
    if args.annotations:
        for annotation in args.annotations:
            key, value = annotation.split("=")
            errors[key] = value
            coefficients[key] = value

    # Save annotated tables.
    errors.to_csv(args.annotated_errors_by_timepoint, sep="\t", header=True, index=False)
    coefficients.to_csv(args.annotated_coefficients_by_timepoint, sep="\t", header=True, index=False)
