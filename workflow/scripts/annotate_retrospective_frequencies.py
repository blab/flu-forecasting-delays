"""Annotate retrospective frequencies to a data frame of tip attribute
"""
import argparse
import csv
import json
import numpy as np

from augur.dates import numeric_date


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description="Annotate retrospective frequencies to a data frame of tip attributes",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument("--tip-attributes", required=True, help="TSV table of tip attributes with at least 'strain' and 'timepoint' columns")
    parser.add_argument("--frequencies", required=True, help="tip frequencies JSON with all strains in the given tip attributes table represented")
    parser.add_argument("--output", required=True, help="TSV table of tip attributes with 'retrospective_frequency' column added for each strain at each timepoint")
    args = parser.parse_args()

    # Load frequencies.
    with open(args.frequencies, "r") as fh:
        frequencies_json = json.load(fh)

    frequencies_by_strain = {
        strain: strain_frequencies["frequencies"]
        for strain, strain_frequencies in frequencies_json.items()
        if "frequencies" in strain_frequencies
    }

    pivots = list(np.array(frequencies_json.pop("pivots")).round(2))
    print(pivots)
    pivot_by_timepoint = {}

    with open(args.tip_attributes, "r", encoding="utf-8", newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        fieldnames = reader.fieldnames + ["retrospective_frequency"]

        with open(args.output, "w", encoding="utf-8", newline="") as oh:
            writer = csv.DictWriter(oh, fieldnames=fieldnames, delimiter="\t", lineterminator="\n")
            writer.writeheader()

            for row in reader:
                if row["timepoint"] not in pivot_by_timepoint:
                    pivot_by_timepoint[row["timepoint"]] = np.round(numeric_date(row["timepoint"]), 2)
                    print(f"Processing {row['timepoint']} with pivot of {pivot_by_timepoint[row['timepoint']]}")

                pivot = pivot_by_timepoint[row["timepoint"]]
                pivot_index = pivots.index(pivot)
                row["retrospective_frequency"] = frequencies_by_strain[row["strain"]][pivot_index]
                writer.writerow(row)
