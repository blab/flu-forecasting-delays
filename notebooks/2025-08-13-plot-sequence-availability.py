import marimo

__generated_with = "0.13.15"
app = marimo.App(width="medium")


@app.cell
def _():
    import altair as alt
    import marimo as mo
    from matplotlib import pyplot as plt
    import numpy as np
    import pandas as pd
    return alt, mo, np, pd


@app.cell
def _(mo):
    mo.md(r"""## Prepare data""")
    return


@app.cell
def _(mo):
    mo.md(r"""### All available data""")
    return


@app.cell
def _(pd):
    full_metadata = pd.read_csv(
        "data/natural/h3n2/annotated_metadata.tsv",
        sep="\t",
        parse_dates=["date"],
        usecols=["strain", "region", "date"],
    ).query(
        "(region != '?') & (date >= '2005-04-01') & (date <= '2019-10-01') & ('X' not in date)"
    )

    full_metadata["year"] = full_metadata["date"].dt.year
    full_metadata["month"] = full_metadata["date"].dt.month

    #full_metadata["year_month_date"] = full_metadata.apply(lambda row: f"{row['year']}-{row['month']:02}-01", axis=1)
    full_metadata["year_month_date"] = full_metadata.apply(lambda row: f"{row['year']}-01-01", axis=1)
    return (full_metadata,)


@app.cell
def _(full_metadata):
    full_metadata.shape[0]
    return


@app.cell
def _(full_metadata):
    full_metadata.head()
    return


@app.cell
def _(full_metadata, np):
    full_count_by_region_date = full_metadata.groupby(["region", "year_month_date"])["strain"].count().reset_index().rename(
        columns={"strain": "count"}
    )
    full_count_by_region_date["log_count"] = np.log2(full_count_by_region_date["count"])
    return (full_count_by_region_date,)


@app.cell
def _(full_count_by_region_date):
    full_count_by_region_date
    return


@app.cell
def _(mo):
    mo.md(r"""### Original subsampling""")
    return


@app.cell
def _(pd):
    subsampled_metadata = pd.read_csv(
        "data/natural/h3n2/strains_metadata.tsv",
        sep="\t",
        parse_dates=["date"],
        usecols=["strain", "region", "date"],
    ).query(
        "(region != '?') & (date >= '2005-04-01') & (date <= '2019-10-01') & ('X' not in date)"
    )

    subsampled_metadata["year"] = subsampled_metadata["date"].dt.year
    subsampled_metadata["month"] = subsampled_metadata["date"].dt.month

    #subsampled_metadata["year_month_date"] = subsampled_metadata.apply(lambda row: f"{row['year']}-{row['month']:02}-01", axis=1)
    subsampled_metadata["year_month_date"] = subsampled_metadata.apply(lambda row: f"{row['year']}-01-01", axis=1)
    return (subsampled_metadata,)


@app.cell
def _(subsampled_metadata):
    subsampled_metadata.shape[0]
    return


@app.cell
def _(np, subsampled_metadata):
    subsampled_count_by_region_date = subsampled_metadata.groupby([
        "region",
        "year_month_date"
    ])["strain"].count().reset_index().rename(
        columns={"strain": "count"}
    )
    subsampled_count_by_region_date["log_count"] = np.log2(subsampled_count_by_region_date["count"])
    return (subsampled_count_by_region_date,)


@app.cell
def _(subsampled_count_by_region_date):
    subsampled_count_by_region_date
    return


@app.cell
def _(mo):
    mo.md(r"""### High-density subsampling""")
    return


@app.cell
def _(pd):
    high_subsampled_metadata = pd.read_csv(
        "data/natural/h3n2_high_density/strains_metadata.tsv",
        sep="\t",
        parse_dates=["date"],
        usecols=["strain", "region", "date"],
    ).query(
        "(region != '?') & (date >= '2005-04-01') & (date <= '2019-10-01') & ('X' not in date)"
    )

    high_subsampled_metadata["year"] = high_subsampled_metadata["date"].dt.year
    high_subsampled_metadata["month"] = high_subsampled_metadata["date"].dt.month

    #subsampled_metadata["year_month_date"] = subsampled_metadata.apply(lambda row: f"{row['year']}-{row['month']:02}-01", axis=1)
    high_subsampled_metadata["year_month_date"] = high_subsampled_metadata.apply(lambda row: f"{row['year']}-01-01", axis=1)
    return (high_subsampled_metadata,)


@app.cell
def _(high_subsampled_metadata, np):
    high_subsampled_count_by_region_date = high_subsampled_metadata.groupby([
        "region",
        "year_month_date"
    ])["strain"].count().reset_index().rename(
        columns={"strain": "count_high_subsampled"}
    )
    high_subsampled_count_by_region_date["log_count_high_subsampled"] = np.log2(
        high_subsampled_count_by_region_date["count_high_subsampled"]
    )
    return (high_subsampled_count_by_region_date,)


@app.cell
def _(high_subsampled_count_by_region_date):
    high_subsampled_count_by_region_date
    return


@app.cell
def _(mo):
    mo.md(r"""## Plot total and subsampled records per region""")
    return


@app.cell
def _(
    full_count_by_region_date,
    high_subsampled_count_by_region_date,
    subsampled_count_by_region_date,
):
    merged_count_by_region_date = full_count_by_region_date.merge(
        subsampled_count_by_region_date,
        on=["region", "year_month_date"],
        how="outer",
        suffixes=["_full", "_subsampled"],
    ).merge(
        high_subsampled_count_by_region_date,
        on=["region", "year_month_date"],
        how="outer",
    ).fillna(0)

    merged_count_by_region_date["proportion_subsampled"] = (
        merged_count_by_region_date["count_subsampled"] / merged_count_by_region_date["count_full"]
    )

    merged_count_by_region_date["proportion_high_subsampled"] = (
        merged_count_by_region_date["count_high_subsampled"] / merged_count_by_region_date["count_full"]
    )
    return (merged_count_by_region_date,)


@app.cell
def _(merged_count_by_region_date):
    merged_count_by_region_date
    return


@app.cell
def _(alt, merged_count_by_region_date):
    _full_chart = alt.Chart(merged_count_by_region_date).mark_line().encode(
        x=alt.X("year_month_date:T", title="date"),
        y=alt.Y("count_full:Q", title="number of records per year"),
        color="region:N",
    )

    _subsampled_chart = alt.Chart(merged_count_by_region_date).mark_line(
        strokeDash=[4, 2],
    ).encode(
        x=alt.X("year_month_date:T", title="date"),
        y=alt.Y("count_subsampled:Q", title="number of records per year"),
        color=alt.ColorValue("black"),
        opacity=alt.OpacityValue(0.7),
    )

    _high_subsampled_chart = alt.Chart(merged_count_by_region_date).mark_line(
        strokeDash=[4, 2],
    ).encode(
        x=alt.X("year_month_date:T", title="date"),
        y=alt.Y("count_high_subsampled:Q", title="number of records per year"),
        color=alt.ColorValue("gray"),
        opacity=alt.OpacityValue(0.7),
    )

    (_full_chart + _subsampled_chart + _high_subsampled_chart).facet(
        "region:N",
        columns=5,
    ).resolve_scale(y="independent")
    return


@app.cell
def _(alt, full_count_by_region_date, merged_count_by_region_date):
    width = 800
    height= 200

    _full_chart = alt.Chart(full_count_by_region_date).mark_line(point=True).encode(
        x=alt.X("year_month_date:T", title="date"),
        y=alt.Y("count:Q", title="number of records per year"),
        color="region:N",
    ).properties(
        width=width,
        height=height,
        title="All sequences",
    )

    _proportion_chart = alt.Chart(merged_count_by_region_date).mark_line(point=True).encode(
        x=alt.X("year_month_date:T", title="date"),
        y=alt.Y("proportion_subsampled:Q", title="proportion of records sampled per year"),
        color="region:N",
    ).properties(
        width=width,
        height=height,
        title="Subsampled sequences (90 per month)",
    )

    _high_proportion_chart = alt.Chart(merged_count_by_region_date).mark_line(point=True).encode(
        x=alt.X("year_month_date:T", title="date"),
        y=alt.Y("proportion_high_subsampled:Q", title="proportion of records sampled per year"),
        color="region:N",
    ).properties(
        width=width,
        height=height,
        title="Subsampled sequences (270 per month)",
    )

    _full_chart & _proportion_chart & _high_proportion_chart
    return


@app.cell
def _(merged_count_by_region_date):
    counts_by_region = merged_count_by_region_date.groupby(["region"]).aggregate(
        count_full=("count_full", "sum"),
        count_subsampled=("count_subsampled", "sum"),
        count_high_subsampled=("count_high_subsampled", "sum"),
    )
    counts_by_region["proportion_subsampled"] = counts_by_region["count_subsampled"] / counts_by_region["count_full"]
    counts_by_region["proportion_high_subsampled"] = counts_by_region["count_high_subsampled"] / counts_by_region["count_full"]
    return (counts_by_region,)


@app.cell
def _(counts_by_region):
    counts_by_region["proportion_subsampled"].round(2).sort_values()
    return


@app.cell
def _(counts_by_region):
    counts_by_region.loc[
        :,
        ["proportion_subsampled", "proportion_high_subsampled"]
    ].round(2).sort_values(
        "proportion_subsampled"
    )
    return


@app.cell
def _(counts_by_region):
    counts_by_region["proportion_subsampled"].mean()
    return


@app.cell
def _(counts_by_region):
    counts_by_region["proportion_subsampled"].median()
    return


@app.cell
def _(counts_by_region):
    counts_by_region["proportion_high_subsampled"].mean()
    return


@app.cell
def _(counts_by_region):
    counts_by_region["proportion_high_subsampled"].median()
    return


@app.cell
def _():
    return


if __name__ == "__main__":
    app.run()
