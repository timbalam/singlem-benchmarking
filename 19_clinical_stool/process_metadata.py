import polars as pl

def read_clinical_stool_benchmark_metadata(file):
    return (
        pl.read_excel(file, sheet_name = "general data",
                      engine = "calamine",
                      read_options = {"header_row": 1},
                      columns = ["StudyID", "Submitter", "CultureSummary", "EBP_summary", "BF_summary"])
        .with_columns(
            Campy = pl.col("CultureSummary").is_in(["Campy", "Campy and STEC"]),
            Salm = pl.col("CultureSummary").is_in(["Salm"]),
            healthy = pl.col("CultureSummary").is_in(["healthy"]),
            neg = pl.col("CultureSummary").is_in(["neg"])
        )
    )