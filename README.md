# Cancer Classification with High‑Dimensional and Imbalanced Gene Expression Data

This repository contains R source code for the manuscript *Cancer Classification with High‑Dimensional and Imbalanced Gene Expression Data*.

## Paper abstract
Recognizing cancer types and molecular phenotypes is important for characterizing disease heterogeneity and informing patient assessment. Expression data support this task, but unequal class sizes and numerous molecular features complicate identification of less frequent groups. This study examined three cancer expression datasets with progressively greater class imbalance, comparing six sampling conditions and three ensemble classifiers. Fisher and ReliefF were combined to select features, and the effects of sampling and candidate feature count were evaluated.

**Journal Manuscript**: Manuscript_Corrected.docx

## Datasets
All datasets are publicly available from NCBI GEO:
- D1 breast HER2: [GSE96058](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE96058)
- D2 lung serum: [GSE106817](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE106817)
- D3 sarcoma serum: [GSE211692](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE211692)

> Note: This repo **does not store raw GEO expression data**. You need to download raw data from GEO website.
> `Supplementary_Data.xlsx` includes sample labels, train‑validation‑test partitions, hyper‑parameters and full experiment tables.

## Environment & Dependencies
R version >=4.2.0
Required packages:
```r
# install packages
install.packages(c("randomForest","gbm","xgboost","FSelector","CORElearn","pROC","openxlsx"))
