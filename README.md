## Script to perfom TP53-PADI4 expression analyses through public data cohorts

## 1. datasets:
### 1.1. TCGA
tcga_data/ contains TCGA datasets except for the RNA sequencing file **which should be downloaded (as the file size exceeds GitHub max)** to tcga_data/ from UCEC Xena (https://xenabrowser.net) through:
'https://toil-xena-hub.s3.us-east-1.amazonaws.com/download/tcga_RSEM_Hugo_norm_count.gz'
name this file: tcga_data/tcga_RSEM_Hugo_norm_count.txt

### 1.2 melanoma immunotherapy
ici_data/ contains anti-PD1 datasets


## 2. script to reproduce the results
tp53_padi4_analysis.py

System: tested with python 3.7 Linux and MACOS

Dependencies: pickle, pandas, numpy, lifelines, scipy, matplotlib

Run:

```python
import tp53_padi4_analysis.main()
```

to reproduce analyses. Results wll be saved to outs/ directory.
