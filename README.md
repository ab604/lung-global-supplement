
Last Updated on 2024-07-25

# Comparative analysis of transcriptomic and proteomic expression reveals distinct molecular signatures between two non-small cell lung cancer subtypes: Supplementary Material

## Gene counts

### Hisat2

                     File                                    
                     Table-S1-Hisat-LUAD-vs-PBMC-counts.csv  
                     Table-S2-Hisat-LUSC-vs-PBMC-counts.csv  
                     Table-S3-Hisat-LUSC-vs-LUAD-counts.csv  

Column names: File

| Column name | Description |
|--------------------------------------------------|----------------------|
| `name` | Ensembl gene identifier |
| `gene` | HGNC gene symbol |
| `sample_id` the donor id or donor id suffixed with `T` for tumour or `N` for PBMC samples. | mapped read counts from featureCounts |

Hisat2 Counts Tables Information

### Salmon

                    File                                     
                    Table-S4-Salmon-LUAD-vs-PBMC-counts.csv  
                    Table-S5-Salmon-LUSC-vs-PBMC-counts.csv  
                    Table-S6-Salmon-LUSC-vs-LUAD-counts.csv  

Column names: File

| Column name | Description |
|-------------------------------------------------|-----------------------|
| `name` | Ensembl gene identifier |
| `gene` | HGNC gene symbol or Ensembl if missing |
| `length` | the length of the target transcript |
| `sample_id` the donor id or donor id suffixed with `T` for tumour or `N` for PBMC samples. | mapped reads counts by Salmon |

Salmon Counts Table Information

## EdgeR

Robinson, McCarthy, and Smyth (2009)

Tables S7-9. These results were filtered for common genes DE from both
Salmon and Hisat2 counts:

                      File                                 
                      Table-S7-edgeR-DEG-LUAD-vs-PBMC.csv  
                      Table-S8-edgeR-DEG-LUSC-vs-PBMC.csv  
                      Table-S9-edgeR-DEG-LUSC-vs-LUAD.csv  

Column names: File

| Column name | Description |
|-------------------------------------------------|-----------------------|
| `name` | Ensembl gene identifier |
| `gene` | HGNC gene symbol or Ensembl if missing |
| `baseMean` | mean read counts |
| `baseMeanA` | mean read count group A |
| `baseMeanB` | mean read count group B |
| `foldChange` | fold change B/A |
| `log2FoldChange` | log2 fold change B/A |
| `PValue` | p-value |
| `PAdj` | Benjamini-Hochbergadjusted p-value |
| `FDR` | False discovery rate |
| `falsePos` | false discovery counts |
| `sample_id` the donor id or donor id suffixed with `T` for tumour or `N` for PBMC samples. | sample Hisat2 read count |

edgeR Table information

## Peaks normalised Top 3 peptide intensities

Tables S10-12

           File                                                        
           Table-S10-Peaks-top3-peptides-intensities-LUAD-vs-NAT.csv   
           Table-S11-Peaks-top3-peptides-intensities-LUSC-vs-NAT.csv   
           Table-S12-Peaks-top3-peptides-intensities-LUSC-vs-LUAD.csv  

Column names: File

| Column name | Description |
|-----------------------------------------------|-------------------------|
| `protein` | protein short name |
| `gene` | HGNC gene symbol |
| `sample_id` the donor id or donor id suffixed with `T` for tumour or `N` for NAT samples | Normalised top 3 peptide intensity from Peaks |

Peaks normalised Top 3 peptide intensities Table information

## DEqMS

Tables S13-15 “DEqMS” (n.d.)

                      File                                  
                      Table-S13-DEqMS-DEP-LUAD-vs-NAT.csv   
                      Table-S14-DEqMS-DEP-LUSC-vs-NAT.csv   
                      Table-S15-DEqMS-DEP-LUSC-vs-LUAD.csv  

Column names: File

| Column name    | Description                                    |
|----------------|------------------------------------------------|
| `logFC`        | log2 fold change between two groups            |
| `AveExpr`      | the mean of the log2 ratios across all samples |
| `t`            | Limma t-values                                 |
| `P.Value`      | Limma p-values                                 |
| `adj.P.Val`    | BH method adjusted Limma p-values              |
| `B`            | Limma B values                                 |
| `gene`         | HGNC gene symbol                               |
| `count`        | peptide count values                           |
| `sca.t`        | DEqMS t-statistics                             |
| `sca.P.Value`  | DEqMS p-values                                 |
| `sca.adj.pval` | BH method adjusted DEqMS p-values              |
| `protein`      | protein short name                             |

DEqMS Table information

## g:Profiler

Tables S16-23 Liis Kolberg (2019)

                    File                                      
                    Table-S16-gprofiler-DEG-LUAD-vs-PBMC.csv  
                    Table-S17-gprofiler-DEG-LUSC-vs-PBMC.csv  
                    Table-S18-gprofiler-DEG-LUAD-vs-LUSC.csv  
                    Table-S19-gprofiler-DEG-LUSC-vs-LUAD.csv  
                    Table-S20-gprofiler-DEP-LUAD-vs-NAT.csv   
                    Table-S21-gprofiler-DEP-LUSC-vs-NAT.csv   
                    Table-S22-gprofiler-DEP-LUAD-vs-LUSC.csv  
                    Table-S23-gprofiler-DEP-LUSC-vs-LUAD.csv  

Column names: File

<table>
<caption>g:Profiler Table information</caption>
<colgroup>
<col style="width: 15%" />
<col style="width: 84%" />
</colgroup>
<thead>
<tr class="header">
<th>Column name</th>
<th>Description</th>
</tr>
</thead>
<tbody>
<tr class="odd">
<td><code>query</code></td>
<td>the name of the input query</td>
</tr>
<tr class="even">
<td><code>significant</code></td>
<td>indicator for statistically significant results</td>
</tr>
<tr class="odd">
<td><code>p_value</code></td>
<td>hypergeometric p-value after correction for multiple testing</td>
</tr>
<tr class="even">
<td><code>term_size</code></td>
<td>number of genes that are annotated to the term</td>
</tr>
<tr class="odd">
<td><code>query_size</code></td>
<td>number of genes that were included in the query</td>
</tr>
<tr class="even">
<td><code>intersection_size</code></td>
<td>the number of genes in the input query that are annotated to the
corresponding term</td>
</tr>
<tr class="odd">
<td><code>precision</code></td>
<td>the proportion of genes in the input list that are annotated to the
function (defined as intersection_size/query_size)</td>
</tr>
<tr class="even">
<td><code>recall</code></td>
<td>the proportion of functionally annotated genes that the query
recovers (defined as intersection_size/term_size)</td>
</tr>
<tr class="odd">
<td><code>term_id</code></td>
<td>unique term identifier</td>
</tr>
<tr class="even">
<td><p><code>source</code></p>
<p><code>term_name</code></p></td>
<td>the abbreviation of the data source for the term (e.g. GO:BP)</td>
</tr>
<tr class="odd">
<td><code>effective_domain_size</code></td>
<td>the total number of genes “in the universe” used for the
hypergeometric test</td>
</tr>
<tr class="even">
<td><code>source_order</code></td>
<td>numeric order for the term within its data source</td>
</tr>
<tr class="odd">
<td><code>parents</code></td>
<td>list of term IDs that are hierarchically directly above the term.
For non-hierarchical data sources this points to an artificial root
node.</td>
</tr>
<tr class="even">
<td><code>evidence_codes</code></td>
<td>a lists of all evidence codes for the intersecting genes between
input and the term. The evidences are separated by comma for each
gene.</td>
</tr>
<tr class="odd">
<td><code>intersection</code></td>
<td>a comma separated list of genes from the query that are annotated to
the corresponding term</td>
</tr>
</tbody>
</table>

g:Profiler Table information

## References

“DEqMS.” n.d. <http://bioconductor.org/packages/DEqMS/>.

Liis Kolberg, Uku Raudvere. 2019. “Gprofiler2: Interface to the
’g:profiler’ Toolset.” <https://biit.cs.ut.ee/gprofiler/page/r>.

Robinson, Mark D., Davis J. McCarthy, and Gordon K. Smyth. 2009. “edgeR:
A Bioconductor Package for Differential Expression Analysis of Digital
Gene Expression Data.” *Bioinformatics* 26 (1): 139–40.
<https://doi.org/10.1093/bioinformatics/btp616>.
