
```

                       888 888    d8b          888 d8b          d8b 888             
                       888 888    Y8P          888 Y8P          Y8P 888             
                       888 888                 888                  888             
88888b.d88b.  888  888 888 888888 888 88888b.  888 888  .d8888b 888 888888 888  888 
888 "888 "88b 888  888 888 888    888 888 "88b 888 888 d88P"    888 888    888  888 
888  888  888 888  888 888 888    888 888  888 888 888 888      888 888    888  888 
888  888  888 Y88b 888 888 Y88b.  888 888 d88P 888 888 Y88b.    888 Y88b.  Y88b 888 
888  888  888  "Y88888 888  "Y888 888 88888P"  888 888  "Y8888P 888  "Y888  "Y88888 
                                      888                                       888 
                                      888                                  Y8b d88P 
                                      888                                   "Y88P"  

```
## <font color=black> Introduction </font>


A simplistic transformation of counts to multiplicity as outlined in our
published works: JaBbA v1 [Choo et
al.](https://www.nature.com/articles/s41588-023-01540-6) & JaBbA v0 [Hadi et
al.](https://www.cell.com/cell/fulltext/S0092-8674(20)30997-1). Please see
Supplementary Note 2 and STAR Methods, respectively, in these publications for
full details on this derivation. 

`multiplicity` is an R library to convert counts from mutation-caller derived .vcfs
to total copies given known purity, ploidy, and consitutional discerned copy
number as ideally determined by [JaBbA](https://github.com/mskilab-org/JaBbA).

## <font color=black> Usage </font>

Either `somatic_snv` or `germline_snv` must be provided. Typically, mutation
callers will deposit alt/ref counts per tumor/normal pair in a single file. 

As such, we suggest you provide the expected define column names for the tumor and
normal sample (`tumor_name` & `normal_name`, respectively). By default, if a
single sample is contained in the provided vcf, that sample will be assumed to be
tumor. If two samples are contained in the provided vcf, the second column will
be defined to be the tumor whereas the first will be defined as the normal.

A hg19 mask is provided in testdata of this package. Please see our dryclean github for 
masks for hg38.


<table style="border: 1px solid black; border-collapse: collapse;">
  <tbody>
    <tr>
      <th style="border: 1px solid black; padding: 5px;">Parameter</th>
      <th style="border: 1px solid black; padding: 5px;">Default value</th>
      <th style="border: 1px solid black; padding: 5px;">Description/notes</th>
    </tr>
    <tr>
      <td style="border: 1px solid black; padding: 5px;"><code>somatic_snv</code></td>
      <td style="border: 1px solid black; padding: 5px;"><code>NULL</code></td>
      <td style="border: 1px solid black; padding: 5px;">Path to somatic SNV file.</td>
    </tr>
    <tr>
      <td style="border: 1px solid black; padding: 5px;"><code>germline_snv</code></td>
      <td style="border: 1px solid black; padding: 5px;"><code>NULL</code></td>
      <td style="border: 1px solid black; padding: 5px;">Path to germline SNV file.</td>
    </tr>
    <tr>
      <td style="border: 1px solid black; padding: 5px;"><code>het_pileups_wgs</code></td>
      <td style="border: 1px solid black; padding: 5px;"><code>NULL</code></td>
      <td style="border: 1px solid black; padding: 5px;">Path to het_pileups_wgs file.</td>
    </tr>
    <tr>
      <td style="border: 1px solid black; padding: 5px;"><code>jabba_rds</code></td>
      <td style="border: 1px solid black; padding: 5px;"><code>NULL</code></td>
      <td style="border: 1px solid black; padding: 5px;">Path to JaBbA file.</td>
    </tr>
    <tr>
      <td style="border: 1px solid black; padding: 5px;"><code>tumor_cbs</code></td>
      <td style="border: 1px solid black; padding: 5px;"><code>NULL</code></td>
      <td style="border: 1px solid black; padding: 5px;">Path to segmented drycleaned coverage file. Overrides <code>tumor_dryclean</code> if provided.</td>
    </tr>
    <tr>
      <td style="border: 1px solid black; padding: 5px;"><code>tumor_dryclean</code></td>
      <td style="border: 1px solid black; padding: 5px;"><code>NULL</code></td>
      <td style="border: 1px solid black; padding: 5px;">Path to drycleaned coverage file. If provided, REF/ALT values will be rescaled based on expected base coverage.</td>
    </tr>
    <tr>
      <td style="border: 1px solid black; padding: 5px;"><code>dryclean_field</code></td>
      <td style="border: 1px solid black; padding: 5px;"></td>
      <td style="border: 1px solid black; padding: 5px;">Field in tumor_dryclean gRanges used as rescaled binned coverage.</td>
    </tr>
    <tr>
      <td style="border: 1px solid black; padding: 5px;"><code>read_size</code></td>
      <td style="border: 1px solid black; padding: 5px;"><code>151</code></td>
      <td style="border: 1px solid black; padding: 5px;">Estimated average read size (default for Illumina sequencers).</td>
    </tr>
    <tr>
      <td style="border: 1px solid black; padding: 5px;"><code>tumor_name</code></td>
      <td style="border: 1px solid black; padding: 5px;"><code>NULL</code></td>
      <td style="border: 1px solid black; padding: 5px;">Expected tumor sample name as annotated in the VCF.</td>
    </tr>
    <tr>
      <td style="border: 1px solid black; padding: 5px;"><code>normal_name</code></td>
      <td style="border: 1px solid black; padding: 5px;"><code>NULL</code></td>
      <td style="border: 1px solid black; padding: 5px;">Expected normal sample name as annotated in the VCF.</td>
    </tr>
    <tr>
      <td style="border: 1px solid black; padding: 5px;"><code>mask</code></td>
      <td style="border: 1px solid black; padding: 5px;"><code>NULL</code></td>
      <td style="border: 1px solid black; padding: 5px;">Loci within provided mask gRanges will be excluded from analysis.</td>
    </tr>
    <tr>
      <td style="border: 1px solid black; padding: 5px;"><code>inferred_sex</code></td>
      <td style="border: 1px solid black; padding: 5px;"><code>NA</code></td>
      <td style="border: 1px solid black; padding: 5px;">Enum of ["M", "F"] indicating the inferred/true patient sex. If not provided, it will be inferred from coverage/JaBbA graph.</td>
    </tr>
    <tr>
      <td style="border: 1px solid black; padding: 5px;"><code>tau_in_gamma</code></td>
      <td style="border: 1px solid black; padding: 5px;"></td>
      <td style="border: 1px solid black; padding: 5px;">Use tau (TRUE) versus tau_hat (FALSE) in computation of gamma. Note: tau_hat uses the average copy number of loci, while tau is simply ploidy.</td>
    </tr>
    <tr>
      <td style="border: 1px solid black; padding: 5px;"><code>filterpass</code></td>
      <td style="border: 1px solid black; padding: 5px;"></td>
      <td style="border: 1px solid black; padding: 5px;">Process only variants with FILTER == PASS?</td>
    </tr>
    <tr>
      <td style="border: 1px solid black; padding: 5px;"><code>purity</code></td>
      <td style="border: 1px solid black; padding: 5px;"></td>
      <td style="border: 1px solid black; padding: 5px;">Purity of the input jabba_rds. Optional if gGraph metadata includes purity.</td>
    </tr>
    <tr>
      <td style="border: 1px solid black; padding: 5px;"><code>ploidy</code></td>
      <td style="border: 1px solid black; padding: 5px;"></td>
      <td style="border: 1px solid black; padding: 5px;">Ploidy of the input jabba_rds. Optional if gGraph metadata includes ploidy.</td>
    </tr>
    <tr>
      <td style="border: 1px solid black; padding: 5px;"><code>snpeff_path</code></td>
      <td style="border: 1px solid black; padding: 5px;"><code>/path/to/snpeff</code></td>
      <td style="border: 1px solid black; padding: 5px;">Path to the unzipped SnpEff toolkit.</td>
    </tr>
    <tr>
      <td style="border: 1px solid black; padding: 5px;"><code>verbose</code></td>
      <td style="border: 1px solid black; padding: 5px;"><code>TRUE</code></td>
      <td style="border: 1px solid black; padding: 5px;">Enable verbose output.</td>
    </tr>
  </tbody>
</table>


## License

[MIT](https://choosealicense.com/licenses/mit/)
