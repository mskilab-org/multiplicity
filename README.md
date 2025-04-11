
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
      <td style="border: 1px solid black; padding: 5px;">Path to somatic vcf. </td>
    </tr>
    <tr>
      <td style="border: 1px solid black; padding: 5px;"><code>germline_snv</code></td>
      <td style="border: 1px solid black; padding: 5px;"><code>NULL</code></td>
      <td style="border: 1px solid black; padding: 5px;">Path to germline vcf.</td>
    </tr>
    <tr>
      <td style="border: 1px solid black; padding: 5px;"><code>jabba_rds</code></td>
      <td style="border: 1px solid black; padding: 5px;"></td>
      <td style="border: 1px solid black; padding: 5px;">JaBbA-derived
  gGraph or jabba.simple.rds</td>
    </tr>
    <tr>
      <td style="border: 1px solid black; padding: 5px;"><code>snpeff_path</code></td>
      <td style="border: 1px solid black; padding: 5px;"><code>/path/to/snpeff</code></td>
      <td style="border: 1px solid black; padding: 5px;">Path to unzipped/untarred
  SnpEff directory, containing discrete scripts; importantly, SnpSift & vcfEffOnePerLine.pl</td>
    </tr>
    <tr>
      <td style="border: 1px solid black; padding: 5px;"><code>tumor_name</code></td>
      <td style="border: 1px solid black; padding: 5px;"><code>NULL</code></td>
      <td style="border: 1px solid black; padding: 5px;"></td>
    </tr>
    <tr>
      <td style="border: 1px solid black; padding: 5px;"><code>normal_name</code></td>
      <td style="border: 1px solid black; padding: 5px;"><code>NULL</code></td>
      <td style="border: 1px solid black; padding: 5px;"></td>
    </tr>
    <tr>
      <td style="border: 1px solid black; padding: 5px;"><code>verbose</code></td>
      <td style="border: 1px solid black; padding: 5px;"><code>TRUE</code></td>
      <td style="border: 1px solid black; padding: 5px;">Informative logging of function?</td>
    </tr>
  </tbody>
</table>


## License

[MIT](https://choosealicense.com/licenses/mit/)
