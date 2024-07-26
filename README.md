# VCF Trimmer for UK Biobank RAP
#### Developed by Andrew Wood. University of Exeter
#### Re-purposed by Sean Jurgens. Broad Institute / Amsterdam UMC

This applet performs parallel processing of VCFs through [bcftools](https://samtools.github.io/bcftools/bcftools.html). This applet has been designed to reduce file sizes down when required for merging of data by:
* Splitting multiallelic sites (note left-alignment and normalisation not presently set)
* Removing fields within `INFO` and/or `FORMAT`
* Applies inclusion quality control thresholding based on field in `INFO`

---
### Obtaining and installing the applet
Clone this github repo to a local directory:
```
git clone https://github.com/drarwood/vcf_trimmer
```
Navigate to a relevant directory within the project directory on the DNAnexus platform
```
dx cd /path/to/install/apps
```
Now you are ready to build and upload the applet to the DNAnexus platform directory
```
dx build -f vcf_trimmer
```
---
### Inputs
This applet takes one file in as input that simply lists the VCFs to process, one per line. VCFs will be processed in the order as the appear in the file. For example:
```
/Bulk/Whole genome sequences/Population level WGS variants, pVCF format - interim 200k release/ukb24304_c1_b0_v1.vcf.gz
/Bulk/Whole genome sequences/Population level WGS variants, pVCF format - interim 200k release/ukb24304_c1_b1_v1.vcf.gz
/Bulk/Whole genome sequences/Population level WGS variants, pVCF format - interim 200k release/ukb24304_c1_b2_v1.vcf.gz
/Bulk/Whole genome sequences/Population level WGS variants, pVCF format - interim 200k release/ukb24304_c1_b3_v1.vcf.gz
/Bulk/Whole genome sequences/Population level WGS variants, pVCF format - interim 200k release/ukb24304_c1_b4_v1.vcf.gz
/Bulk/Whole genome sequences/Population level WGS variants, pVCF format - interim 200k release/ukb24304_c1_b5_v1.vcf.gz
/Bulk/Whole genome sequences/Population level WGS variants, pVCF format - interim 200k release/ukb24304_c1_b6_v1.vcf.gz
/Bulk/Whole genome sequences/Population level WGS variants, pVCF format - interim 200k release/ukb24304_c1_b7_v1.vcf.gz
/Bulk/Whole genome sequences/Population level WGS variants, pVCF format - interim 200k release/ukb24304_c1_b8_v1.vcf.gz
/Bulk/Whole genome sequences/Population level WGS variants, pVCF format - interim 200k release/ukb24304_c1_b9_v1.vcf.gz
```
---
### Command line usage
Strings that define inclusion criteria and fields to exclude should be consistent with input expected by bcftools and placed within quotes. See the [bcftools](https://samtools.github.io/bcftools/bcftools.html) manual for more information.

### IMPORTANT:
As opposed to the original version, this one is made for analyzing the DRAGEN 500k WGS data, and also performs the following QC by default:
`FILTER==PASS`
`QC>20`
`FT==PASS`

#### Example 1:
Removing all fields within `FORMAT` except for `GT`, `GQ` and `FT` (which are needed for genotypes and basic filtering):
```
dx run vcf_trimmer \
  -ivcf_file_list=/path/to/vcf_file_list.txt \
  -ifile_label=trimmed3 \
  -ioutput_dir=/path/to/output/dir \
  -iqc_thresholds="NA" \
  -ifields_to_remove="FORMAT/LAD,FORMAT/LPL,FORMAT/LAA,FORMAT/LAF,FORMAT/QL" \
  -y
```

#### Example 3:
Removing all fields within `FORMAT` (except for `GT`, `GQ` and `FT`) and removing all from `INFO`:

```
dx run vcf_trimmer \
  -ivcf_file_list=/path/to/vcf_file_list.txt \
  -ifile_label=trimmed2 \
  -ioutput_dir=/path/to/output/dir \
  -iqc_thresholds="NA" \
  -ifields_to_remove="FORMAT/GQ,FORMAT/LAD,FORMAT/FT,FORMAT/LPL,FORMAT/LAA,FORMAT/LAF,FORMAT/QL,INFO/AC,INFO/AN,INFO/NS,INFO/NS_GT,INFO/NS_NOGT,INFO/NS_NODATA,INFO/IC,INFO/HWE,INFO/ExcHet,INFO/HWE_CHISQ" \
  -y
```

#### Example 3:
Removing all fields within `FORMAT` (except for `GT`, `GQ` and `FT`) and removing all from `INFO`, except for `ExcHet` which must remain as as inclusion critera of an `ExcHet` > 0.5:

```
dx run vcf_trimmer \
  -ivcf_file_list=/path/to/vcf_file_list.txt \
  -ifile_label=trimmed2 \
  -ioutput_dir=/path/to/output/dir \
  -iqc_thresholds="INFO/ExcHet>0.5" \
  -ifields_to_remove="FORMAT/GQ,FORMAT/LAD,FORMAT/FT,FORMAT/LPL,FORMAT/LAA,FORMAT/LAF,FORMAT/QL,INFO/AC,INFO/AN,INFO/NS,INFO/NS_GT,INFO/NS_NOGT,INFO/NS_NODATA,INFO/IC,INFO/HWE,INFO/HWE_CHISQ" \
  -y
```

#### Additional parameters
Compute costs not yet intenstively testes; for previous version the below applied:

To assist in the throughput of this applet, multiple VCFs will be processed at the same time on a given workstation. 
The default instance type is `mem1_ssd1_v2_x36`. Benchmarking has been based on this server type for the 200K WGS data on chromosome 17 and the maximum number of concurrent VCFs processed at a time to `20` to avoid resource issues. 
For other datasets, you may need to lower this limit if jobs fail due to errors raised because of server response timeouts (e.g. 12). 
```
-iconcurrent_processes (default: 20) : maximum number of VCFs to process at a given time.
-ithreads (default: 4) : number of threads bcftools should use when writing compressed output
```


