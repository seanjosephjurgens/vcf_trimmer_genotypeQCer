#!/bin/bash
# vcf_trimmer 0.0.1
# Andrew R Wood
# University of Exeter
# No warrenty provided!

# $1 = string containing full file path of VCF on DNANexus RAP
# $2 = string containing fields to remove from VCF
# $3 = string containing fields and threholds to apply for purposes of QC
# $4 = number of threads for compressed file writing by bcftools
# $5 = label to add to vcf processed
# $6 = output directory

# Set error catching
set -euo pipefail

# Check user has provided either fields to remove or metric thresholding
if [ "$2" != "NA" ] || [ "$3" != "NA" ]
then

    # 1. Get the file onto the worker
    echo "Processing $1"
    dx download -a -f "$1"

    # 2. Strip original path and define local input and output files
    IFS='/' read -r -a f <<< "$1"
    FILEIN=${f[-1]}
    #FILEINTER="${FILEIN%.vcf.gz}_inter.bcf"
    FILEOUT="${FILEIN%.vcf.gz}_"$5".vcf.gz"

    # 3. Run bcftools accounting for options to include
    if [ "$2" != "NA" ] && [ "$3" != "NA" ]
    then
        export BCFTOOLS_PLUGINS=plugins/
        bcftools annotate -x $2 $FILEIN -Ou | bcftools +setGT -Ou -- -t q -i 'FORMAT/FT!="PASS" & FORMAT/FT!="."' -n . | bcftools annotate -x FORMAT/FT -Ou | bcftools filter -i $3 -Oz -o $FILEOUT --threads $4
    elif [ "$2" != "NA" ] && [ "$3" == "NA" ]
    then
        export BCFTOOLS_PLUGINS=plugins/
        bcftools annotate -x $2 $FILEIN -Ou | bcftools +setGT -Ou -- -t q -i 'FORMAT/FT!="PASS" & FORMAT/FT!="."' -n . | bcftools annotate -x FORMAT/FT -Oz -o $FILEOUT --threads $4 
    else
        export BCFTOOLS_PLUGINS=plugins/
        bcftools +setGT $FILEINTER -Ou -- -t q -i 'FORMAT/FT!="PASS" & FORMAT/FT!="."' -n . | bcftools filter -i $3 -Oz -o $FILEOUT --threads $4
        #bcftools norm -m - $FILEIN -Ob | bcftools view -f PASS -Ob | bcftools filter -i $3 -Oz -o $FILEINTER --threads $4
    fi
    #export BCFTOOLS_PLUGINS=plugins/
    #bcftools +setGT $FILEINTER -Oz -o $FILEOUT --threads $4 -- -t q -i 'FORMAT/FT!="PASS" & FORMAT/FT!="."' -n .

    # 4. Upload trimmed VCF here:
    dx upload $FILEOUT --path "$6/" --brief

    # 5. Clear the worker of input and output files to reserve storage
    #rm $FILEIN $FILEINTER $FILEOUT
    rm $FILEIN $FILEOUT

else
    echo "Ignoring VCF as no filtering or removal of data fields requested!"
fi


