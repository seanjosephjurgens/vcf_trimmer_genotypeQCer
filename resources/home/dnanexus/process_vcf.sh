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
    dx download "$1"

    # 2. Strip original path and define local input and output files
    IFS='/' read -r -a f <<< "$1"
    FILEIN=${f[-1]}
    FILEOUT="${FILEIN%.vcf.gz}_"$5".vcf.gz"

    # 3. Run bcftools accounting for options to include
    if [ "$2" != "NA" ] && [ "$3" != "NA" ]
    then
        export BCFTOOLS_PLUGINS=plugins/
        bcftools annotate -x $2 $FILEIN -Ou --max-mem 8G | bcftools norm -m - -Ou --max-mem 8G | bcftools view -f PASS -Ou  --max-mem 8G | bcftools +setGT -Ou --max-mem 8G -- -t q -i 'FORMAT/FT!="PASS" & FORMAT/FT!="."' -n . | bcftools filter -i $3 -Oz --max-mem 8G -o $FILEOUT --threads $4
        #bcftools annotate -x $2 $FILEIN | bcftools norm -m - | bcftools view -f PASS | bcftools filter -i $3 -Oz -o $FILEOUT --threads $4
        #export BCFTOOLS_PLUGINS=plugins/
        #bcftools +setGT $FILEOUT --threads $4 -- -t q -i 'FORMAT/FT!="PASS"' -n .
    elif [ "$2" != "NA" ] && [ "$3" == "NA" ]
    then
        export BCFTOOLS_PLUGINS=plugins/
        bcftools annotate -x $2 $FILEIN -Ou --max-mem 8G | bcftools norm -m - -Ou  --max-mem 8G | bcftools view -f PASS -Ou  --max-mem 8G | bcftools +setGT -Oz -o $FILEOUT --threads $4 --max-mem 8G -- -t q -i 'FORMAT/FT!="PASS" & FORMAT/FT!="."' -n .
        #bcftools annotate -x $2 $FILEIN | bcftools norm -m - | bcftools view -f PASS | bcftools -Oz -o $FILEOUT --threads $4
        #export BCFTOOLS_PLUGINS=plugins/
        #bcftools +setGT $FILEOUT --threads $4 -- -t q -i 'FORMAT/FT!="PASS"' -n .
    else
        export BCFTOOLS_PLUGINS=plugins/
        bcftools norm -m - $FILEIN -Ou --max-mem 8G | bcftools view -f PASS -Ou --max-mem 8G | bcftools +setGT -Ou --max-mem 8G -- -t q -i 'FORMAT/FT!="PASS" & FORMAT/FT!="."' -n . | bcftools filter -i $3 -Oz -o $FILEOUT --threads $4 --max-mem 8G 
        #bcftools norm -m - $FILEIN | bcftools view -f PASS | bcftools filter -i $3 -Oz -o $FILEOUT --threads $4
        #export BCFTOOLS_PLUGINS=plugins/
        #bcftools +setGT $FILEOUT --threads $4 -- -t q -i 'FORMAT/FT!="PASS"' -n .
    fi

    # 4. Upload trimmed VCF here:
    dx upload $FILEOUT --path "$6/" --brief

    # 5. Clear the worker of input and output files to reserve storage
    rm $FILEIN $FILEOUT

else
    echo "Ignoring VCF as no filtering or removal of data fields requested!"
fi


