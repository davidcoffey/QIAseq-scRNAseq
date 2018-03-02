#!/bin/bash
# Extract antigen receptor sequences using MiXCR
# Written by David Coffey dcoffey@fhcrc.org
# Updated December 15, 2017

## Prerequisites (see Software_installation.sh)
# Install MiXCR

## Variables
# export SAMPLE="..."
# export MIXCR_DIRECTORY=".../MiXCR"
# export READ1="..."
# export READ2="..."

## Reference: https://mixcr.readthedocs.io/en/latest/rnaseq.html

START=`date +%s`
echo Begin MiXCR.sh for sample $SAMPLE cell $CELL on `date +"%B %d, %Y at %r"`

module load mixcr/2.1.5-foss-2016b

LOCUS="IGH IGK IGL TRB TRA"

for L in ${LOCUS}; do
    echo ${L} $SAMPLE.$CELL
    mkdir -p $MIXCR_DIRECTORY/${L}

    # Align
    mixcr align \
    --chains ${L} \
    --species hsa \
    -f \
    -p rna-seq \
    -OallowPartialAlignments=true \
    $READ1 \
    $MIXCR_DIRECTORY/${L}/$SAMPLE.$CELL.${L}.vdjca

    # Assemble
    mixcr assemblePartial \
    $MIXCR_DIRECTORY/${L}/$SAMPLE.$CELL.${L}.vdjca \
    $MIXCR_DIRECTORY/${L}/$SAMPLE.$CELL.${L}.1.vdjca

    mixcr assemblePartial \
    $MIXCR_DIRECTORY/${L}/$SAMPLE.$CELL.${L}.1.vdjca \
    $MIXCR_DIRECTORY/${L}/$SAMPLE.$CELL.${L}.2.vdjca

    mixcr assemble \
    $MIXCR_DIRECTORY/${L}/$SAMPLE.$CELL.${L}.2.vdjca \
    $MIXCR_DIRECTORY/${L}/$SAMPLE.$CELL.${L}.clns

    # Export
    mixcr exportClones \
    --chains ${L} \
    -cloneId \
    -aaFeature CDR3 \
    -nFeature CDR3 \
    -count \
    -fraction \
    -vHit \
    -dHit \
    -jHit \
    -cHit \
    -vFamily \
    -dFamily \
    -jFamily \
    -cFamily \
    -vHitScore \
    -dHitScore \
    -jHitScore \
    -cHitScore \
    -vAlignment \
    -dAlignment \
    -jAlignment \
    -cAlignment \
    -vBestIdentityPercent \
    -dBestIdentityPercent \
    -jBestIdentityPercent \
    -cBestIdentityPercent \
    -qFeature CDR3 \
    -avrgFeatureQuality CDR3 \
    -minFeatureQuality CDR3 \
    -lengthOf CDR3 \
    $MIXCR_DIRECTORY/${L}/$SAMPLE.$CELL.${L}.clns $MIXCR_DIRECTORY/${L}/$SAMPLE.$CELL.${L}.txt

    # Rename columns
    sed -i 's/Clone ID/cloneID/g' $MIXCR_DIRECTORY/${L}/$SAMPLE.$CELL.${L}.txt
    sed -i 's/aaSeqCDR3/aminoAcid/g' $MIXCR_DIRECTORY/${L}/$SAMPLE.$CELL.${L}.txt
    sed -i 's/nSeqCDR3/nucleotide/g' $MIXCR_DIRECTORY/${L}/$SAMPLE.$CELL.${L}.txt
    sed -i 's/cloneCount/count/g' $MIXCR_DIRECTORY/${L}/$SAMPLE.$CELL.${L}.txt
    sed -i 's/cloneFraction/frequencyCount/g' $MIXCR_DIRECTORY/${L}/$SAMPLE.$CELL.${L}.txt
    sed -i 's/bestVHit/vGeneName/g' $MIXCR_DIRECTORY/${L}/$SAMPLE.$CELL.${L}.txt
    sed -i 's/bestDHit/dGeneName/g' $MIXCR_DIRECTORY/${L}/$SAMPLE.$CELL.${L}.txt
    sed -i 's/bestJHit/jGeneName/g' $MIXCR_DIRECTORY/${L}/$SAMPLE.$CELL.${L}.txt
    sed -i 's/bestCHit/cGeneName/g' $MIXCR_DIRECTORY/${L}/$SAMPLE.$CELL.${L}.txt
    sed -i 's/bestVFamily/vFamilyName/g' $MIXCR_DIRECTORY/${L}/$SAMPLE.$CELL.${L}.txt
    sed -i 's/bestDFamily/dFamilyName/g' $MIXCR_DIRECTORY/${L}/$SAMPLE.$CELL.${L}.txt
    sed -i 's/bestJFamily/jFamilyName/g' $MIXCR_DIRECTORY/${L}/$SAMPLE.$CELL.${L}.txt
    sed -i 's/bestCFamily/cFamilyName/g' $MIXCR_DIRECTORY/${L}/$SAMPLE.$CELL.${L}.txt

    mv $MIXCR_DIRECTORY/${L}/$SAMPLE.$CELL.${L}.txt $MIXCR_DIRECTORY/${L}/$SAMPLE.$CELL.${L}.tsv

    # Remove temporary files
    rm $MIXCR_DIRECTORY/${L}/$SAMPLE.$CELL.${L}.clns
    rm $MIXCR_DIRECTORY/${L}/$SAMPLE.$CELL.${L}.vdjca
    rm $MIXCR_DIRECTORY/${L}/$SAMPLE.$CELL.${L}.1.vdjca
    rm $MIXCR_DIRECTORY/${L}/$SAMPLE.$CELL.${L}.2.vdjca
    
    # Remove files with less than 2 lines
    find $MIXCR_DIRECTORY/${L}/$SAMPLE.$CELL.${L}.tsv -type f -exec awk -v x=2 'NR==x{exit 1}' {} \; -exec rm -f {} \;
    
done

END=`date +%s`
MINUTES=$(((END-START)/60))
echo End MiXCR.sh for sample $SAMPLE cell $CELL.  The run time was $MINUTES minutes.