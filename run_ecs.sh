#!/bin/bash
# Qun Liu, Univ of Exeter

calc="ULI_MEDIUM_SAMPLE"
inpath='./data'
outpath='./data'
ref_papers="sherwood myers"

for ref_paper in $ref_papers
do
    echo $ref_paper
    log_fn=$outpath/out_"$ref_paper"_"$calc".txt 

    if [[ -f $log_fn ]]
    then
        echo "$log_fn already present!"
    else
        python -u ./scripts/ecs-baseline-ec2.py $calc $outpath $inpath $ref_paper &> $log_fn
        ./makelinks.sh $outpath/$ref_paper
    fi
    # ./makelinks.sh $outpath/$ref_paper
done

python -u ./scripts/reproduce_myers2021.py

