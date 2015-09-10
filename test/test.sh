#export PYTHONPATH=$PYTHONPATH:path/the/TFFM/package;
export PYTHONPATH=$PYTHONPATH:/raid6/amathelier/TFFM+DNAshape/bin/TFFM/;
araTha=/sharedro/RESOURCES/GBshape/araTha10;
helt=$araTha/araTha10.HelT.bigWig;
mgw=$araTha/araTha10.MGW.bigWig;
prot=$araTha/araTha10.ProT.bigWig;
roll=$araTha/araTha10.Roll.bigWig;
helt2=$araTha/araTha10.HelT.2nd.wig.bw;
mgw2=$araTha/araTha10.MGW.2nd.wig.bw;
prot2=$araTha/araTha10.ProT.2nd.wig.bw;
roll2=$araTha/araTha10.Roll.2nd.wig.bw;

echo "Training a TFFM + DNA shape classifier.";
time python2.7 ../DNAshapedTFBS.py trainTFFM -T TFFM_trained.xml \
    -i foreground/train.fa -I foreground/train.bed \
    -b background/train.fa -B background/train.bed \
    -o DNAshapedTFFM_classifier \
    -H $helt -E $helt2 -M $mgw -G $mgw2 -P $prot -O $prot2 -R $roll -L $roll2 \
    -n;

echo "Training a PSSM + DNA shape classifier.";
time python2.7 ../DNAshapedTFBS.py trainPSSM -f MA0563.1.pfm \
    -i foreground/train.fa -I foreground/train.bed \
    -b background/train.fa -B background/train.bed \
    -o DNAshapedPSSM_classifier \
    -H $helt -E $helt2 -M $mgw -G $mgw2 -P $prot -O $prot2 -R $roll -L $roll2 \
    -n;

echo "Applying the trained TFFM + DNA shape classifier on foreground sequences.";
time python2.7 ../DNAshapedTFBS.py applyTFFM -T TFFM_trained.xml \
    -i foreground/test.fa -I foreground/test.bed \
    -c DNAshapedTFFM_classifier.pkl -o DNAshapedTFFM_fg_predictions.txt \
    -H $helt -E $helt2 -M $mgw -G $mgw2 -P $prot -O $prot2 -R $roll -L $roll2 \
    -n;

echo "Applying the trained PSSM + DNA shape classifier on foreground sequences.";
time python2.7 ../DNAshapedTFBS.py applyPSSM -f MA0563.1.pfm \
    -i foreground/test.fa -I foreground/test.bed \
    -c DNAshapedPSSM_classifier.pkl -o DNAshapedPSSM_fg_predictions.txt \
    -H $helt -E $helt2 -M $mgw -G $mgw2 -P $prot -O $prot2 -R $roll -L $roll2 \
    -n;

echo "Applying the trained TFFM + DNA shape classifier on background sequences.";
time python2.7 ../DNAshapedTFBS.py applyTFFM -T TFFM_trained.xml \
    -i background/test.fa -I background/test.bed \
    -c DNAshapedTFFM_classifier.pkl -o DNAshapedTFFM_bg_predictions.txt \
    -H $helt -E $helt2 -M $mgw -G $mgw2 -P $prot -O $prot2 -R $roll -L $roll2 \
    -n;

echo "Applying the trained PSSM + DNA shape classifier on background sequences.";
time python2.7 ../DNAshapedTFBS.py applyPSSM -f MA0563.1.pfm \
    -i background/test.fa -I background/test.bed \
    -c DNAshapedPSSM_classifier.pkl -o DNAshapedPSSM_bg_predictions.txt \
    -H $helt -E $helt2 -M $mgw -G $mgw2 -P $prot -O $prot2 -R $roll -L $roll2 \
    -n;
