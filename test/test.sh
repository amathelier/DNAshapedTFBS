export PYTHONPATH=$PYTHONPATH:[MODIFY HERE IF NECESSARY]
# Define where are located the DNA shape data from GBshape
araTha=[MODIFY HERE]
helt=$araTha/araTha10.HelT.bigWig;
mgw=$araTha/araTha10.MGW.bigWig;
prot=$araTha/araTha10.ProT.bigWig;
roll=$araTha/araTha10.Roll.bigWig;
helt2=$araTha/araTha10.HelT.2nd.wig.bw;
mgw2=$araTha/araTha10.MGW.2nd.wig.bw;
prot2=$araTha/araTha10.ProT.2nd.wig.bw;
roll2=$araTha/araTha10.Roll.2nd.wig.bw;

echo "Training a first order TFFM + DNA shape classifier.";
time python2.7 ../DNAshapedTFBS.py trainTFFM -T TFFM_trained.xml \
    -i foreground/train.fa -I foreground/train.bed \
    -b background/train.fa -B background/train.bed \
    -o DNAshapedTFFM_fo_classifier -t first_order \
    -1 $helt $mgw $prot $roll -2 $helt2 $mgw2 $prot2 $roll2 -n;

echo "Training a detailed TFFM + DNA shape classifier.";
time python2.7 ../DNAshapedTFBS.py trainTFFM -T TFFM_trained.xml \
    -i foreground/train.fa -I foreground/train.bed \
    -b background/train.fa -B background/train.bed \
    -o DNAshapedTFFM_d_classifier -t detailed \
    -1 $helt $mgw $prot $roll -2 $helt2 $mgw2 $prot2 $roll2 -n;

echo "Training a PSSM + DNA shape classifier.";
time python2.7 ../DNAshapedTFBS.py trainPSSM -f MA0563.1.pfm \
    -i foreground/train.fa -I foreground/train.bed \
    -b background/train.fa -B background/train.bed \
    -o DNAshapedPSSM_classifier \
    -1 $helt $mgw $prot $roll -2 $helt2 $mgw2 $prot2 $roll2 -n;

echo "Training a 4-bits + DNA shape classifier.";
time python2.7 ../DNAshapedTFBS.py train4bits -f MA0563.1.pfm \
    -i foreground/train.fa -I foreground/train.bed \
    -b background/train.fa -B background/train.bed \
    -o DNAshaped4bits_classifier \
    -1 $helt $mgw $prot $roll -2 $helt2 $mgw2 $prot2 $roll2 -n;

echo "Applying the trained first order TFFM + DNA shape classifier on foreground sequences.";
time python2.7 ../DNAshapedTFBS.py applyTFFM -T TFFM_trained.xml \
    -i foreground/test.fa -I foreground/test.bed \
    -c DNAshapedTFFM_fo_classifier.pkl -o DNAshapedTFFM_fo_fg_predictions.txt \
    -1 $helt $mgw $prot $roll -2 $helt2 $mgw2 $prot2 $roll2 -n -t first_order;

echo "Applying the trained detailed TFFM + DNA shape classifier on foreground sequences.";
time python2.7 ../DNAshapedTFBS.py applyTFFM -T TFFM_trained.xml \
    -i foreground/test.fa -I foreground/test.bed \
    -c DNAshapedTFFM_d_classifier.pkl -o DNAshapedTFFM_d_fg_predictions.txt \
    -1 $helt $mgw $prot $roll -2 $helt2 $mgw2 $prot2 $roll2 -n -t detailed;

echo "Applying the trained PSSM + DNA shape classifier on foreground sequences.";
time python2.7 ../DNAshapedTFBS.py applyPSSM -f MA0563.1.pfm \
    -i foreground/test.fa -I foreground/test.bed \
    -c DNAshapedPSSM_classifier.pkl -o DNAshapedPSSM_fg_predictions.txt \
    -1 $helt $mgw $prot $roll -2 $helt2 $mgw2 $prot2 $roll2 -n;

echo "Applying the trained 4-bits + DNA shape classifier on foreground sequences.";
time python2.7 ../DNAshapedTFBS.py apply4bits -f MA0563.1.pfm \
    -i foreground/test.fa -I foreground/test.bed \
    -c DNAshaped4bits_classifier.pkl -o DNAshaped4bits_fg_predictions.txt \
    -1 $helt $mgw $prot $roll -2 $helt2 $mgw2 $prot2 $roll2 -n;

echo "Applying the trained first order TFFM + DNA shape classifier on background sequences.";
time python2.7 ../DNAshapedTFBS.py applyTFFM -T TFFM_trained.xml \
    -i background/test.fa -I background/test.bed \
    -c DNAshapedTFFM_fo_classifier.pkl -o DNAshapedTFFM_fo_bg_predictions.txt \
    -1 $helt $mgw $prot $roll -2 $helt2 $mgw2 $prot2 $roll2 -n -t first_order;

echo "Applying the trained detailed TFFM + DNA shape classifier on background sequences.";
time python2.7 ../DNAshapedTFBS.py applyTFFM -T TFFM_trained.xml \
    -i background/test.fa -I background/test.bed \
    -c DNAshapedTFFM_d_classifier.pkl -o DNAshapedTFFM_d_bg_predictions.txt \
    -1 $helt $mgw $prot $roll -2 $helt2 $mgw2 $prot2 $roll2 -n -t detailed;

echo "Applying the trained PSSM + DNA shape classifier on background sequences.";
time python2.7 ../DNAshapedTFBS.py applyPSSM -f MA0563.1.pfm \
    -i background/test.fa -I background/test.bed \
    -c DNAshapedPSSM_classifier.pkl -o DNAshapedPSSM_bg_predictions.txt \
    -1 $helt $mgw $prot $roll -2 $helt2 $mgw2 $prot2 $roll2 -n;

echo "Applying the trained 4-bits + DNA shape classifier on background sequences.";
time python2.7 ../DNAshapedTFBS.py apply4bits -f MA0563.1.pfm \
    -i background/test.fa -I background/test.bed \
    -c DNAshaped4bits_classifier.pkl -o DNAshaped4bits_bg_predictions.txt \
    -1 $helt $mgw $prot $roll -2 $helt2 $mgw2 $prot2 $roll2 -n;
