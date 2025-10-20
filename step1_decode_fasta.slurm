###generate fasta files
awk -F'\t' 'NR>1 {
    fname = $1 "_IgG.fasta"
    print ">" $1 "_heavy" > fname
    print $2 > fname
    print ">" $1 "_light" >> fname
    print $3 >> fname
}' IgG.tsv 


awk -F'\t' 'NR>1 {
    fname = $1 "_Fab.fasta"
    print ">" $1 "_heavy" > fname
    print $2 > fname
    print ">" $1 "_light" >> fname
    print $3 >> fname
}' Fab.tsv 

###decode scfv in supplemental data
cd /group/ll010/hgao/chromatin_remodel/scFv/histone_train/scfv_WT
awk -F'\t' 'NR>1 {filename=$1 "_scfv.fasta"; print ">" $1 "_scfv" > filename; print $2 ":GGGGGGGGGG" >> filename}' Histone_WT.txt

cd /group/ll010/hgao/chromatin_remodel/scFv/Fab

for i in SMARCC2 SMARCA4 CHD7 CHD4 CHD2 SMARCA2 BRD7
do
cp ${i}_Fab.fasta ../train_Fab_scFv
done
####before run the step2, we have to move each antibody.fasta into individual files!
