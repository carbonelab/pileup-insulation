# Pileup insulation

This script is used to pileup the HIC insulation score around a set of features in a bed file. 


## Usage

```
conda env create -f environment.yml -n pileup-insulation
conda activate pileup-insulation
```

```
./pileup_insulation.R -i score.bedgraph -f features -l featureLabal -s hg38 -o outfile_name.pdf
```

