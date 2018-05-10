# AMPEL Preliminary

### Set up sample data

Create samples text file

```
prev=/groups/cbi/Projects/AMPEL_Temple1
for f in $prev/Nick/*1.fq; do echo $(basename $f); done | sed 's/_L1_1.fq//' > samples.txt
```

Create sample directories

```
cat samples.txt | xargs mkdir
````

Symlink sample files

```
cat samples.txt | while read d; do
    ln -s $prev/Nick/${d}_L1_1.fq $d/read_1.fq
    ln -s $prev/Nick/${d}_L1_2.fq $d/read_2.fq
done
```

### Setup references

Create directory for references

```bash
mkdir -p refs
cd refs
```

Link to human reference genome hg38:

```bash
ln -s /lustre/groups/cbi/shared/References/Homo_sapiens/UCSC/hg38
```

Download HERV annotation:

```bash
wget -O HERV_rmsk.hg38.gtf https://github.com/mlbendall/telescope_annotation_db/raw/master/builds/HERV_rmsk.hg38.v2/transcripts.gtf
```

Go back to root directory:

```bash
cd ..
```

### Adapter Sequences

Save adapter sequences for paired-end as `refs/adapters_PE.fa`.

### Run `snakemake`

The snakemake command blocks until everything finishes so you want to run this in tmux.

```
cd /lustre/groups/cbi/Projects/AMPEL_Temple1/telescope
module unload python
module load miniconda3
source activate teleESC

ccmd='sbatch {cluster.args} -N {cluster.N} -p {cluster.p} -t {cluster.t} -o {cluster.out} -J {cluster.name}'
snakemake \
  -j 10000 \
  -c "$ccmd" \
  -u cluster.yaml \
  --cluster-status 'slurm-profile/slurm_status.py' \
  -T -k --ri \
  all
```

### Make count matrix

The following script loads HERV information, alignment metrics, and telescope
reports and creates a count matrix for DESeq2 or other DE software.

```bash
Rscript makeresults.R
```


```
