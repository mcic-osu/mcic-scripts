
- Example:

```sh
indir=data/fastq/
outdir=data/fastq_subsample
prop_reads=0.1 # 10% of reads
log=slurm-subsample-fastq-%j.out
sbatch -o $log admin/subsample_fastq_dir.sh -i $indir -o $outdir -p $prop_reads
```
