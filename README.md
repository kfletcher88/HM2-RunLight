This is a fork of the original [HaploMerger2](https://github.com/mapleforest/HaploMerger2 "Haplomerger2 - GitHub Repo") software.

The aim of this fork is to provide users with a single line bash script to run the basic Haplomerger2 processes, referred to as hm.batchA, hm.batchB and hm.batchD in the [original documentation](https://github.com/mapleforest/HaploMerger2/blob/master/manual.v3.4.pdf "Haplomerger2 User Manual").

## Citation
I recommend that you cite the original [Haplomerger2 paper](https://www.ncbi.nlm.nih.gov/pubmed/28407147 "Huang et al. 2017") if you use this script.
Depending on the soft-masking strategy deployed (see Inputs below), the masking software should also be cited.

## Inputs
HaploMerger2 runs optimally on a soft-masked assembly. This should be provided to the script with the option `-f`.
Soft-masking may be performed with any software, Haplomerger2 reccomends [WindowMasker](ftp://ftp.ncbi.nih.gov/toolbox/ncbi_tools++/CURRENT/ "WindowMasker: latest version"), described [here](https://www.ncbi.nlm.nih.gov/pubmed/16287941 "WindowMasker citation")
Please refer to the [Haplomerger2 manual](https://github.com/mapleforest/HaploMerger2/blob/master/manual.v3.4.pdf "Haplomerger2 User Manual") for more information and instructions.

The number of threads may also be specified with `-t`.

Finally, by default [HaploMerger2](https://github.com/mapleforest/HaploMerger2 "Haplomerger2 - GitHub Repo") and this bash script will only perform the final stage (remove tanem assemblies) on the reference assembly.
`-a` modifies this and will perform the final stages on the alternative assembly as well as the reference assembly.

## Other Considerations
The original [HaploMerger2](https://github.com/mapleforest/HaploMerger2 "Haplomerger2 - GitHub Repo") software also provides options to run: 

a scaffolding stage with mate-pair reads, using [SSPACE](https://academic.oup.com/bioinformatics/article/27/4/578/197626 "SSPACE citation")
a gap-filling stage with paired-end reads using [GapCloser](https://www.ncbi.nlm.nih.gov/pubmed/23587118 "SOAPdenovo2 citation")

Should users want to make use of these work-flows then the original [HaploMerger2](https://github.com/mapleforest/HaploMerger2 "Haplomerger2 - GitHub Repo") software should be used.

If users want to scaffold or gap-fill with other types of reads (e.g. single molecule), within the workflow, as prescribed by the [HaploMerger2 documentation](https://github.com/mapleforest/HaploMerger2/blob/master/manual.v3.4.pdf "Haplomerger2 User Manual") then the original [HaploMerger2](https://github.com/mapleforest/HaploMerger2 "Haplomerger2 - GitHub Repo") workflow should be modified.
