# sbx_mbmixture

Lev Litichevskiy  
2022-06-23  

Karl Broman's group published a [pipeline](https://github.com/kbroman/Paper_MBmixups) that identifies microbiome sample mix-ups in metagenomic sequencing datasets. It does this by mapping host reads (which would normally be discarded) to host genotypes. If the host reads are most concordant with the expected host genotype, then we are confident the microbiome sample comes from the expected host. This pipeline can also determine whether a microbiome sample is a mixture of multiple samples. I adapted this pipeline into a [Sunbeam extension](https://sunbeam.readthedocs.io/en/latest/extensions.html) in order to efficiently run it on many samples on HPC.
