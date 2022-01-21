# sbx_mbmixture

Lev Litichevskiy
2022-01-21

Karl Broman's group published a [pipeline](https://github.com/kbroman/Paper_MBmixups) that uses discarded host reads from metagenomic sequencing in order to determine whether a microbiome sample belongs to the host that it should. In other words, the pipeline helps to identify sample mix-ups. It can also determine whether a microbiome sample is a mixture of multiple samples. I adapted this pipeline into a [Sunbeam extension](https://sunbeam.readthedocs.io/en/latest/extensions.html) in order to efficiently run it on many samples on HPC.
