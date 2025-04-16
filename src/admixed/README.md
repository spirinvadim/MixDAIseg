# DAIseg.mex
Highly accurate method for detecting archaic segments in the modern admixed genomes 


![Demography](https://github.com/Genomics-HSE/DAIseg.mex/blob/main/Mex.svg)


DAIseg.mex method is created to detect ancient introgressed segments using unadmixed outgroup population and several reference archaic genomes. 





The description of the used files is [here][1]




# Preparations of files

Prepare file with ancestral alleles, use guide in the directory./preparations/hg19(grch38) OR download prepared files from [here][6]   



To read more details for files preparation see [readme][2]. To avoid details use script 
```bash
 ./full.preparation.Linux(MacOS).sh hg19 22 path.to/file.bed n1 n2 n3 1000GP path.to/mex  path.to/europeans path.to/americans path.to/africans path.to/archaic name.out_vcf name.out_txt path.to/ancestral.allele.file
```
where 1000GP, n1, n2, n3 are [vcf with 1000 Genome project][3], [neanderthal vcfs][4] and [chagyrskaya][5] (needed to be splited). 

Write full path to the list of samples path.to/outgroup.list,  path.to/obserables.list. path.to/archaic.list

name.out_vcf name.out_txt are the names of the resulting files(be saved in the hg19 directory).





There are two options without EM-algorithm and with EM algorithm. 


# Run DAI.seg without EM algorithm



```bash
python3 daiseg.mex.2.py --obs_samples path.to/obserables.list --bed path.to/file.bed   --HMM_par par.file.txt --EM no --prepared_file ./hg19/name.out_txt --out_prefix out.chr --arch_cover ./hg19/arch.covering.chr22.txt --transition_matrix full --obs_type simple
```


# Run DAI.seg using EM algorithm

```bash
python3 daiseg.mex.2.py --obs_samples path.to/obserables.list --bed path.to/file.bed   --HMM_par par.file.txt --EM yes --EM_steps 20 --EM_samples 10 --prepared_file ./hg19/allels.ref.and.obs.chr22.txt --out_prefix out.chr22 --arch_cover ./hg19/arch.covering.chr22.txt --transition_matrix full --obs_type simple
```
to obtain estimations of the  coalescent times and run DAIseg. Here par.file.txt is used as the initial guess for EM algorithm.
--EM_samples is the number of sampes for training.




[1]: [here](https://github.com/Genomics-HSE/DAIseg/blob/main/File.types.md)
[2]: https://github.com/Genomics-HSE/DAIseg/blob/main/hg19/README.md
[6]: https://drive.google.com/drive/folders/1_zE9eaV3psFPRdFatkq-R1yGluvjgiX6?usp=sharing

[3]: http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chr22.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz 
[4]: http://cdna.eva.mpg.de/neandertal/Vindija/VCF/
[5]: http://ftp.eva.mpg.de/neandertal/ChagyrskayaOkladnikov/
[6]: https://drive.google.com/drive/folders/1_zE9eaV3psFPRdFatkq-R1yGluvjgiX6?usp=drive_link
