## Pipeline for genome version conversion

Here, we demonstrated how to convert a methylation BED file based on GRCz10 to GRCz11 genome.

1. Download the BED files from NCBI. 

   ```shell
   wget https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM4662075&format=file&file=GSM4662075%5FYueLab%2DWGBS%2Dkidney%2ECG%2Ebw
   ```

2.  Translate the BW file to BedGraph type.

   ```shell
   bigWigToBedGraph GSM4662075_YueLab-WGBS-kidney.CG.bw GSM4662075_YueLab-WGBS-kidney.CG.bed
   ```

3.  Convert the methylation BED file from GRCz10 to GRCz11.

   ```shell
   liftOver GSM4662075_YueLab-WGBS-kidney.CG.bed \
   	danRer10ToDanRer11.over.chain \
   	GSM4662075_YueLab-WGBS-kidney.CG.\
   	GSM4662075_YueLab-WGBS-kidney.GRCz11.bed \
   	unmapped.txt
   ```

   note: the **danRer10ToDanRer11.over.chain** is avialable [here](http://hgdownload.soe.ucsc.edu/downloads.html), the tools including bigWigToBedGraph and liftOver are available [here](http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/).

