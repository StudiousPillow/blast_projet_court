
To download the SARS-CoV datasets, you can go to the website, https://www.ncbi.nlm.nih.gov/labs/virus/vssi/#/virus?SeqType_s=Nucleotide&VirusLineage_ss=Severe%20acute%20respiratory%20syndrome%20coronavirus%202,%20taxid:2697049&utm_source=data-hub.
Or you can use the command-line tool, downloadable at https://www.ncbi.nlm.nih.gov/datasets/docs/v2/download-and-install/ and use CLI v16+ (API v2alpha) not the old version v13.x (API v1)
You can also use conda :
conda create -n pp-blast
conda activate pp-blast
conda install -c conda-forge ncbi-datasets-cli 

You can then use the download_sequence.sh file.

Ressources : 
Stephen Altschul, Warren Gish, Webb Miller; Eugene Myers; David J. Lipman (1990). "Basic local alignment search tool". Journal of Molecular Biology. 215 (3): 403– 410. doi:10.1016/S0022-2836(05)80360-2. 
Tatiana A. Tatusova, Thomas L. Madden, BLAST 2 Sequences, a new tool for comparing protein and nucleotide sequences, FEMS Microbiology Letters, Volume 174, Issue 2, May 1999, Pages 247–250, https://doi.org/10.1111/j.1574-6968.1999.tb13575.x
Altschul SF, Madden TL, Schäffer AA, Zhang J, Zhang Z, Miller W, Lipman DJ. Gapped BLAST and PSI-BLAST: a new generation of protein database search programs. Nucleic Acids Res. 1997 Sep 1;25(17):3389-402. doi: 10.1093/nar/25.17.3389. PMID: 9254694; PMCID: PMC146917.