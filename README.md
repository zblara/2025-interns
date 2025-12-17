To do list: 

1. Run script #1 on raw metagenome reads (fastq) from the Guiuan Dataset: t.ly/b_utO
      Please check if you have the necessary prerequisites via Conda, and list the program versions you used. I will add you as contributors to this git so we can have a repository for all our actions.
      Please edit and adapt any necessary portions in the script such as directory paths, filenames, etc.\
  
3. Edit script #1 (mapping portion) but repalce the mus musculus genome with the human genome. Goal is to remove/filter human DNA contaminants from the reads.
      You can get the human genome here: https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/GCF_000001405.26_GRCh38/GCF_000001405.26_GRCh38_genomic.fna.gz
      Notice that there is a step that requires the "sam file." You will use bowtie2-build to index the human genome and that will be your sam file and reference database for filtering.

5. Run script #3, again please check if you have the prerequisites needed.

3. Script #4: Create a workflow/pipeline that annotates the cleaned/filtered reads for antimicrobial resistance genes (AMR).
	(Hint: https://github.com/zblara/DugesiaAMR)
4. You will use the output of script #3 (taxonomic annotation) and #4 (AMR annotation) as input to make the bipartite graph.
      Feel free to use any tool you would like, such as gephi, igraph, or whatever package.
