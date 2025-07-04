# get fasta and mutate a specific sequence depending on bed file

insilico.mutagenesis = 
  function (reference.region,
            mutation.target.region,
            mutation.character = "-"
    
  ) {
    
  }


species = "Homo sapiens"
assembly = "Hg38"
fastaFolder = "~/fasta_results"
export 	

Foldername.
fileName 	

Filena


myBed <- data.frame(chr=c(1,2),
                    start=c(235265,12356742),
                    end=c(435265,12386742),
                    gene=c("LOC1", "LOC2"))

myFA <- hoardeR::getFastaFromBed(myBed, species="Homo sapiens", fastaFolder="/home/user/fasta/", export=F)