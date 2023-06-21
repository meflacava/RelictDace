#### Split barcodes on FARM ####

#This script will create a list of commands to demultiplex individuals in every combination of lane and plate
# - Relict dace had 15 plates each sequenced on 3 lanes, so this script creates 45 files, each with up to 96 commands listed

#Import barcode info: SampleID, Plate, Well, WellBarcode, PlateIndex, IndexCode
file <- read.table("/Users/melanielacava/Library/CloudStorage/Box-Box/Postdoc/RelictDace/Bioinformatics/IndBarcodes_RD.txt",header=T,stringsAsFactors=F)
head(file)

#Separate command list for each plate and each lane
for (i in c(1,2,8)){ #lane
  for (j in unique(file$Plate)){ #plate
    x <- matrix(data=NA,nrow=length(unique(file$SampleID[file$Plate==j])),ncol=5)
    
    #fill consistent columns
    x[,1] <- "/home/mlacava/scripts/BarcodeSplitListBestRadPairedEnd.pl"
    x[,2] <- paste0("/home/mlacava/raw_data/RelictDace_BMAG070/demulti",i,"/plate",formatC(j,width=2,flag="0"),"_R1.fastq")
    x[,3] <- paste0("/home/mlacava/raw_data/RelictDace_BMAG070/demulti",i,"/plate",formatC(j,width=2,flag="0"),"_R2.fastq")
    count <- 1
    
    for (k in unique(file$SampleID[file$Plate==j])){ #individuals
      x[count,4] <- paste0("GG",file$WellBarcode[file$SampleID==k]) #add GG to beginning of each barcode
      x[count,5] <- paste0("/home/mlacava/raw_data/RelictDace_BMAG070/demulti",i,"/inds/",k)
      count <- count + 1
    }
    
    #save command list for each lane/plate combo
    write.table(x,paste0("/Users/melanielacava/Library/CloudStorage/Box-Box/Postdoc/RelictDace/Bioinformatics/split_ind/lane",i,"_plate",formatC(j,width=2,flag="0"),".sh"),quote=F,row.names=F,col.names=F)
  }
}