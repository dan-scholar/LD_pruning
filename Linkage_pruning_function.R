LD_elimination_iterations<- function(dataframe, file_path, rep_LD_matrix, r2) {
  d=dataframe
  file_path = file_path
  g=dataframe
  rep_LD_matrix=rep_LD_matrix
  j=0
  r=r2
  dchunk_P <- d[ order(d[,"chunk"],d[,"PValue"]), ]
  leading_SNPs <- dchunk_P[!duplicated(dchunk_P$chunk), ][,"SNP"]
  repeat { 
    j= j + 1
    assign(paste("d_IT", j, sep = ""), subset (d, R2 < r)[,c("SNP","PValue","SnpType", "Chr", "Position", "chunk")])
    d$chunk <-factor(d$chunk)
    assign(paste("d_IT", j,"_chunk_P", sep = ""),get(paste("d_IT", j, sep = ""))[ order(get(paste("d_IT", j, sep = ""))[,"chunk"],get(paste("d_IT", j, sep = ""))[,"PValue"]), ])
    assign(paste("d_IT", j,"_chunk_P_unique", sep = ""),get(paste("d_IT", j,"_chunk_P", sep = ""))[!duplicated(get(paste("d_IT", j,"_chunk_P", sep = ""))$chunk), ][,"SNP"])
    l=1
    z=NULL
    for(i in 1:nrow(get(paste("d_IT", j, sep = "")))) { 
      if(as.numeric(as.factor(get(paste("d_IT", j, sep = ""))$chunk[i])) == l )  { 
        z[i]<- c(rep_LD_matrix[as.character(as.factor(get(paste("d_IT", j,"_chunk_P_unique", sep = ""))[l])),as.character(as.factor(get(paste("d_IT", j, sep = ""))[i,"SNP"]))])
      } else { 
        l <- l+ 1
        z[i] <- c(rep_LD_matrix[as.character(as.factor(get(paste("d_IT", j,"_chunk_P_unique", sep = ""))[l])),as.character(as.factor(get(paste("d_IT", j, sep = ""))[i,"SNP"]))])
      } 
    }
    assign(paste("d_IT", j, sep = ""), `[[<-`(get(paste("d_IT", j, sep = "")), 'R2', value = z))
    d <- get(paste("d_IT", j, sep = ""))
    if(j==1){
      assign(paste("Loop_list_of_rep_SNP_R2_", j,"_iterations", sep = ""),as.character(get(paste("d_IT", j,"_chunk_P_unique",sep = "")))) 
    }
    else{
      assign(paste("Loop_list_of_rep_SNP_R2_", j,"_iterations" ,sep = ""), c(as.character(get(paste("d_IT", j,"_chunk_P_unique",sep = ""))), get(paste("Loop_list_of_rep_SNP_R2_", j-1 ,"_iterations" ,sep = ""))) )
    }
    assign(paste("test", j, sep = ""), subset (d, R2 < r)[,c("SNP","PValue","SnpType", "Chr", "Position", "chunk")])
    if(nrow(get(paste("test", j, sep = ""))) ==0) break }
  j= j + 1
  assign(paste("Loop_list_of_rep_SNP_R2_", j,"_iterations" ,sep = ""), c(as.character(leading_SNPs), as.character(d$SNP), get(paste("Loop_list_of_rep_SNP_R2_", j-1 ,"_iterations" ,sep = ""))) )
  print(length(unique(get(paste("Loop_list_of_rep_SNP_R2_", j,"_iterations" ,sep = "")))) )
  write.table(g[g$SNP %in%  as.character(unique(get(paste("Loop_list_of_rep_SNP_R2_", j,"_iterations" ,sep = "")))),], file_path, quote=F, row.names=FALSE, sep="\t")
  #
  list(number.iterations=j)
}