DNAs2codons <- function(transcripts){
  # x = DNAStringSet
  my_list <- list()
  for(i in 1:length(transcripts)){
    element <- DNAStringSet(transcripts[[i]],
                            start = seq(from = 1, by = 3,  to = length(transcripts[[i]])-2),
                            width = 3)
    element <- as.data.frame(element)
    colnames(element) <- "seq"
    element$codon <- 1:nrow(element)
    
    my_list <- c(my_list, list(element))
  }
  names(my_list) <- names(transcripts)
  return(my_list)
}


free2bind2 <- function(sequence1, sequence2 = "UCCUCC", window = 8, loc_free2bind = "/home/piotr/bin/free2bind"){
  require(Biostrings)
  require(stringr)
  
  sequence2 <- RNAString(sequence2)
  
  dir.create("temp")
  on.exit(unlink("temp", recursive = T))
  
  my_ranges <- data.frame(start = 1:(length(sequence1) - window +1), 
                          end = window:length(sequence1),
                          deltaG = NA,
                          sequence = NA)
  
  
  pb <- txtProgressBar(min = 0, max = nrow(my_ranges), style = 3) 
  for(i in 1:nrow(my_ranges)){
    temp_sequence1 <- sequence1[(my_ranges[i,1]):(my_ranges[i,2])]
    temp_rna <- RNAString(temp_sequence1)
    system(paste("cd ", loc_free2bind, "; ./free_align.pl ", temp_rna, " ", sequence2, " > out.txt", sep = ""))
    a <- readLines(paste(loc_free2bind, "/out.txt", sep = ""))[8]
    deltaG <- as.numeric(str_sub(string = a, start = str_locate(a, "=")[1]+2, end = str_length(a)))
    my_ranges[i,3] <- deltaG
    my_ranges[i,4] <- as.character(temp_rna)
    
    setTxtProgressBar(pb, i)
    
  }
  close(pb)
  unlink("temp", recursive = T)
  return(my_ranges)
}

seq2mfe <- function(sequence, window = 51){
  # as input takes DNA sequence ("DNAString") and calculates secondary structure (mfe in kcal/mol)
  # for all possible sequences of length equal to window
  # returns data.frame object 
  
  require(Biostrings)
  require(stringr)
  
  dir.create("temp")
  on.exit(unlink("temp", recursive = T))
  
  my_ranges <- data.frame(start = 1:(length(sequence) - window +1), 
                          end = window:length(sequence),
                          mfe = NA,
                          mid = NA)
  
  for(i in 1:nrow(my_ranges)){
    my_seq <- sequence[my_ranges[i,1]:my_ranges[i,2]]
    writeXStringSet(RNAStringSet(list(seq1 = my_seq)),
                    filepath = "temp/seqs.fa")
    system("RNAfold < temp/seqs.fa > temp/out.txt")
    my_lines <- readLines("temp/out.txt")[3]
    mfe <- as.numeric(str_sub(my_lines, 
                              start = window+3,
                              end = str_length(my_lines)-1))
    my_ranges[i,3] <- mfe
    my_ranges[i,4] <- ((my_ranges[i,1]):(my_ranges[i,2]))[round(length(((my_ranges[i,1]):(my_ranges[i,2])))/2)]
  }
  unlink("temp", recursive = T)
  unlink(x = c("seq1_ss.ps"))
  
  return(my_ranges)
}

mutate.codons <- function(seq, codon_mut, codon_ext, max_comb = 1000, max_res = 10, codons_file = "data/codons"){
  require(dplyr)
  
  
  codons_table <- read.table(codons_file, header = T)
  my_codons <- DNAs2codons(DNAStringSet(seq))[[1]]
  
  ## makink a list with all possible codons at each position to be mutated
  selected_codons <- my_codons[codon_mut,]
  alt_codons <- vector("list", length = nrow(selected_codons))
  for(i in 1:nrow(selected_codons)){
    my_aa <- codons_table[codons_table$codon == selected_codons[i,1],2]
    alt_codons[[i]] <- as.vector(codons_table[codons_table$aa == my_aa,1])
  }
  names(alt_codons) <- codon_mut
  message(paste("There are ", Reduce("*", sapply(alt_codons, length)), " sequence combinations.", sep = ""))
  
  
  # 1000000 is a good limit, if the value is higher expand.grid function is very slow.
  if(Reduce("*", sapply(alt_codons, length)) > 5000000){
    repeat{
      selected_codon <- as.character(sample(codon_mut, 1)) # Random selection !!
      alt_codons[[selected_codon]] <- alt_codons[[selected_codon]][1]
      if(Reduce("*", sapply(alt_codons, length)) < 5000000) break
    }
  }  
  
  # creating all aviable codon combinations
  codons_comb <- expand.grid(alt_codons)
  
  if(Reduce("*", sapply(alt_codons, length)) > max_comb){
    message(paste("Selecting", max_comb, "random sequence combinations."))
    #set.seed(1) # ??
    codons_comb <- codons_comb[sample(1:Reduce("*", sapply(alt_codons, length)), size = max_comb),]
    codons_comb_STR <- apply(codons_comb, 1, paste, collapse = "")
  }else{
    message("Selecting all combinations.")
    codons_comb_STR <- apply(codons_comb, 1, paste, collapse = "")
  }
  seq_variants <- DNAStringSet(codons_comb_STR)
  
  message("Predicting structure of mutated regions")
  
  my_list <- vector("list", length = length(seq_variants))
  pb <- txtProgressBar(min = 0, max = length(seq_variants), style = 3)
  for(i in 1:length(seq_variants)){
    my_temp <- seq2mfe(sequence = seq_variants[[i]])
    my_temp <- my_temp[,4:3]
    my_temp$name <- names(seq_variants)[i]
    my_list[[i]] <- my_temp
    setTxtProgressBar(pb, i)
  }
  close(pb)
  my_data <- do.call(rbind.data.frame, my_list)
  
  selected_mutations        <- my_data %>% arrange(-mfe) %>% head(max_res)
  selected_mutations_seqs   <- seq_variants[names(seq_variants) %in% selected_mutations$name]
  
  left_b <- paste(my_codons[(codon_mut[1]-codon_ext):(codon_mut[1]-1),1], collapse = "")
  right_b <- paste(my_codons[(last(codon_mut)+1):(last(codon_mut)+codon_ext),1], collapse = "")
  
  selected_expanded_seqs <- DNAStringSet(paste(left_b, as.character(selected_mutations_seqs), right_b, sep = ""))
  names(selected_expanded_seqs) <- names(selected_mutations_seqs)
  
  message("Predicting structure for expanded regions")
  
  my_list <- vector("list", length = length(selected_expanded_seqs)+1)
  pb <- txtProgressBar(min = 0, max = length(selected_expanded_seqs), style = 3)
  for(i in 1:length(selected_expanded_seqs)){
    my_temp <- seq2mfe(sequence = selected_expanded_seqs[[i]])
    my_temp <- my_temp[,4:3]
    my_temp$name <- names(selected_expanded_seqs)[i]
    my_temp$type <- 0.5
    my_list[[i]] <- my_temp
    setTxtProgressBar(pb, i)
  }
  
  WT_sequence <- DNAString(paste(my_codons[(codon_mut[1]-codon_ext):(last(codon_mut)+codon_ext),1], collapse = ""))
  my_temp <- seq2mfe(WT_sequence)
  my_temp <- my_temp[,4:3]
  my_temp$name <- "WT"
  my_temp$type <- 1
  my_list[[i+1]] <- my_temp
  
  close(pb)
  my_data <- do.call(rbind.data.frame, my_list)
  
  return(list(numeric_data = my_data, sequences = selected_mutations_seqs))
} 



mutate.codons2 <- function(seq, codon_mut, codon_ext, max_comb = 1000, max_res = 10, codons_file = "data/codons"){
  require(dplyr)
  
  
  codons_table <- read.table(codons_file, header = T)
  my_codons <- DNAs2codons(DNAStringSet(seq))[[1]]
  
  ## makink a list with all possible codons at each position to be mutated
  selected_codons <- my_codons[codon_mut,]
  alt_codons <- vector("list", length = nrow(selected_codons))
  for(i in 1:nrow(selected_codons)){
    my_aa <- codons_table[codons_table$codon == selected_codons[i,1],2]
    alt_codons[[i]] <- as.vector(codons_table[codons_table$aa == my_aa,1])
  }
  names(alt_codons) <- codon_mut
  message(paste("There are ", Reduce("*", sapply(alt_codons, length)), " sequence combinations.", sep = ""))
  
  
  # 1000000 is a good limit, if the value is higher expand.grid function is very slow.
  if(Reduce("*", sapply(alt_codons, length)) > 1000000){
    repeat{
      selected_codon <- as.character(sample(codon_mut, 1)) # Random selection !!
      alt_codons[[selected_codon]] <- alt_codons[[selected_codon]][1]
      if(Reduce("*", sapply(alt_codons, length)) < 1000000) break
    }
  }  
  
  # creating all aviable codon combinations
  codons_comb <- expand.grid(alt_codons)
  
  if(Reduce("*", sapply(alt_codons, length)) > max_comb){
    message(paste("Selecting", max_comb, "random sequence combinations."))
    #set.seed(1) # ??
    codons_comb <- codons_comb[sample(1:Reduce("*", sapply(alt_codons, length)), size = max_comb),]
    codons_comb_STR <- apply(codons_comb, 1, paste, collapse = "")
  }else{
    message("Selecting all combinations.")
    codons_comb_STR <- apply(codons_comb, 1, paste, collapse = "")
  }
  seq_variants <- DNAStringSet(codons_comb_STR)
  
  message("Predicting structure of mutated regions")
  
  my_list <- vector("list", length = length(seq_variants))
  pb <- txtProgressBar(min = 0, max = length(seq_variants), style = 3)
  for(i in 1:length(seq_variants)){
    my_temp <- seq2mfe(sequence = seq_variants[[i]])
    my_temp <- my_temp[,4:3]
    my_temp$name <- names(seq_variants)[i]
    my_list[[i]] <- my_temp
    setTxtProgressBar(pb, i)
  }
  close(pb)
  my_data <- do.call(rbind.data.frame, my_list)
  
  selected_mutations        <- my_data %>% arrange(mfe) %>% head(max_res)
  selected_mutations_seqs   <- seq_variants[names(seq_variants) %in% selected_mutations$name]
  
  left_b <- paste(my_codons[(codon_mut[1]-codon_ext):(codon_mut[1]-1),1], collapse = "")
  right_b <- paste(my_codons[(last(codon_mut)+1):(last(codon_mut)+codon_ext),1], collapse = "")
  
  selected_expanded_seqs <- DNAStringSet(paste(left_b, as.character(selected_mutations_seqs), right_b, sep = ""))
  names(selected_expanded_seqs) <- names(selected_mutations_seqs)
  
  message("Predicting structure for expanded regions")
  
  my_list <- vector("list", length = length(selected_expanded_seqs)+1)
  pb <- txtProgressBar(min = 0, max = length(selected_expanded_seqs), style = 3)
  for(i in 1:length(selected_expanded_seqs)){
    my_temp <- seq2mfe(sequence = selected_expanded_seqs[[i]])
    my_temp <- my_temp[,4:3]
    my_temp$name <- names(selected_expanded_seqs)[i]
    my_temp$type <- 0.5
    my_list[[i]] <- my_temp
    setTxtProgressBar(pb, i)
  }
  
  WT_sequence <- DNAString(paste(my_codons[(codon_mut[1]-codon_ext):(last(codon_mut)+codon_ext),1], collapse = ""))
  my_temp <- seq2mfe(WT_sequence)
  my_temp <- my_temp[,4:3]
  my_temp$name <- "WT"
  my_temp$type <- 1
  my_list[[i+1]] <- my_temp
  
  close(pb)
  my_data <- do.call(rbind.data.frame, my_list)
  
  return(list(numeric_data = my_data, sequences = selected_mutations_seqs))
} 