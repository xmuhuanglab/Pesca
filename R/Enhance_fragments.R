
enhance_fragment <- function(
  neighbors_pairs, # the neighbor pairs identified in enhancement step
  organs, # organs/tissues data for enhancement
  fragment_path # path storing the original fragment files, whose names should contain the organs' name
){
  neigbor <- neighbors_pairs
  # rename cell name
  colnames(neigbor) <- sapply(colnames(neigbor), function(x){
    unlist(strsplit(x, "[.]"))[2]
  }) %>% unlist()
  
  sccell <- c()
  for ( i in names(neigbor)) {
    cells <- neigbor[[i]]
    sncell <- append(sccell, cells)
  }
  
  fragment <- data.frame()
  enhanced_frag <- data.frame()
  for (i in organ) {
    print(i)
    cells <- grep(i, sccell, value = T)
    print(cells)
    data <- fread(grep(i, list.files("/fragment_path", include.dirs = T, full.names = T), value = T)[1]) %>% as.data.frame()
    data <- subset(data, data$V4 %in% cells)

    for (a in names(neigbor)) {
      if(neigbor[[a]] %in% cells){
        temp <- data[which(data$V4 %in% neigbor[[a]]),]
        temp$V4 <- a
        enhanced_frag <- rbind(enhanced_frag, temp)
      }
    }
    fragment <- rbind(fragment, enhanced_frag)
  }
  write.table(fragment, file = paste0( "Enhanced_frag", ".bed"), row.names = F, col.names = F, sep = "\t", quote = F)
}
