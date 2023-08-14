
################ load data
load("./Enhancement/Step3_3.Neighbors_Pairs.rds")

colnames(neigbor)=sapply(colnames(neigbor), function(x){
  unlist(strsplit(x,"[.]"))[2]
}) %>% unlist()

sccell=c()
for ( i in names(neigbor)) {
  cells = neigbor[[i]]
  sncell=append(sccell,cells)
}

organ=c("eye","forebrain","midbrain","hindbrain","limb")
fragment=data.frame()
enhanced_frag=data.frame()
for (i in organ) {
  print(i)
  cells=grep(i,sccell,value = T)
  print(cells)
  data=fread(grep(i,list.files("/fragment_path",include.dirs = T,full.names =T),value = T)[1]) %>% as.data.frame()
  data=subset(data,data$V4 %in% cells)
  
  for (a in names(neigbor)) {
    if(neigbor[[a]] %in% cells){
      temp=data[which(data$V4 %in% neigbor[[a]]),]
      temp$V4=a
      enhanced_frag=rbind(enhanced_frag,temp)
    }
  }
  fragment=rbind(fragment,enhanced_frag)
}
write.table(fragment,file =paste0( "D1.enhanced_frag",".bed"),row.names = F,col.names = F,sep = "\t",quote = F)

######### linux
# cd PATH/bcftools-1.2/htslib-1.2.1/
# export PATH=PATH/bcftools-1.2/htslib-1.2.1:$PATH
# 
# cd PATH/tabix-0.2.6/tabix/tabix
# export PATH/tabix-0.2.6/tabix/tabix:$PATH
# sort -V PATH/D1_enhanced_frag.bed -o D1_enhanced_frag_sorted.bed
# bgzip PATH/D1_enhanced_frag_sorted.bed
# tabix -p bed PATH/D1_enhanced_frag_sorted.bed.gz

