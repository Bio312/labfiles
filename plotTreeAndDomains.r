#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

#mkdir /home/ec2-user/tools/Rlib4
#wget https://raw.githubusercontent.com/Bio312/lab4files/main/plotTreeAndDomains.r
#Rscript --vanilla plotTreeAndDomains.r hox.bs.mid.suptree hox.rps-blast.out


# test if there is at least one argument: if not, return an error
if (length(args)<2) {
  stop("You must provide a tree and then domain file.", call.=FALSE)
} 

# .libPaths(c("/home/ec2-user/tools/Rlib4","/usr/local/lib64/R/library"))

#install.packages("drawProteins", dependencies = TRUE,repos = "http://cran.us.r-project.org",lib="/home/ec2-user/tools/Rlib4")



if (!require(devtools)) {
    install.packages('devtools')
}
dev_mode(on=TRUE)
install.packages("https://cran.r-project.org/src/contrib/Archive/rlang/rlang_0.4.10.tar.gz", repos = NULL, type="source",lib="/home/ec2-user/R-dev")
devtools::install_github('brennanpincardiff/drawProteins')

library(rlang,lib.loc="/home/ec2-user/R-dev")
library(ggtree)
library(data.table)
#library(drawProteins,lib.loc="/home/ec2-user/tools/Rlib4")
library(drawProteins)
library(ggplot2)



hoxt <- read.tree(args[1]) #"hox.bs.mid.suptree"

t <- ggtree(hoxt, aes(x, y)) + geom_tree() + theme_tree()  + geom_tiplab(cex=0.9) + geom_text2(size=2,aes(subset = !isTip, label=label))
torder <- data.frame("entryName" = get_taxa_name(t), "order" = length(get_taxa_name(t)):1)

hoxd <- fread(args[2]) #"hox.rps-blast.out"

rel_data0 <- hoxd
rel_data0$type <- "DOMAIN"
reldatachain0 <- rel_data0[!(duplicated(rel_data0$V1)),]
reldatachain <- data.frame("V1"="","V2"=reldatachain0$V2,"V3"=1,"V4"=reldatachain0$V2,"V5"=reldatachain0$V1,"type"="CHAIN")

rel_data1 <- rbind(rel_data0,reldatachain)
rel_data2 <- as.data.frame(rel_data1[,c(6,5,3,4,2,1,1)])
rel_data2$taxid <- 1
colnames(rel_data2) <-  c( "type", "description", "begin","end","length","accession","entryName","taxid")
rel_data2$description <- sapply(rel_data2$description, function(dx) substr(dx, 1, 20))
rel_data3 <- merge(rel_data2,torder)

draw_canvas(rel_data3) -> p

p <- draw_chains(p, rel_data3,label_chains = FALSE)


p <- draw_domains(p, rel_data3,label_domains = FALSE) + theme_bw(base_size = 20) + # white background
    theme(panel.grid.minor=element_blank(), 
        panel.grid.major=element_blank()) +
    theme(axis.ticks = element_blank(), 
        axis.text.y = element_blank()) +
    theme(panel.border = element_blank())+
 theme(legend.text=element_text(size=2))+
  theme(legend.key.size = unit(0.15, 'cm'))+
    theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())+
        theme(legend.title= element_blank())+
     theme(plot.margin=unit(c(-5,5,-5,-5), "mm"))+
 theme(legend.position = c(0.1, 0.15))

pdf(paste0(args[1],args[2],".pdf"))
multiplot(t, p,ncol=2)
dev.off()

print(paste0(args[1],args[2],".pdf has been outputted"))


