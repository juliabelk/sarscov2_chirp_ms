rfmt <- function(x) {
    return(unlist(lapply(x,function(i){
        tmp <- strsplit(i,"\\.")[[1]]
        return(paste(tmp[3:length(tmp)],sep = "",collapse = "."))
    })))
}

parse_uniprot <- function(x){
    tmp <- strsplit(strsplit(x,";")[[1]][1],"\\|")[[1]][3]
    
    if (endsWith(tmp,"HUMAN") | endsWith(tmp,"CHLSB")) return(tmp)
    
    if (strsplit(strsplit(x,";")[[1]][1],"\\|")[[1]][2] == "SARSCOV2") {
        return(paste0("SARS-",strsplit(tmp,"-")[[1]][2]))
    }
    return(paste0(strsplit(tmp,"_")[[1]][1],strsplit(tmp,"_")[[1]][2]))
}

parse_species <- function(x){
    tmp <- strsplit(strsplit(x,";")[[1]][1],"\\|")[[1]][3]

    if (endsWith(tmp,"HUMAN")) return("HUMAN")
    if (endsWith(tmp,"CHLSB")) return("CHLSB")

    return(strsplit(strsplit(x,";")[[1]][1],"\\|")[[1]][2])
}

make_bar_plot <- function(df, sw="HNRNP",outpth=NA,genelist=NA,fixord=F,barout=NA) {
    
    options(repr.plot.width=24, repr.plot.height=10) 
    if (is.na(genelist)) genelist <- df$name[which(startsWith(df$name,sw))]
    
    new <- data.frame("gene"=c(),"virus"=c(),"enrich"=c())
    
    a <- c(
        "huhcov.d1.enrich.mean","huhcov.d2.enrich.mean","vero.d1.enrich.mean","vero.d2.enrich.mean"#,
        #"huhcov.d2.enrich.mean","vero.d2.enrich.mean",
        #"zika.enrich.mean","deng.enrich.mean",'rv.enrich.1'
    )#,"rv.enrich") "huhcov.d1.enrich",
    
    for (g in genelist) {
        for (i in a) {
            if (!(g %in% df$name)) {
                new <- rbind(new,data.frame("gene"=g,"virus"=i,"enrich"=NA))
                next
            }
            new <- rbind(new,
                data.frame(
                    "gene"=g,
                    "virus"=i,
                    "enrich"=df[which(df$name==g),i]
                )
            )
        }
    }
    
    new[which(is.na(new$enrich)),"enrich"] <- 0
    new[which(new$enrich < 0),"enrich"] <- 0
    
    if (!fixord) {
        new$gsort <- unlist(lapply(new$gene,function(x) 
            sum(new[which(new$gene == x),"enrich"]) 
        ))
        new <- new[order(new$gsort,decreasing = T),]# %>% arrange(gsort)

        new$gene <- factor(new$gene, levels=unique(new$gene))#c(sort(as.character(unique(new$gene)))))
    }
    else new$gene <- factor(new$gene, levels=genelist)
        #new$gene <- as.character(new$gene)
    
    plt <- ggplot(new, aes(x=gene, y=enrich, fill=virus)) + 
      theme_classic() +
      theme(axis.text.x = element_text(angle = 270, hjust = 0, vjust = 0.5)) +
      theme(panel.grid.minor = element_blank()) +
      theme(text = element_text(size=24,face="bold")) +
      theme(axis.ticks = element_line(colour = "black", size = 2)) +
      #theme(panel.border = element_rect(size=1.5)) +
      scale_fill_brewer(palette="Paired") +
      geom_bar(stat="identity", position=position_dodge(), color="black")
    
    if (min(new$enrich,na.rm=T) >=0) {
        plt <- plt + scale_y_continuous(expand = expansion(mult = c(0, .05)))
    }
    if (!is.na(barout)) {
        ggsave(paste0(barout,sw,"-allsars.eps"),plot = plt,width = 16, height = 8, dpi = 300)
    }
    return(plt)
}

make_gene_set_heatmap <- function(df, sw, colsel="everything", glist=c(),
                                  fixord=F,rmempty=T,returnMat=F,clip=100000,repna=F,setmax=NA,outs=NA) {

    options(repr.plot.width=32, repr.plot.height=12)
    a <- c(
        'huhcov.d1.enrich.1','huhcov.d1.enrich.2','huhcov.d1.enrich.3' ,
        'huhcov.d2.enrich.1', 'huhcov.d2.enrich.2', 'huhcov.d2.enrich.3' ,
        'zika.enrich.1', 'zika.enrich.2', 'zika.enrich.3',
        'deng.enrich.1', 'deng.enrich.2', 'deng.enrich.3', 
        'rv.enrich.1'
    )
    b <- c(
        'huhcov.d1.enrich.1','huhcov.d1.enrich.2','huhcov.d1.enrich.3' ,
        'huhcov.d2.enrich.1', 'huhcov.d2.enrich.2', 'huhcov.d2.enrich.3' ,
        'vero.d1.enrich.1','vero.d1.enrich.2','vero.d1.enrich.3' ,
        'vero.d2.enrich.1', 'vero.d2.enrich.2', 'vero.d2.enrich.3' 
    )
    d <- c(
        'huhcov.d1.enrich.1','huhcov.d1.enrich.2','huhcov.d1.enrich.3' ,
        'huhcov.d2.enrich.1', 'huhcov.d2.enrich.2', 'huhcov.d2.enrich.3',
        'vero.d1.enrich.1','vero.d1.enrich.2','vero.d1.enrich.3' ,
        'vero.d2.enrich.1', 'vero.d2.enrich.2', 'vero.d2.enrich.3',
        'zika.enrich.1', 'zika.enrich.2', 'zika.enrich.3',
        'deng.enrich.1', 'deng.enrich.2', 'deng.enrich.3', 
        'rv.enrich.1'
    )
    ave <- c(
        'huhcov.d1.enrich.mean',
        'huhcov.d2.enrich.mean',
        'zika.enrich.mean',
        'deng.enrich.mean',
        'rv.enrich.1'
    )
    allave <- c(
        'huhcov.d1.enrich.mean',
        'huhcov.d2.enrich.mean',
        'vero.d1.enrich.mean',
        'vero.d2.enrich.mean',
        'zika.enrich.mean',
        'deng.enrich.mean',
        'rv.enrich.1'
    )
    if (colsel == "all_human") cname <- a
    else if (colsel == "all_sars") cname <- b
    else if (colsel == "averages") cname <- ave
    else if (colsel == "allaverages") cname <- allave
    else cname <- d

    if (length(glist) >= 1) {
        tmp <- match(glist,df$name)
        tmp <- tmp[which(!is.na(tmp))] #glist #which(df$name %in% glist)
    } else {
        tmp <- which(startsWith(df$name,sw))
    }
    #print(tmp)
    mat <- as.matrix(df[tmp,cname])
    rownames(mat) <- df$name[tmp]

    #clp <- 2
    mat[mat < 0] <- 0 
    mat[mat > clip] <- clip 
    if (repna) mat[is.na(mat)] <- 0
    
    if (!fixord) mat <- mat[order(rowSums(mat,na.rm=T),decreasing=T),]
    
    if (rmempty) mat <- mat[rowSums(is.na(mat)) != ncol(mat), ]
    
    mat <- t(mat)
    
    if (returnMat) return(mat)
    
    mat_breaks <- seq(0, max(mat,na.rm = T), length.out = 10)
    if (!is.na(setmax)) mat_breaks <-  seq(0, setmax, length.out = 10)
    
    clust = !fixord
    if (length(which(is.na(mat))) > 0) clust=FALSE
    
    plt <- pheatmap(
        mat,cellwidth=35, cellheight=35, border_color=NA, fontsize=20,
        cluster_rows = F,cluster_cols = clust,na_col = "grey70",
        color = viridis(length(mat_breaks) - 1),#,direction = -1),
        breaks = mat_breaks
    )
    
    if (!is.na(outs)) {
      pdf(paste0(outs,colsel,"-",sw,".pdf"),width = 36,height=12)
      print(plt)
      dev.off()
    }
    return(plt)
}



