
pca_plt <- function(df,colsel,nm,group1="sample",group2=NA,offset=9,rmneg=F,corrout=NA) {
    
    options(repr.plot.width=10, repr.plot.height=8)   
    
    x <- as.matrix(df[,colsel])
    
    if (rmneg) {
        x[x < 0] <- 0
        x <- x[which(rowSums(x) > 0),]
    }
    pca <- prcomp(t(x),scale=TRUE)
    var_explained <- pca$sdev^2/sum(pca$sdev^2)

    d <- as.data.frame(pca$x)
    d$name  <- rownames(d)
    d$item1 <- unlist(lapply(d$name,function(x) strsplit(x,"_")[[1]][1] ))
    d$item2 <- unlist(lapply(d$name,function(x) strsplit(x,"_")[[1]][2] ))
    d$sample <- substr(d$name,1,nchar(d$name)-offset)
    d$virus <- unlist(lapply(d$sample,function(x){
        return(strsplit(x,"\\.")[[1]][1])
    }))
    xpc <- 1
    ypc <- 2

    plt <- ggplot(d,aes_string(x=paste0("PC",xpc),y=paste0("PC",ypc))) 

    if (is.na(group2)) plt <- plt + geom_point(aes_string(color=group1),size=6)
    else plt <- plt + geom_point(aes_string(color=group1,shape=group2),size=6)

    plt <- plt +
        theme_bw(base_size=24) + 
        labs(x=paste0("PC",xpc,": ",round(var_explained[xpc]*100,1),"%"),
             y=paste0("PC",ypc,": ",round(var_explained[ypc]*100,1),"%")) +
        theme(panel.grid.minor = element_blank()) +
        theme(panel.grid.major = element_blank()) +
        theme(aspect.ratio=1) +
        theme(panel.border = element_rect(size=1.5))

    if (!is.na(corrout)) {
      ggsave(paste0(corrout,nm,"-pca.eps"),plot = plt,width = 8, height = 8, dpi = 300)
    }
    return(plt)
}

make_rugplot <- function(f1, glist, name, fixord=F, out=NA) {
    
    if (length(glist) == 0) return
    #print(glist)

    options(repr.plot.width=6, repr.plot.height=2.5+0.5*length(glist))

    if (!fixord) {
        f1 <- f1 %>% arrange(meanres)
        glist <- unique(f1[sort(which(f1$gene %in% glist)),"gene"])
        #print(glist)
    }
    print(glist)
    
    height <- 0.4
    spacing <- 0.5
    offset <- height / 2 - 0.1

    #xlim <- max(abs(

    plt <- ggplot(f1) +
      geom_density(aes(x = residual,y=..scaled..)) +
      #scale_x_continuous(limits=c(-3,3)) +
      theme_bw(base_size=12) +
      theme(text = element_text(size=20)) +
      theme(panel.grid.minor = element_blank()) +
      ylab(NULL) +
      ggtitle(name) +
      #scale_y_continuous(expand = c(0, 0)) +
      theme(panel.border = element_rect(colour = "black", fill=NA, size=2))+
      scale_y_continuous(
          breaks = c(0,seq(offset - spacing, offset - spacing*(length(glist)), by = -spacing)),
          labels=c("",glist),expand=c(0,0.5)
      )

    xlim <- 3
    for (i in 1:length(glist)) {
      center <- offset - spacing*i
      top <- center + height/2
      bot <- center - height/2
        
      sell <- which(f1$gene==glist[i])
        
      plt <- plt + geom_segment(dat=f1[sell,], 
                                aes_string(x="residual", xend="residual", y=top, yend=bot), 
                                size=0.5, colour="grey30")
      plt <- plt + geom_segment(dat=data.frame(name="mean",residual=f1$meanres[sell[1]]), 
                                aes_string(x="residual", xend="residual", y=top, yend=bot), 
                                size=1, colour="red")
      xlim <- max(xlim, max(abs(f1$residual[sell])))
    }
    plt <- plt + scale_x_continuous(limits=c(-1*xlim,xlim))
    
    if (!is.na(out)) {
        ggsave(paste0(out,name,".eps"),plot = plt,width = 6, height = 2.5+0.5*length(glist), dpi = 300)
    }
    return(plt)
}

crispr_volcano <- function(df,name="",x="z",c1=c(),c2=c(),sig_t=0.05,labelc1=F) {
    #df <- genes
    #id <- "Gene"
    
    options(repr.plot.width=10, repr.plot.height=8)

    df$sig <- df$fdr <= sig_t#0.05
    df$fdr[which(df$fdr == 0)] <- min(df$fdr[which(df$fdr != 0)])

    proviral <- which(df$sig & df[[x]] > 0)
    antiviral <- which(df$sig & df[[x]] < 0)

    df$label <- ""

    if (labelc1) {
      df$label[which(df$Gene %in% c1)] <- df$Gene[which(df$Gene %in% c1)]
    } else {
      df <- df %>% arrange(z)
      df$label[1:10] <- df$Gene[1:10]
      df <- df %>% arrange(desc(z))
      df$label[1:10] <- df$Gene[1:10]
    }    

    df$color <- unlist(lapply(1:nrow(df),function(x){
        if (!df$sig[x]) return("ns")
        if (df$Gene[x] %in% c1) return("sars2 antiviral")
        if (df$Gene[x] %in% c2) return("sars2 proviral")
        return("sig")
    }))
    colorlevels <- c("ns","sig","sars2 antiviral","sars2 proviral")
    df$color <- factor(df$color, levels=colorlevels)#,"sig"))
    
    #red <- c("grey80",'#FB9A99','#E31A1C')
    #orange <- c("grey80",'#FDBF6F','#FF7F00')
    #'#313695' '#436EB0' '#6AA1CB' '#9ACBE1' '#CAE8F2' '#EFE9C3' '#FDCB7D' '#FA9A57' '#EE603D' '#D12B26' '#A50026'
    colors <- list("grey80","black",'#313695','#A50026')
    #colors <- list("grey80","black",'#FF7F00','#E31A1C')
    names(colors) <- colorlevels
    
    plt <- ggplot(df,aes_string(x=x,y="-log10(fdr)",label="label",color="color")) +
        theme_bw() +
        geom_point(size=3)+#data=df[which(df$sig),],size=3) +
        #geom_point(data=df[which(!df$sig),],size=3,color="grey80") +
        scale_y_continuous(expand = expansion(mult = c(0, .15))) +
        ggtitle(paste0(name,"\n",length(proviral)," proviral\n",length(antiviral)," antiviral")) +
        #ggtitle(paste0(f,"\n",length(proviral)," proviral\n",length(antiviral)," antiviral")) +
        xlab("infected vs mock LFC zscore") +
        scale_color_manual(values=colors) +
        theme(panel.border = element_rect(colour = "black", fill=NA, size=2))+
        theme(text = element_text(size=20)) +
        geom_label_repel() +
        theme(aspect.ratio=1)
    
    return(plt)
}

normalize_counts <- function(df) {
    # log2((#of reads for guide/total reads in condition) x 1e6 + 1)
    df[,3:ncol(df)] <- log2( (df[,3:ncol(df)]*1000000 / colSums(df[,3:ncol(df)])) + 1 )
    return(df)
}

create_condition_list <- function(btable,merge_reps=TRUE) {
    
    tr <- list()
    
    # to merge reps
    if (merge_reps) {
        allvirus <- unique(unlist(lapply(colnames(btable)[3:ncol(btable)],
                                         function(x) substr(x,1,nchar(x)-2) )))
    } else {
        allvirus <- colnames(btable)[3:ncol(btable)]
    }
    mockprefix <- paste0(strsplit(allvirus,"_")[[1]][1],"_mock")
    allvirus <- allvirus[which(!startsWith(allvirus,mockprefix))]
                                     
    for (virus in allvirus) {
        
        tr[[virus]] <- data.frame("sgRNA"=c(),"Gene"=c(),"treatgrp"=c(),"ctrlgrp"=c())
        
        if (!merge_reps) {
            if (endsWith(virus,"_1")) suff <- "_1"
            if (endsWith(virus,"_2")) suff <- "_2"
            #tr[[virus]] <- btable[,c("sgRNA","Gene",virus,paste0(mockprefix,suff))]
            
            tr[[virus]] <- data.frame(
                "sgRNA"=btable$sgRNA,
                "Gene"=btable$Gene,
                "treatgrp"=btable[[virus]],
                "ctrlgrp"=rowMeans(btable[,c(paste0(mockprefix,"_1"),paste0(mockprefix,"_2"))]),
                stringsAsFactors=F
            ) 
            
            next
        }
        
        for (suff in c("_1","_2")) {
            t <- paste0(virus,suff)
            if (t %in% colnames(btable)) {
                #print(t)
                tr[[virus]] <- rbind(
                    tr[[virus]], 
                    data.frame(
                        "sgRNA"=paste0(btable$sgRNA,suff),
                        "Gene"=btable$Gene,
                        "treatgrp"=btable[[t]],
                        "ctrlgrp"=rowMeans(btable[,c(paste0(mockprefix,"_1"),paste0(mockprefix,"_2"))]),
                        stringsAsFactors=F
                    )
                )
            }
        }
    }
    
    return(tr)
}

guide_zscores <- function(df_list,treat="treatgrp",ctrl="ctrlgrp") {
    
    for (virus in names(df_list)) {

        df_list[[virus]]$ctrl <- startsWith(df_list[[virus]]$Gene,"CTRL")

        df_list[[virus]]$lfc <- df_list[[virus]][[treat]] - df_list[[virus]][[ctrl]]

        ave_ctrl <- mean(df_list[[virus]]$lfc[which(df_list[[virus]]$ctrl)])
        std_ctrl <- sd(df_list[[virus]]$lfc[which(df_list[[virus]]$ctrl)])

        df_list[[virus]]$z <- (df_list[[virus]]$lfc - ave_ctrl) / std_ctrl

        df_list[[virus]]$pvalue <- 2*pnorm(-abs(df_list[[virus]]$z))
        df_list[[virus]]$fdr <- p.adjust(df_list[[virus]]$pvalue)
    }
    return(df_list)
}

gene_zscores <- function(df_list) {
    genes <- list()
    for (virus in names(df_list)) {
        guides <- df_list[[virus]]
        
        genes[[virus]] <- data.frame(Gene=unique(guides$Gene),stringsAsFactors=F)
        genes[[virus]]$z <- unlist(lapply(1:nrow(genes[[virus]]),function(i){
            mtch <- which(guides$Gene == genes[[virus]]$Gene[i])
            return( mean(guides$z[mtch])*length(mtch) )
        }))
        genes[[virus]]$pvalue <- 2*pnorm(-abs(genes[[virus]]$z))
        genes[[virus]]$fdr <- p.adjust(genes[[virus]]$pvalue)
    }
    return(genes)
}
                                     
create_zscore_table <- function(genes,mb="Gene") {
    z_df_list <- list()
    
    for (virus in names(genes)) {
        z_df_list[[virus]] <- genes[[virus]][,c(mb,"z","fdr")]
        colnames(z_df_list[[virus]]) <- c(mb,paste0(virus,"_z"),paste0(virus,"_fdr"))
    }
    
    z_all <- Reduce(function(x, y) merge(x, y,by=mb, all=TRUE), z_df_list)
    return(z_all)
}


