## ---------------------------
##
## Script name: MRMStaR/MRMnormalizeR
##
## Purpose of script: MRM data second level normalization and group comparison
##
## Author: Merve Alp
##
## Date Created: 2020-09-13
##
## Email: kezibanmerve.alp@mdc-berlin.de
##
## ---------------------------
##
## Notes: 
##   
##
## ---------------------------

## set working directory 

setwd(choose.dir()) 

## install and import required packages

pacman::p_load("pheatmap", 
               "RColorBrewer", 
               "ggplot2",
               "grid", 
               "gridExtra",
               "dplyr",
               "tidyr",
               "reshape2",
               "stringr",
               "broom",
               "miscset",
               "data.table",
               "purrr")



library(ggplot2)
library(dplyr)
library(tidyr)
library(reshape2)
library(stringr)
library(broom)
library(grid)
library(gridExtra)
library(miscset)
library(tibble)
library(gplots)
library(corrplot)
library(RColorBrewer)
library(ggsci)
library(ggpubr)
library(data.table)
library(purrr)

## Import input files

data              <- read.csv(file= "Input/Input_Group_Comparison_Normalization.csv", 
                              header = TRUE, 
                              sep = ",", 
                              stringsAsFactors = F, 
                              as.is = T,
                              na.strings = 'NA')

annotation.list   <- read.delim("Input/Annotation_Paper_woClass.txt", stringsAsFactors = F)
class.ann.list    <- read.delim("Input/Annotation_Paper.txt", stringsAsFactors = F)

## Modify column names of data for convenience

data$Replicate.Name    <- gsub(pattern = "^[^_]*_", replacement = "",data$Replicate.Name)###removes the date from the beginning of the names, optional
data$Replicate.Name    <- gsub(pattern = "_", replacement = ".",data$Replicate.Name)###optional
data$Condition         <- gsub(pattern = "_", replacement = ".",data$Condition)
data$Protein.Name      <- gsub(pattern = "^.*\\|", replacement = "", data$Protein.Name)###should be  modified depending on how protein name look like


## Filter truncated, non-quantitave and iRT transitions

data <- data %>%  
        filter(Truncated == "False" & Quantitative == "True" &  Standard.Type != "iRT") %>%
        select(-c(Truncated, Quantitative, Standard.Type))

## Sum up all areas available for each protein in each wiff file (Total Area)
        
data <- data %>%        
        group_by(Protein.Name, Replicate.Name, Condition, BioReplicate, Isotope.Label.Type) %>%
        summarise(Total.Area = sum(as.numeric(Area, na.rm = TRUE)))

## Log transformation    

data$Total.Area <- log2(data$Total.Area)

## Calculate light to heavy ratios

data <- data %>% 
        dcast(Protein.Name + Replicate.Name + Condition + BioReplicate ~ Isotope.Label.Type, 
              value.var = "Total.Area", 
              fill = 0) %>%
        mutate(Ratio = light - heavy) %>%
        select(Protein.Name, Replicate.Name, Condition, BioReplicate, Ratio)

## Normalize against houskeeping protein

data <- data %>% 
        dcast(Replicate.Name + Condition + BioReplicate ~ Protein.Name, 
              value.var = "Ratio", 
              fill = 0) 

select          <- grep("_HUMAN", names(data)) ## selecting only columns with protein values
#select         <- data[,-1:-2]  ## if your protein names do not have _HUMAN, continue with this select command

ctrl            <- "SEC63_HUMAN"
data[, select]  <- data[,select] - data[, ctrl]


## Create a new data frame with normalized values

data.normalized <- data.frame(Replicate.Name = data$Replicate.Name, 
                              Condition = data$Condition, 
                              Technical.Rep = data$BioReplicate, 
                              data[,select])

## Extract replicate number from file name (Replicate.Name) and condition name from Condition columns respectively

data.normalized$Replicate <- str_extract(data.normalized$Replicate.Name, pattern = "R[0-9]")
data.normalized$Condition <- gsub(pattern = ".*\\.", replacement = "", data.normalized$Condition)

## Order conditions in a specific order in a way that control condition is the first one   --> needed for group comparison

order <- c("Hek293T", "HeLa", "Fibroblasts")

## Wrangle data in order to combine technical replicates by averaging if there is any 
## Split each condition 

data.normalized <- data.normalized %>% 
                   gather(key = Protein.Name, 
                          value = Log2Ratio, 
                          -c(Replicate.Name, Condition, Replicate, Technical.Rep)) %>%
                   group_by(Protein.Name, Replicate, Condition, Technical.Rep) %>%  ###Tech rep. part
                   summarise(Replicate.Name = first(Replicate.Name), Log2Ratio = mean(Log2Ratio)) %>%
                   ungroup() %>%
                   mutate_at(vars(Condition), 
                             funs(factor)) %>% 
                   mutate(Condition = factor(Condition, levels = order))  %>% 
                   group_split(Condition)


## Linear regression 

control <- data.normalized[[1]]   

coef.list <- list()
dof.list  <- list()

## Start with 2 to eliminate control to control comparison

for (i in 2:length(order)){ 
        
        control.condition <- rbind(control, data.normalized[[i]])
        fitted.models     <- control.condition %>% 
                             group_by(Protein.Name) %>% 
                             do(model = lm(Log2Ratio ~ Condition, data = .)) 
        
        coef.list[[i]]    <- fitted.models %>% 
                             ungroup %>% 
                             transmute(Protein.Name, Coef = map(model, tidy)) %>% 
                             unnest(Coef) %>% 
                             spread(key = "term", value = "estimate")
        
        
        dof.list[[i]]     <- fitted.models %>% 
                             ungroup %>% 
                             transmute(Protein.Name, Dof = map(model, glance)) %>% 
                             unnest(Dof)              
}

## Store results of linear regression out of a tibble in df

fc.results       <- data.frame(Protein.Name = unique(coef.list[[ i ]]$Protein.Name), stringsAsFactors=FALSE)
adj.pval.results <- data.frame(Protein.Name = unique(coef.list[[ i ]]$Protein.Name), stringsAsFactors=FALSE)
pval.results     <- data.frame(Protein.Name = unique(coef.list[[ i ]]$Protein.Name), stringsAsFactors=FALSE)
conf.int.results <- data.frame(Protein.Name = unique(coef.list[[ i ]]$Protein.Name), stringsAsFactors=FALSE)
sig.results      <- data.frame(Protein.Name = unique(coef.list[[ i ]]$Protein.Name), stringsAsFactors=FALSE)

## Combine protein values from different comparisons in one dataframe

for (i in 2:length(order)){
        
        if( is.null( coef.list[[i]] ) ){next()}
        
        coef.table             <-  coef.list[[i]]
        dof.table              <-  dof.list[[i]]
        coef.table$adj.p.value <-  p.adjust(coef.table$p.value, method = "BH", n = length(coef.table$p.value))
        coef.table$sig         <-  ifelse(coef.table$adj.p.value <= 0.05, "*", "")
        
        coef.all            <-  data.frame(cbind(coef.table[paste0("Condition", order[i])],
                                                 coef.table["std.error"] ,
                                                 coef.table["adj.p.value"],
                                                 coef.table["p.value"],
                                                 coef.table["sig"]))
        
        coef.all            <-  subset(coef.all, is.na(coef.all[paste0("Condition", order[i])]) == FALSE) 
        row.names(coef.all) <-  unique(coef.table$Protein.Name)
        
        fc.results[paste("log2FC", order[i], sep = ".")]          <-  coef.all[,paste0("Condition", order[i])]
        adj.pval.results[paste("adj.p.val", order[i], sep = ".")] <-  coef.all[,"adj.p.value"]
        pval.results[paste("p.val", order[i], sep = ".")]         <-  coef.all[,"p.value"]
        conf.int.results[paste("conf.int", order[i], sep = ".")]  <-  coef.all[,"std.error"]*qt(.975, as.numeric(unique(dof.table["df.residual"])))
        sig.results[paste("is.sign", order[i], sep = ".")]        <-  coef.all[,"sig"]
}


results              <- Reduce(merge, list(fc.results, adj.pval.results, pval.results, sig.results, conf.int.results ))
#results             <- results[!(results$Protein.Name=="SEC63_HUMAN"),] ##removing SEC63 from the results
results$Protein.Name <- gsub(pattern = "\\_.*", replacement = "", results$Protein.Name)
results              <- results[ order(match(results$Protein.Name, annotation.list$Protein.Name)), ]

## Export Log2FC, Adj-p value and CI results in txt file

write.table( results, file = 'Output/Stat_results.txt',
             quote = F, sep = '\t', dec = '.', row.names = F ) 


## ---------------------------
##
## Visualization with plots
##   
##
## ---------------------------

## Convert necessary columns to matrix for heatmap

data.matrix              <-  as.matrix(results[,grep("log2FC.", names(results))])
row.names(data.matrix)   <-  gsub(pattern = "\\_.*", replacement = "", results$Protein.Name) 
colnames(data.matrix)    <-  gsub(pattern = ".*\\.", replacement = "", colnames(data.matrix))

## Row annotation

group.ann            <- as.matrix(row.names(data.matrix))
colnames(group.ann)  <- "id" 
group.ann            <- merge(group.ann, class.ann.list, by.x = "id", by.y = "Protein.name", all.x = F) 
group.ann            <- group.ann %>% 
                        column_to_rownames(var="id")


## Heat map

paletteLength <- 100
breaks        <- c(seq(-5, 0, length.out=ceiling(paletteLength/2)-1 ), 
                   seq(2/paletteLength, 5, length.out=floor(paletteLength/2)))


plot.hm <- pheatmap(
        data.matrix,
        main              = "Comparison_to_control",
        color             = colorRampPalette(c("blue","white","red"))(paletteLength),
        breaks            = breaks,
        annotation_row    = group.ann,
        cluster_rows      = T,
        cluster_cols      = F,
        show_colnames     = T,
        show_rownames     = T,
        fontsize_row      = 11, 
        cellwidth         = 15,
        cellheight        = 15,
        treeheight_row    = 20,
        treeheight_col    = 20
)


save_pheatmap_pdf <- function(x, filename, width=7, height=7) {
        stopifnot(!missing(x))
        stopifnot(!missing(filename))
        pdf(filename, width=width, height=height)
        grid::grid.newpage()
        grid::grid.draw(x$gtable)
        dev.off()
}
save_pheatmap_pdf(plot.hm, "Output/Heatmap.pdf")

## Barplots with log2FC and 95% CI as error bars for individual proteins across all conditions
## in single pdfs which are saved in working directory

plot <- list()
protein.bar.plot <- function(data, condition, value, class, errbar, pn){
        ggplot(data, aes(x = condition, y= value, fill = class)) + 
                geom_bar(stat = "identity", 
                         color = "black",
                         position = position_dodge())  +
                geom_errorbar(aes(ymin = value - errbar, 
                                  ymax = value + errbar), 
                              width = .2,
                              position = position_dodge(.9))+
                geom_hline(yintercept = 0) +
                theme_bw() +
                theme(axis.line.x = element_blank()) +
                scale_fill_jco()+
                #ylim(-2, 2)+
                labs(  fill = "adj p-value <= 0.05",
                       x = "Condition", 
                       y = "Log2FC") +
                ggtitle(pn)
}

protein.list <- unique(results$Protein.Name)
for (p in protein.list){
        
        
        protein.results      <- results[results$Protein.Name == p,]
        dm                   <- melt( protein.results, id.vars = "Protein.Name")
        dm$condition         <- gsub(pattern = ".*\\.", replacement = "", dm$variable)
        dm$variable          <- gsub(pattern = "\\..*", replacement = "", dm$variable)
        dm                   <- spread(dm, key = "variable", value = "value")
        dm$log2FC            <- as.numeric(as.character(dm$log2FC))
        dm$conf              <- as.numeric(as.character(dm$conf))
        dm$condition         <- factor(dm$condition,
                                       levels = order) 
        
        
        plot[[p]]  <- protein.bar.plot(data     = dm, 
                                       condition= dm$condition, 
                                       value    = dm$log2FC,
                                       class    = dm$is,
                                       errbar   = dm$conf, 
                                       pn       = p)
        
        
        pdf(sprintf("Output/%s.pdf", p),
            width = 6, height = 4, onefile = T)
        plot(plot[[p]])
        dev.off()
        
}

## Barplots with log2FC and 95% CI as error bars with all proteins per condition as separate single pdfs

p = list()
#results <- results[ order(match(results$Protein.Name, annotation.list$Protein.Name)), ]
plot_data_column = function (data, variable, value, class, errbar){
        ggplot(data = data, 
               aes(x = variable, 
                   y = value, 
                   fill = class)) +
                geom_bar(stat     = "identity", 
                         color    = "black",
                         position = position_dodge()) +
                scale_fill_manual(values = c('#999999','#E69F00')) + 
                geom_errorbar(aes(ymin = value - errbar, 
                                  ymax = value + errbar), 
                              width    = .2,
                              position = position_dodge(.9)) +
                geom_hline(yintercept  = 0) +
                theme_bw() + 
                theme(axis.line.x = element_blank(), 
                      axis.text.x = element_text(angle = 60, hjust = 1)) +
                #ylim(-10, 10) +
                labs(title = order[i],
                     fill  = "adj p-value <= 0.05",
                     x     = "Protein name", 
                     y     = "Log2FC") 
}



for (i in 2:length(order)){
        
        proteinNames  <- factor(results$Protein.Name, levels=unique(results$Protein.Name))
        Log2FCs       <- results[,paste0("log2FC", sep = ".", order[i])]
        significance  <- results[,paste0("is.sign", sep = ".", order[i])]
        confIntervals <- results[,paste0("conf.int", sep = ".", order[i])]
        
        p[[i]]  <- plot_data_column(data     = results, 
                                    variable = proteinNames, 
                                    value    = Log2FCs,
                                    class    = significance,
                                    errbar   = confIntervals)
        
        pdf(sprintf("Output/plot_%s.pdf", order[i]),
            width = 6, height = 4, onefile = T)
        plot(p[[i]])
        dev.off()
}



## Barplots with log2FC and 95% CI as error bars with all proteins per condition in one page pdf

dm              <- melt(results, id.vars = "Protein.Name")
#dm             <- dm[order(match(dm$Protein.Name, annotation.list$Protein.Name)),]
dm$condition    <- gsub(pattern = ".*\\.", replacement = "", dm$variable)
dm$variable     <- gsub(pattern = "\\..*", replacement = "", dm$variable)
dm              <- spread(dm, key = "variable", value = "value")
dm              <- dm[order(match(dm$Protein.Name, annotation.list$Protein.Name)),]
dm$log2FC       <- as.numeric(as.character(dm$log2FC))
dm$conf         <- as.numeric(as.character(dm$conf))
dm$condition    <- factor(dm$condition,
                          levels = order)
dm$Protein.Name <- factor(dm$Protein.Name, 
                          levels =  unique(dm$Protein.Name))

conditon_plot <- ggplot(dm, aes(x = Protein.Name, y = log2FC, fill = is)) + 
                        geom_bar(stat     = "identity", 
                                 color    = "black",
                                 position = position_dodge())  +
                        geom_errorbar(aes(ymin = log2FC - conf, 
                                 ymax = log2FC + conf), 
                                 width    = .4,
                                 position = position_dodge(.9))+
                        geom_hline(yintercept  = 0) +
                        scale_fill_jco()+
                        theme_bw() +
                        theme(axis.text.x     = element_text(angle = 90, hjust = 0.5, size = 6),
                              axis.line.x     = element_blank(),
                              legend.title    = element_text(size = 9),
                              legend.key.size = unit(0.2, "cm"), 
                              plot.title      = element_text(hjust = 0.5, size = 12 )) +
                        #ylim(-4, 6.5) + 
                        labs(  title = "Cell_Line_Comparisons",
                               fill  = "adj p-value \n<= 0.05",
                               x     = "Proteins", 
                               y     = "Log2FC")+
                        facet_wrap(~ condition, ncol=2) 

pdf("Output/Proteins in all conditions.pdf",
    width = 6, height = 4, onefile = T)
plot(conditon_plot)
dev.off()

## Barplots with log2FC and 95% CI as error bars with all proteins grouped 


offset_asterisk <- .2

theme_ms <- function(base_size=12, base_family="Helvetica") {
                (theme_bw(base_size   = base_size, 
                          base_family = base_family)+
                 theme(text              = element_text(color = "black"),
                       axis.title        = element_text(face  = "bold", 
                                                        size  = rel(1.3)),
                       axis.text         = element_text(size  = rel(1), 
                                                        color = "black"),
                       legend.title      = element_text(face  = "bold"),
                       legend.text       = element_text(face  = "bold"),
                       legend.background = element_rect(fill  = "transparent"),
                       legend.key.size   = unit(0.8, 'lines'),
                       panel.border      = element_rect(color = "black",
                                                        size  = 1),
                       panel.grid        = element_blank()
                 ))
}

plot <- ggplot(dm, aes(x = Protein.Name, y = log2FC, fill = condition)) + 
        geom_bar(stat     = "identity", 
                 color    = "black",
                 width    = .7,
                 position = position_dodge(.7) )  +
        geom_errorbar(aes(ymin = log2FC - conf, 
                          ymax = log2FC + conf), 
                      width    = .3,
                      position = position_dodge(.7))+
        geom_hline(yintercept  = 0) +
        geom_text( aes(y = ifelse(log2FC > 0, (log2FC + conf) + offset_asterisk, (log2FC - conf ) - offset_asterisk), label = is),
                   colour ="black",
                   vjust = 0,
                   position = position_dodge(width= .7),
                   #hjust = .75,
                   size = 5)+
        scale_fill_jco()+
        theme_ms() +
        theme(axis.text.x     = element_text(angle = 0, hjust = .5, size = 8),
              axis.line.x     = element_blank(),
              legend.title    = element_text(size = 9),
              legend.key.size = unit(0.2, "cm")) +
        # ylim(-3, 5) +
        labs(  fill  = "Cell lines",
               x     = "Proteins", 
               y     = "Log2FC")


pdf("Output/Group comparison in one plot with significance asterisk.pdf",
    width = 15, height = 10, onefile = T)
plot(plot)
dev.off()




