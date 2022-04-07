###R studio version 1.3.1093
###R version 4.0.3 (2020-10-10) "Bunny-Wunnies Freak Out"

rm(list=ls())

#installing and importing required packages
if   (!require("pacman")) install.packages("pacman")
pacman::p_load("RColorBrewer", 
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
library(RColorBrewer)
library(data.table)
library(purrr)

#setting working directory --> where the input csv file is
#setwd("~/Desktop/Projects/Other projects/Glycoproject/") #Provide the path to the working directory 
#setwd(choose.dir()) #Interactively navigate to the working directory 

#Importing csv file and modifying column names; Adjust the csv file name
data              <- read.csv(file= "Control_MRMNormalizeR.csv", 
                              header=TRUE, 
                              sep=",", 
                              stringsAsFactors = F, 
                              as.is = T, 
                              na.strings = 'NA')

data$Condition         <- gsub(pattern = "_", replacement = ".",data$Condition)
data$Protein.Name      <- gsub(pattern = "^.*\\|", replacement = "", data$Protein.Name)###should be  modified depending on how protein name look like
class.ann.list         <- read.delim("Annotation.txt", stringsAsFactors = F)
#dim(data)

#Filtering data and sum up all areas available for each protein in each measurement (Total Area)
filtered.data <-   data %>%  
              filter(Truncated != "" & Quantitative == "True" & Standard.Type != "iRT") %>% #use if you want to keep quantitative but truncated areas
              #filter(Truncated != "" & Quantitative == "True" & Standard.Type != "iRT") %>%   #use if you want to get rid of quantitative but truncated areas
              select(-c(Truncated, Quantitative, Standard.Type))  %>%
              group_by(Protein.Name, Replicate.Name, Condition, TechReplicate, Isotope.Label.Type) %>%
              summarise(Total.Area = sum(as.numeric(Area, na.rm = TRUE)))



#dim(filtered.data)
#log transformation    
filtered.data$Total.Area <- log2(filtered.data$Total.Area)


#calculating light to heavy ratios
ratios.data <-   filtered.data %>% 
              reshape2::dcast(Protein.Name + Replicate.Name + Condition + TechReplicate ~ Isotope.Label.Type, 
              value.var = "Total.Area", 
              fill = 0) %>%
              mutate(Ratio = light - heavy) %>%
              select(Protein.Name, Replicate.Name, Condition, TechReplicate, Ratio)


#dim(ratios.data)

#normalizing against housekeeping / reference protein
wide.data <-  ratios.data %>% 
              reshape2::dcast(Replicate.Name + Condition + TechReplicate ~ Protein.Name, 
              value.var = "Ratio", 
              fill = 0) 

select          <- colnames(wide.data[,-1:-3]) #To artificially get rid of first three columns 

# write.table(wide.data, file = 'Randomname.txt',
#              quote = F, sep = '\t', dec = ',', row.names = F )

wide.data[, select]  <- wide.data[,select] - wide.data[,"SEC63_HUMAN"] #second normalization step with reference protein Sec63 in that example
data.normalized <- data.frame(Replicate.Name = wide.data$Replicate.Name, 
                              Condition = wide.data$Condition, 
                              Technical.Rep = wide.data$TechReplicate,
                              wide.data[,select])

# write.table(data.normalized, file = 'Randomname.txt',
#             quote = F, sep = '\t', dec = ',', row.names = F )

#Extract replicate number from file name and make a new column Replicate, necessary to work for possible technical replicates annotation, otherwise just run regardless of what gets extracted 
data.normalized$Replicate <- str_extract(data.normalized$Replicate.Name, pattern = "[0-9]")
#data.normalized$Condition <- gsub(pattern = ".*\\.", replacement = "", data.normalized$Condition) #optional

#ordering conditions in a specific order   --> needed for group comparison, control condition must be put in first!!
order <- c("Control", "HeLa","Fibroblasts")

data.normalized <- data.normalized %>% 
                   gather(key = Protein.Name, 
                   value = Log2Ratio, 
                   -c(Replicate.Name, Condition, Replicate, Technical.Rep)) %>%
  #group_by(Protein.Name, Replicate, Condition, Technical.Rep) %>%  ###Tech rep. part, needs testing, use at own risk for now J
  #summarise(Replicate.Name = first(Replicate.Name), Log2Ratio = mean(Log2Ratio)) %>%
  #ungroup() %>%
                   mutate_at(vars(Condition), 
                   list(factor)) %>% 
                   mutate(Condition = factor(Condition, levels = order))  %>% #factor makes it to preserve the order, not the string alphabetical one 
                   group_split(Condition) #previous data frame is now a tibble, one data frame per condition


##linear regression 
control <- data.normalized[[1]]   ### according to order which is defined above, only control condition data inside

coef.list <- list()
dof.list  <- list()
##starts with 2 to eliminate control to control comparison, 
##therefore it is important to have the first  condition as control (see above)
for (i in 2:length(order)){ 
        
        control.condition <- rbind(control, data.normalized[[i]]) #combines control data with data from other Conditions, pairwise
        fitted.models     <- control.condition %>% 
                             group_by(Protein.Name) %>% 
                             do(model = lm(Log2Ratio ~ Condition, data = .)) #linear regression

#Now, results from linear regression for pairwise comparisons must be stored in lists, as they would otherwise be overwritten by the next comparison        
        
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

##preparation for storing results of linear regression
fc.results       <- data.frame(Protein.Name = unique(coef.list[[ i ]]$Protein.Name), stringsAsFactors=FALSE)
adj.pval.results <- data.frame(Protein.Name = unique(coef.list[[ i ]]$Protein.Name), stringsAsFactors=FALSE)
pval.results     <- data.frame(Protein.Name = unique(coef.list[[ i ]]$Protein.Name), stringsAsFactors=FALSE)
conf.int.results <- data.frame(Protein.Name = unique(coef.list[[ i ]]$Protein.Name), stringsAsFactors=FALSE)
sig.results      <- data.frame(Protein.Name = unique(coef.list[[ i ]]$Protein.Name), stringsAsFactors=FALSE)
#barplots         <- list()
# combining protein values from different comparisons in one dataframe
for (i in 2:length(order)){
        
        if( is.null( coef.list[[i]] ) ){next()}
#writing results from the linear regression back from lists to data frames       
        coef.table             <-  coef.list[[i]]
        dof.table              <-  dof.list[[i]]
        
        coef.all               <-  data.frame(cbind(coef.table[paste0("Condition", order[i])],
                                                 coef.table["std.error"] ,
                                                 coef.table["p.value"]))
        
        coef.all            <-  subset(coef.all, is.na(coef.all[paste0("Condition", order[i])]) == FALSE) 
        row.names(coef.all) <-  unique(coef.table$Protein.Name)
        
        #Performing the BH correction
        coef.all$adj.p.value <-  p.adjust(coef.all$p.value, method = "BH", n = length(coef.all$p.value))
        coef.all$sig         <-  ifelse(coef.all$adj.p.value <= 0.05, "+", "-")
        
        #populating separate data frames
        fc.results[paste("log2FC", order[i], sep = ".")]          <-  coef.all[,paste0("Condition", order[i])]
        adj.pval.results[paste("adj.p.val", order[i], sep = ".")] <-  coef.all[,"adj.p.value"]
        pval.results[paste("p.val", order[i], sep = ".")]         <-  coef.all[,"p.value"]
        conf.int.results[paste("conf.int", order[i], sep = ".")]  <-  coef.all[,"std.error"]*qt(.975, as.numeric(unique(dof.table["df.residual"]))) #calculating 95 % confidence intervals
        sig.results[paste("is.sign", order[i], sep = ".")]        <-  coef.all[,"sig"]
}


results     <- Reduce(merge, list(fc.results, adj.pval.results, pval.results, sig.results, conf.int.results ))
results     <- results[!(results$Protein.Name=="SEC63_HUMAN"),] ##removing SEC63 from the results

#exporting Log2FC, Adj-p value and CI results in txt file
write.table( results, file = '22.02.07_Result.txt',
             quote = F, sep = '\t', dec = ',', row.names = F )


#############################################################################################
##barplots with log2FC and 95% CI as error bars for individual proteins across all conditions
##in single pdfs which are saved in working directory

plot <- list()
results$Protein.Name <- gsub(pattern = "\\_.*", replacement = "", results$Protein.Name)

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
        dm                   <- reshape2::melt( protein.results, id.vars = "Protein.Name")
        dm$condition         <- gsub(pattern = ".*\\.", replacement = "", dm$variable)
        dm$variable          <- gsub(pattern = "\\..*", replacement = "", dm$variable)
        dm                   <- spread(dm, key = "variable", value = "value")
        dm$log2FC            <- as.numeric(as.character(dm$log2FC))
        dm$conf              <- as.numeric(as.character(dm$conf))
        dm$condition         <- factor(dm$condition,
                                       levels = c("Control", "HeLa", "Fibroblasts")) #order does not matter here
        
  
        plot[[p]]  <- protein.bar.plot(data     = dm, 
                                       condition= dm$condition, 
                                       value    = dm$log2FC,
                                       class    = dm$is,
                                       errbar   = dm$conf, 
                                       pn       = p)


        pdf(sprintf("%s.pdf", p),
             width = 6, height = 4, onefile = T)
             plot(plot[[p]])
        dev.off()
        
}

##################################################################################################
##barplots with log2FC and 95% CI as error bars with all proteins per condition as separate single pdfs

p = list()
results$Protein.Name <- gsub(pattern = "\\_.*", replacement = "", results$Protein.Name) #getting rid of "_HUMAN"
results              <- results[ order(match(results$Protein.Name, class.ann.list$Protein.Name)), ] #custom protein order 
results$Class        <- class.ann.list[ match(results$Protein.Name, class.ann.list$Protein.Name), "Class" ] 
my_colors            <- c('#b2abd2', '#e66101', '#fdb863', '#5e3c99' )
plot_data_column = function (data, variable, value, sign, class, errbar){
        ggplot(data = data, 
               aes(x = variable, 
                   y = value,
                   fill = factor(class))) +
                geom_bar(stat  = "identity", 
                         color = "black",
                         width = .8, 
                         position = position_identity()) +
                scale_fill_manual(values  = my_colors) + 
                scale_x_discrete(limits = c(levels(variable)[1:table(class)[1]], "ABC",
                                            levels(variable)[(1+table(class)[1]):(table(class)[1]+table(class)[2])], "DEF",
                                            levels(variable)[(1+table(class)[1]+table(class)[2]):(table(class)[1]+table(class)[2]+table(class)[3])], "KLM",
                                            levels(variable)[(1+table(class)[1]+table(class)[2]+table(class)[3]):length(variable)]),
                     labels = c("ABC" = "",
                                "DEF" = "", 
                                "KLM" = ""))+
                scale_y_continuous(breaks = seq(-ceiling(max(abs(value))) - 2, 
                                                 ceiling(max(abs(value))) + 2, 
                                                by = 1),
                                   limits = c(-ceiling(max(abs(value))) - 2, 
                                               ceiling(max(abs(value))) + 2)) +
                geom_errorbar(aes(ymin = value - errbar, 
                                  ymax = value + errbar), 
                              width = .2,
                              position = position_dodge(.9)) +
                geom_hline(yintercept = 0) +
                geom_hline(yintercept =  1, linetype = 2) +
                geom_hline(yintercept = -1, linetype = 2) +
                geom_text(aes(label = ifelse(sign == "+", "*", "")), 
                          nudge_y = ifelse(value <0, -(errbar + 0.5), errbar + 0.5),
                          vjust = 1) +
                theme_bw() + 
                theme(plot.title = element_text(face = "bold",
                                                size = rel(1), hjust = 0),
                      text = element_text(),
                      panel.background = element_rect(colour = NA),
                      plot.background  = element_rect(colour = NA),
                      panel.border     = element_rect(colour = NA),
                      axis.title       = element_text(face = "bold", size = rel(1)),
                      axis.title.y     = element_text(angle = 90, vjust = 2),
                      axis.title.x     = element_text(vjust = -.2),
                      axis.text        = element_text(), 
                      axis.text.x      = element_text(angle = 60, face = "bold", vjust = 1, hjust = 1),
                      axis.line        = element_line(colour = "black"),
                      axis.ticks       = element_line(),
                      panel.grid.major = element_blank(),
                      panel.grid.minor = element_blank(),
                      legend.key       = element_rect(colour = NA),
                      legend.position  = c(.8,.9),
                      legend.direction = "vertical",
                      legend.key.size  = unit(.4, "cm"),
                      legend.spacing   = unit(0, "cm"),
                      legend.title     = element_blank(),
                      plot.margin      = unit(c(10,5,5,5),"mm"),
                      strip.background = element_rect(colour = "#f0f0f0",fill = "#f0f0f0"),
                      strip.text       = element_text(face="bold")) +
                labs(title = paste0("Log2FC \n",order[i], " / HEK 293T"),
                     fill = "Class",
                     x = "", 
                     y = "")
}




for (i in 2:length(order)){
        
        proteinNames  <- factor(results$Protein.Name, levels=unique(results$Protein.Name))
        Log2FCs       <- results[,paste0("log2FC", sep = ".", order[i])]
        significance  <- results[,paste0("is.sign", sep = ".", order[i])]
        confIntervals <- results[,paste0("conf.int", sep = ".", order[i])]
        class         <- factor(results$Class, levels=unique(results$Class))
        
        p[[i]]  <- plot_data_column(data     = results, 
                                    variable = proteinNames, 
                                    value    = Log2FCs,
                                    sign     = significance ,
                                    class    = class,
                                    errbar   = confIntervals)
        
        pdf(sprintf("plot_%s.pdf", order[i]),
            width = 4, height = 5, onefile = T)
        plot(p[[i]])
        dev.off()
}

