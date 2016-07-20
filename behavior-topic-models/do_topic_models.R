## This script will use topic models to check how many different behaviors
## are among beetles attracted to lights
library(topicmodels)
library(mcmcse)
library(RColorBrewer)


#first, read data
dados <- read.csv("../all_data.csv")

#now, table counts by species, ignoring date
counts <- xtabs(formula = ~paste(higher_taxon,species,sep="_")+trap, data = dados)


#now, divide data in training and testing datasets
training_rows <- sample(dim(counts)[1], ceiling(dim(counts)[1] * 0.75), replace = F)
counts_training <- counts[training_rows,]
counts_testing <- counts[-training_rows,]

#checking how many mcmc samples needed
K = 3 
LDA_test <- LDA(x = counts_training,k = K ,method = "Gibbs",control = list(verbose = TRUE, alpha=7/K, best=FALSE, burnin=1000, iter=10000, thin = 10))
plot(logLik(LDA_test),type="l")
ess(matrix(logLik(LDA_test), ncol=1))
betas = sapply(LDA_test@fitted,function(x){x@beta})
for (i in 1:12){
  plot(betas[i,],type='l')
}
ess(t(betas))


###estimate number of topics using perplexity
alpha_prior = 15
LDA_results = list()
for (K in 2:10){
  LDA_results[K-1] <- LDA(x = counts_training,k = K ,method = "Gibbs",control = list(verbose = TRUE, alpha=alpha_prior/K, burnin=1000, iter=10000, thin = 100))
}

plot(2:10,unlist(lapply(LDA_results,logLik)))
plot(2:10,unlist(lapply(LDA_results,perplexity,counts_testing)))

#even though perplexity suggests 4 topics, a look at the posteriors reveals overfitting
#three topics
posterior(LDA_results[[2]])$terms
#four topics
posterior(LDA_results[[3]])$terms
#will consider three topics from now on
#running final analysis on all data
K = 3
LDA_final = LDA(x = as.matrix(counts),
                k = K,
                method = "Gibbs",
                control = list(verbose = TRUE, alpha=alpha_prior/K, best=TRUE, burnin=100000, iter=1000000, thin = 10))

#### for each behavor ('topic'), find all species ('documents') with > 0.50 posterior probability
large_pp <- posterior(LDA_final)$topics > 0.5
species_behavior = apply(large_pp,2,function(x){rownames(large_pp)[x]})
behavior_counts = lapply(species_behavior, function(x){table(unlist(strsplit(x,'_[0-9]*')))})

# create color palettes
# beetle familes (ordered by number of species)

familes = sapply(strsplit(names(unlist(behavior_counts)),'[0-9]*\\.'),function(x){x[2]})
fam_counts = unlist(behavior_counts)
names(fam_counts) = sapply(strsplit(names(unlist(behavior_counts)),'[0-9]*\\.'),function(x){x[2]})
fam_counts = sort(tapply(fam_counts,names(fam_counts),sum),decreasing = T)

####The code below may vary depending on the number of species found
####If more than 12 families, another pallete will be needed
behavior_colors = brewer.pal(n = length(fam_counts),name = 'Set3')
names(behavior_colors) = names(fam_counts)

# traps
plot_colors = brewer.pal(4,'Dark2')
names(plot_colors) = c('Hg','Na','Na_F','Control')

# since the order of the results from LDA is arbitrary, the following code will order it for repeatability
hg = order(posterior(LDA_final)$terms[,'Hg'],decreasing = T)[1]
na = order(posterior(LDA_final)$terms[,'Na'],decreasing = T)[1]
order_to_plot = c(hg,na,setdiff(1:3,c(hg,na)))

pdf('behavior_counts.pdf',width = 6,height = 4)
par(mfrow=c(2,3), 
    omi = c(0.3,0.5,0.1,1), 
    mai = c(0,0.2,0,0.25),
    cex.axis = 0.9,
    cex.lab = 0.8)
pie(posterior(LDA_final)$terms[order_to_plot[1],],
    radius = 1,
    labels=NA,
    col = plot_colors[names(posterior(LDA_final)$terms[order_to_plot[1],])],
    edges = 1000)
mtext('Posterior probability\nof inferred behaviors',side = 2,line = 2.5,cex=0.8)

for (i in 2:3){
  pie(posterior(LDA_final)$terms[order_to_plot[i],],
      radius = 1,
      labels=NA,
      col = plot_colors[names(posterior(LDA_final)$terms[order_to_plot[i],])],
      edges=1000)
}


temp = sort(behavior_counts[[order_to_plot[1]]],decreasing = T)
barplot(temp, 
        col = behavior_colors[names(temp)],
        cex.names = 0.5,
        ylim=c(0,max(fam_counts)),
        xaxt='n',
        xlim=c(0,max(sapply(behavior_counts, length))))
mtext('Number of species',side = 2,line = 2.5,cex=0.8)

for (i in 2:3){
  temp = sort(behavior_counts[[order_to_plot[i]]],decreasing = T)
barplot(temp, 
        col = behavior_colors[names(temp)],
        cex.names = 0.5,
        ylim=c(0,max(fam_counts)),
        yaxt='n',
        xaxt='n',
        xlim=c(0,max(sapply(behavior_counts, length))))
}

mtext('Behaviors and beetle families',
      side = 3,
      line = -1,
      outer=T,
      cex=1)

par(fig = c(0, 1, 0, 1), omi = c(0,0,0,0), mar = c(0, 0, 0, 0), new = TRUE)
plot(0, 0, type = "n", bty = "n", xaxt = "n", yaxt = "n", xlim=c(0,1000),ylim=c(0,1000))
legend(x= 830, y = 900 , legend = names(plot_colors),
       xpd = TRUE, inset = c(0, 0), 
       bty = "n", fill = plot_colors, cex = 0.9)
legend(x= 830, y = 500 , legend = names(behavior_colors),
       xpd = TRUE, inset = c(0, 0), 
       bty = "n", fill = behavior_colors, cex = 0.9)
rect(xleft = 35, xright = 320, ybottom = 0, ytop = 920)
rect(xleft = 335, xright = 590, ybottom = 0, ytop = 920)
rect(xleft = 605, xright = 830, ybottom = 0, ytop = 920)
text(labels = 1:3,
     x = c(45,345,615),
     y = 900)

dev.off()
