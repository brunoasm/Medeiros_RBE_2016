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
alpha_prior = 7
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

#### for each behavor ('topic'), find all species ('documents') with > 0.50 posterior probability
large_pp <- posterior(LDA_results[[2]])$topics > 0.5
species_behavior = apply(large_pp,2,function(x){rownames(large_pp)[x]})
behavior_counts = lapply(species_behavior, function(x){table(unlist(strsplit(x,'_[0-9]*')))})

# create color palettes
# beetle familes (ordered by number of species)

familes = sapply(strsplit(names(unlist(behavior_counts)),'[0-9]*\\.'),function(x){x[2]})
fam_counts = unlist(behavior_counts)
names(fam_counts) = sapply(strsplit(names(unlist(behavior_counts)),'[0-9]*\\.'),function(x){x[2]})

behavior_colors = rainbow(length(fam_counts))
#c(brewer.pal(n = 12,name = 'Set3'),gray(.20))
names(behavior_colors) = names(sort(tapply(fam_counts,names(fam_counts),sum),decreasing = T))

# traps
plot_colors = brewer.pal(4,'Dark2')
names(plot_colors) = c('Hg','Na','Na_F','Control')

pdf('behavior_counts.pdf',width = 6,height = 4)
par(mfrow=c(2,3), 
    omi = c(0.5,0.5,0.2,1), 
    mai = c(0,0.2,0.2,0.25),
    cex.axis = 0.9,
    cex.lab = 0.8)
pie(posterior(LDA_results[[2]])$terms[1,],
    labels=NA,
    col = plot_colors[names(posterior(LDA_results[[2]])$terms[1,])],
    edges = 1000)
mtext('Posterior probability',side = 2,line = 2.5)

for (i in 2:3){
  pie(posterior(LDA_results[[2]])$terms[i,],
      labels=NA,
      col = plot_colors[names(posterior(LDA_results[[2]])$terms[i,])],
      edges=1000)
}


temp = sort(behavior_counts[[1]],decreasing = T)
barplot(temp, 
        col = behavior_colors[names(temp)],
        cex.names = 0.5,
        ylim=c(0,8),
        xaxt='n')
mtext('Species associated',side = 2,line = 2.5)

for (i in 2:3){
  temp = sort(behavior_counts[[i]],decreasing = T)
barplot(temp, 
        col = behavior_colors[names(temp)],
        cex.names = 0.5,
        ylim=c(0,8),
        yaxt='n',
        xaxt='n')
}

title(main='Behaviors and beetle families',outer=T)

par(fig = c(0, 1, 0, 1), omi = c(0,0,0,0), mar = c(0, 0, 0, 0), new = TRUE)
plot(0, 0, type = "n", bty = "n", xaxt = "n", yaxt = "n", xlim=c(0,1000),ylim=c(0,1000))
legend(x= 830, y = 800 , legend = names(plot_colors),
       xpd = TRUE, inset = c(0, 0), 
       bty = "n", fill = plot_colors, cex = 0.9)
legend(x= 830, y = 500 , legend = names(behavior_colors),
       xpd = TRUE, inset = c(0, 0), 
       bty = "n", fill = behavior_colors, cex = 0.9)

dev.off()

