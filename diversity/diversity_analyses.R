#### This script will do two analyses:
## 1st - run a linear model to test for differences in diversity attracted to lamps
## 2nd - test for differences in phylogenetic clustering in lamps

library(ape) #needed for trees
library(picante) #needed for phylognetic diversity and clustering
library(RColorBrewer) #needed for plot colors
library(lme4) #needed for mixed models
library(car) #needed for ANOVA for mixed models (Wald chi-sq)
library(TeachingDemos) #needed for inset plot


#establish plot colors
plot_colors = brewer.pal(4,'Dark2')
names(plot_colors) = c('Hg','Na','Na_F','Control')

#1 - read data file and trees
dados <- read.csv('../all_data.csv', row.names = 1)
dados$trap <- factor(dados$trap, levels = c('Hg','Na','Na_F','Control'),ordered = T)
trees <- read.tree('../beetle_tree/final_1000_trees.tre')

## model abundance attracted to lights

abundance = as.data.frame.table(xtabs(formula = ~trap+date, data = dados))
abundance_model = glmer(formula = Freq~trap + (1|date), data=abundance, family = 'poisson')

summary(abundance_model)
car::Anova(abundance_model)
plot(fitted(object=abundance_model),residuals(abundance_model))
qqplot(fitted(abundance_model),abundance$Freq)
hist(residuals(abundance_model))



## model diversity attracted to lights
#A - species richness
names(dados)
species_richness <- as.data.frame.table(tapply(X = dados$species,
                                               INDEX = dados[c('trap','date')], 
                                               FUN = function (x) {length(unique(x))}))

species_richness$Freq[is.na(species_richness$Freq)] <- 0
names(species_richness) <- c('trap','date','richness')



rich_model <- glmer(richness~trap+(1|date),family=poisson,data=species_richness)
summary(rich_model)
car::Anova(rich_model)
plot(fitted(object=rich_model),residuals(rich_model))
qqplot(fitted(rich_model),abundance$Freq)
hist(residuals(rich_model))


#B - phylogenetic diversity
#function to get phylogenetic diversity
get_phylo_diver <- function (species,tree){
  if (length(species) == 0) { return(0) } #if no species, diversity is 0
  #standardize tree to height = 1
  tree$edge.length = tree$edge.length/max(branching.times(tree))
  #get phylogenetic diversity
  div = pd(samp = matrix(data=1,ncol=length(unique(species)),dimnames = list('com',unique(species))),
           tree = tree,
           include.root = T)
  return(div$PD)
}

#initiate a dataframe to record all diversity data
diversity_data <- species_richness

#calculate phylogenetic diversity for a hundred trees
for (i in 1:100){
  tree = trees[[i]]
  phylo_diver <- as.data.frame.table(tapply(X = paste(dados$higher_taxon,dados$species,sep="_"),
                                            INDEX = dados[c('trap','date')], 
                                            FUN = function (x) {get_phylo_diver(x,tree) }))
  phylo_diver$Freq[is.na(phylo_diver$Freq)] <- 0
  names(phylo_diver) <- c('trap','date',paste('pd',i,sep='_'))
  
  diversity_data <- merge(diversity_data,phylo_diver,by=c('trap','date'))
}


#plot boxplots
pdf(file = 'boxplots_diversity.pdf',
    width = 4,
    height = 6)
par(mfrow=c(3,1), 
    omi = c(0.5,0,0.2,0), 
    mai = c(0,0.6,0.2,0.25),
    cex.axis = 0.9,
    cex.lab = 0.8)

#plot abundance
boxplot(Freq~trap,
        data = abundance, 
        border = plot_colors, 
        axes = F)
mtext(side = 2, text = 'Abundance', line = 2.5, cex = 0.8)
axis(2)
axis(1,labels = FALSE)
box()

title(main = 'Beetles collected per day', outer = T)
#plot richness
boxplot(richness~trap,
        data = species_richness, 
        border = plot_colors, 
        axes = F)
mtext(side = 2, text = 'Species richness', line = 2.5, cex = 0.8)

axis(2)
axis(1,labels = FALSE)
box()

#plot phylogenetic diversity
boxplot(diversity_data[[1+3]]~trap,
        data = diversity_data, 
        ylim = range(diversity_data[4:length(diversity_data)]),
        lwd=0.15,
        border = plot_colors)
mtext(side = 2, text = 'Phylogenetic diversity', line = 2.5, cex = 0.8)
for (i in 2:100){
  boxplot(diversity_data[[i+3]]~trap,
          data = diversity_data,
          ylim = range(diversity_data[4:length(diversity_data)]),
          lwd=0.15,
          add = T,
          bty='n',
          axes=F,
          border = plot_colors)
}


dev.off()


#do linear model for phylogenetic diversity for 100 trees
#store coefficients and p values
model_results = array(data=NA,dim=c(4,2,100),dimnames = list(c('(Intercept)','trapNa','trapNa_F', 'trapControl'), c('Estimate','t value'),1:100))

for (i in 1:100){
  pd_model = lmer(diversity_data[[i+3]]~trap+(1|date), data = diversity_data)
  model_results[,,i] = summary(pd_model)$coefficients[c('(Intercept)','trapNa','trapNa_F', 'trapControl'),c('Estimate','t value')]
}

#variation in estimated coefficients due to phylogeny
boxplot(Freq~Var1,data=as.data.frame.table(model_results[,'Estimate',]))
#variation in p-value due to phylogeny
boxplot(Freq~Var1,data=as.data.frame.table(model_results[,'t value',]))
##Variation is negligible, will do model with one of the trees

pd_model = lmer(sqrt(diversity_data[[i+3]])~trap+(1|date), data = diversity_data)
summary(pd_model)
car::Anova(pd_model)
plot(fitted(object=pd_model),residuals(pd_model))
qqplot(fitted(pd_model),abundance$Freq)
hist(residuals(pd_model))

cor.test(~richness+pd_1, diversity_data)
pdf('PD_S_correlation.pdf',width = 4,height = 3)
par(omi = c(0,0,0,0), 
    mai = c(0.7,0.7,0.3,0.25),
    cex.axis = 0.7,
    cex.lab = 0.8,
    mgp = c(0,0.7,0))
plot(diversity_data[c('richness','pd_1')],
     col=plot_colors[as.numeric(diversity_data$trap)],
     axes=F,
     xlab='',
     ylab='')
axis(1)
axis(2)
box()
mtext(text = 'Species richness',side = 1,line = 2)
mtext(text = 'Phylogenetic diversity',side = 2,line = 2)
abline(0,lm(pd_1~0+richness,diversity_data)$coefficients)
dev.off()
###2 - rarefaction curves
# make community matrix
count_matrix <- xtabs(~trap+paste(higher_taxon,species,sep="_"),dados) #create community matrix
total_counts <- tapply(dados$species,dados$trap,length) #count number of specimens
rarefy_results <- array(data = NA, #create an array to save results
                        dim = c(4,2,max(total_counts)),
                        dimnames = list(traps = levels(dados$trap), stats = c('S','se'), N=1:max(total_counts))
                        )
for (i in 1:max(total_counts)){ #do rarefaction
  rarefy_results[total_counts>=i,,i] = t(rarefy(count_matrix[total_counts>=i,], sample = c(i), se = T))
}

pdf('rarefaction.pdf',width = 4,height = 3)
par(omi = c(0,0,0,0), 
    mai = c(0.7,0.7,0.3,0.25),
    cex.axis = 0.7,
    cex.lab = 0.8,
    mgp = c(0,0.7,0))

plot_lines = function(xmax,ymax){
  plot(NULL, #start plot
       xlim = c(0,xmax),#range(as.numeric(dimnames(rarefy_results)[[3]])),
       ylim = c(min(rarefy_results[,1,]-rarefy_results[,2,],na.rm = T),ymax),
       xlab = '',
       ylab = '')
for (i in 1:4){ #plot lines
  lines(as.numeric(dimnames(rarefy_results)[[3]]), #plot estimated S
        rarefy_results[i,1,], col = plot_colors[i])
  lines(as.numeric(dimnames(rarefy_results)[[3]]), #plot S + se
        rarefy_results[i,1,] + rarefy_results[i,2,],
        lty = 2, col = plot_colors[i])
  lines(as.numeric(dimnames(rarefy_results)[[3]]), #plot S -se
        rarefy_results[i,1,] - rarefy_results[i,2,],
        lty = 2, col = plot_colors[i])
}}
plot_lines(xmax = 200,
           ymax = 120)

title(main = 'Rarefaction curves')
mtext(side = 2, text = 'Species richness', line = 2, cex = 0.8)
mtext(side = 1, text = 'Number of individuals', line = 2, cex = 0.8)


text(x = c(150,150,125,100),
     y = c(87,60,42,23),
     labels = levels(dados$trap),
     col = plot_colors)
subplot( 
  plot_lines(xmax = max(as.integer(dimnames(rarefy_results)[[3]])),
             ymax = max(rarefy_results[,1,]+rarefy_results[,2,], na.rm = T)),
  x=c(10,75),
  y=c(66,115))
dev.off()

#proportion of singletons
sum(xtabs(~species, dados) == 1)/length(levels(dados$species))
#over 50%

#####phylogenetic clustering
#1 - modify ses.mpd function to sample a different tree in each randomization trial
new.ses.mpd = function (samp, trees, null.model = c("taxa.labels", "richness", 
                                    "frequency", "sample.pool", "phylogeny.pool", "independentswap", 
                                    "trialswap"), abundance.weighted = FALSE, runs = 999, iterations = 1000) 
{
  dis <- function(tr=trees){as.matrix(cophenetic(tr[[sample(1:length(tr),1)]]))} #modified part, this will sample a new tree in each replication round
  mpd.obs <- apply(replicate(runs,mpd(count_matrix,dis(),abundance.weighted = abundance.weighted)),1,mean)
  null.model <- match.arg(null.model)
  mpd.rand <- switch(null.model, trialswap = t(replicate(runs,mpd(randomizeMatrix(samp, null.model = "trialswap", iterations), dis(), abundance.weighted))))
  mpd.rand.mean <- apply(X = mpd.rand, MARGIN = 2, FUN = mean, 
                         na.rm = TRUE)
  mpd.rand.sd <- apply(X = mpd.rand, MARGIN = 2, FUN = sd, 
                       na.rm = TRUE)
  mpd.obs.z <- (mpd.obs - mpd.rand.mean)/mpd.rand.sd
  mpd.obs.rank <- apply(X = rbind(mpd.obs, mpd.rand), MARGIN = 2, 
                        FUN = rank)[1, ]
  mpd.obs.rank <- ifelse(is.na(mpd.rand.mean), NA, mpd.obs.rank)
  data.frame(ntaxa = specnumber(samp), mpd.obs, mpd.rand.mean, 
             mpd.rand.sd, mpd.obs.rank, mpd.obs.z, mpd.obs.p = mpd.obs.rank/(runs + 
                                                                               1), runs = runs, row.names = row.names(samp))
}

weighted = ses.mpd(count_matrix,cophenetic(trees[[1]]), null.model = 'trialswap',abundance.weighted = T, runs=10000, iterations=100000)
weighted2 = new.ses.mpd(count_matrix,trees,null.model = 'trialswap',abundance.weighted = T, runs=10000, iterations=100000)
unweighted = ses.mpd(count_matrix,cophenetic(trees[[1]]), null.model = 'trialswap',abundance.weighted = F, runs=10000, iterations=100000)
unweighted2 = new.ses.mpd(count_matrix,trees, null.model = 'trialswap',abundance.weighted = F, runs=10000, iterations=100000)

