syllabus <- read.csv('https://raw.githubusercontent.com/PeterCC36/syllables/main/20-7-nw.csv')

library(vegan)

pam.clustering = function(x,k) { # x is a distance matrix and k the number of clusters
  require(cluster)
  cluster = as.vector(pam(as.dist(x), k, diss = TRUE)$clustering)
  return(cluster)
}

#box plot
syllable = read.csv('https://raw.githubusercontent.com/PeterCC36/syllables/main/20-7-nw.csv')
song = read.csv('https://raw.githubusercontent.com/PeterCC36/syllables/main/20-7-location-song.csv') 

min.freq = syllable[, c(10, 5)]
max.freq = syllable[, c(10, 6)]
freq.range = syllable[, c(10, 7)]
duration = syllable[, c(10, 8)]

nsyllable = song[, c(6, 2)]
length = song[, c(6, 3)]
rate = song[, c(6, 4)]

par(mfrow = c(2, 4))
boxplot(min.freq$min.freq ~ min.freq$species, las = 1, xlab = 'Species', ylab = 'Frequency(Hz)'
        , main = 'Minimum Frequency', col = c('grey40', 'grey60', 'grey90'))
boxplot(max.freq$max.freq ~ max.freq$species, las = 1, xlab = 'Species', ylab = 'Frequency(Hz)', main = 'Maximum Frequency', col = c('grey40', 'grey60', 'grey90'))
boxplot(freq.range$freq.range ~ freq.range$species, las = 1, xlab = 'Species', ylab = 'Frequency(Hz)', main = 'Frequency Range', col = c('grey40', 'grey60', 'grey90'))
boxplot(duration$duration ~ duration$species, las = 1, xlab = 'Species', ylab = 'Duration(s)', main = 'Syllable Duration', col = c('grey40', 'grey60', 'grey90'))
boxplot(nsyllable$Number.of.syllable ~ nsyllable$Species, las = 1, xlab = 'Species', ylab = 'Number', main = 'Number of Syllables', col = c('grey40', 'grey60', 'grey90'))
boxplot(length$Length ~ length$Species, las = 1, xlab = 'Species', ylab = 'Time(seconds)', main = 'Song Duration', col = c('grey40', 'grey60', 'grey90'))
boxplot(rate$Note.rate ~ rate$Species, las = 1, xlab = 'Species', ylab = 'syllables/second', main = 'Syllable rate', col = c('grey40', 'grey60', 'grey90'))

#ANOVA
summary(aov(min.freq$min.freq ~ min.freq$species))
summary(aov(max.freq$max.freq ~ max.freq$species))
summary(aov(freq.range$freq.range ~ freq.range$species))
summary(aov(duration$duration ~ duration$species))
summary(aov(nsyllable$Number.of.syllable ~ nsyllable$Species))
summary(aov(length$Length ~ length$Species))
summary(aov(rate$Note.rate ~ rate$Species))

#Tukey
par(mfrow = c(2, 4))
plot(TukeyHSD(aov(min.freq$min.freq ~ min.freq$species)))
plot(TukeyHSD(aov(max.freq$max.freq ~ max.freq$species)))
plot(TukeyHSD(aov(freq.range$freq.range ~ freq.range$species)))
plot(TukeyHSD(aov(duration$duration ~ duration$species)))
plot(TukeyHSD(aov(nsyllable$Number.of.syllable ~ nsyllable$Species)))
plot(TukeyHSD(aov(length$Length ~ length$Species)))
plot(TukeyHSD(aov(rate$Note.rate ~ rate$Species)))

#draw a cluster analysis dendrogram
types <- syllabus[, 5:8]
eu.dist <- dist(types)
cluster <- hclust(eu.dist)
par(mfrow = c(1, 1))
plot(cluster, cex=0.5)
rect.hclust(cluster, 8)

#record the classification value of k-medoids
k8 <- pam.clustering(eu.dist, 8)
syllabus$cluster8 <- k8

#NMDS 
NMDS <- metaMDS(types, group = as.factor(syllabus$cluster8), distance = 'euclidian')
ordiplot(NMDS, display = 'sites', type = 'n')
points (NMDS, pch = 20, col = k8, cex = 2)
legend('right', c('1', '2', '3', '4', '5', '6', '7', '8'), fill = c(1:8))

#PCoA
pcoa <- cmdscale (eu.dist, eig = TRUE)
ordiplot(pcoa, display = 'sites', type = 'n')
points (pcoa$points, pch = 20, col = k8, cex = 2)
legend('right', c('1', '2', '3', '4', '5', '6', '7', '8'), fill = c(1:8))

#CA
data <- syllabus
location <- as.data.frame(table(data$File, data$location))
location <- location[location$Freq > 0, ]
data <- table(as.data.frame(data)$cluster8, data$File)
class(data) <- "matrix"

ca <- cca(t(data))
scores(ca)$site# position of each point

loc.comb <- factor(location$Var2,labels = c("H","B","H","W","W","B","H")) #lable the species, H = hybrid, B = Styan's, W = light-vented

par(mar = c(5.1, 4.1, 4.1, 10.1))
plot(scores(ca)$site,pch=c(19,19,1)[loc.comb], col=c("#ffa801","#000000","#000000")[loc.comb], xlab=paste("CA1 (", round((prop.table(ca$CA$eig)*100)[1], 2), "%)"), ylab=paste("CA2 (", round((prop.table(ca$CA$eig)*100)[2],2), "%)"))
legend("topright", legend = levels(loc.comb), bty="n", pch=c(19,19,1), col=c("#ffa801","#000000","#000000"), xpd=T, inset = c(-0.15, 0))

ordiellipse(ca, loc.comb, lty=c(1, 1, 4), col = c("#b15928","#000000","#000000"), label = T, draw = "line", kind = "sd", lwd=1.5)
#Functions to add convex hulls, ¡§spider¡¨ graphs, ellipses or cluster dendrogram to ordination diagrams. The ordination diagrams can be produced by vegan plot.cca, plot.decorana or ordiplot.

#chi-square test
chisquare <- read.csv('https://raw.githubusercontent.com/PeterCC36/syllables/main/qisquare.csv')
chisq = chisq.test(chisquare)#p ~ 0, there are significant differences among three populations.

mosaicplot(chisquare, las=1, color = c('white', 'grey80', 'black'), main = '')#draw a mosaic plot for visualization
