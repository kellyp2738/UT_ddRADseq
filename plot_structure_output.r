## Attepts to visualize structure output:

# This is the structure output that denotes %membership in each of the 15 possible host clusters for an individual tick
# If each individual host has a unique set of tick genotypes, we expect to see a correlation between cluster representation and host identity.
# Alternatively, if the tick population is admixed and host doesn't predict tick genotypes, then we should not see a correlation between cluster membership and host identity
# It's important to remember also that I told Structure there were 15 populations, and structure inferred membership in those populations accordingly.
# Just because there are different color clusters that emerge does not mean that there is anything meaningful going on.

clusters<-read.table('~/Desktop/distruct1.1/second_try_tick.indivq')
head(clusters)

# each row is an individual tick, and columns 6-20 represent host cluster membership percentages
cluster.totals<-rowSums(clusters[,6:20]) # everything more or less adds up to 1

?order
cluster.ordered<-clusters[order(clusters[,6],
                                clusters[,7],
                                clusters[,8],
                                clusters[,9],
                                clusters[,10],
                                clusters[,11],
                                clusters[,12],
                                clusters[,13],
                                clusters[,14],
                                clusters[,15],
                                clusters[,16],
                                clusters[,17],
                                clusters[,18],
                                clusters[,19]),]
head(cluster.ordered)
?barplot
barplot(t(as.matrix(cluster.ordered[,6:20])))

other.clusters<-read.table('~/Desktop/distruct1.1/first_try_tick.indivq')
other.ordered<-other.clusters[order(other.clusters[,6]),]
barplot(t(as.matrix(other.ordered[,6:7])))

other.ordered[1:5,]
IDs<-read.table('~/Desktop/distruct1.1/first_try_tick.ids.csv')
high.cluster.2<-IDs[other.ordered[1:5,1],]
high.cluster.2[,2]
