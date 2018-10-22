################################################################################
## Script Name: Network analysis with DNA-hyperauthor data
## Author: Ly Dinh
## Date Created: April 24, 2018
################################################################################

## Data Sources
## Network data: Co-Authorship of DNA Chromosome 1 paper 
## Nodes: 1901 (original n = 166)
## Edges: if 2 nodes have co-authored at least 5 papers together (e=14088)
## Attribute Data: Gender, affiliation

################################################################################
## Data Pre-Processing
## Set working directory
## Install appropropriate network & statistics packages
## Converting data files from .txt to graph data frame
## Ensure no nodes or edges are missing

################################################################################
## Data analysis
## Run basic metrics & centrality measures at network-level and other levels of analysis if need be
## Plot degree distributions, and correlation among centrality measures 

################################################################################

#clear all in the 'global environment' 
rm(list=ls())

#Set working directory
setwd("/Users/lydinh/Documents/MAI LY/U of Illinois/Spring 2018/Jessica & Ly")

install.packages("igraph")
install.packages("plyr")
install.packages("calibrate")
install.packages("poweRlaw")

install.packages("bibliometrix")
library(bibliometrix)

library(igraph)
library(plyr) ##Tools for Splitting, Applying and Combining Data
library(calibrate) ##also for formatting 
library(poweRlaw)



dat=read.csv(file.choose('Edgelist_n166_e14088.csv'),header=TRUE) # choose an edgelist in .csv file format
el=as.matrix(dat) # coerces the data into a two-column matrix format that igraph likes
el[,1]=as.character(el[,1])
el[,2]=as.character(el[,2])
hypernet=graph.edgelist(el,directed=FALSE) # turns the edgelist into a 'graph object'

V(hypernet)$name

## remove self loops
is_simple(simplify(hypernet, remove.loops=TRUE))

#check header to see dataframe
head(hypernet)

# Counts for number of nodes
vcount(hypernet)
# Counts for number of edges
ecount(hypernet)

# Plot the network using the plot() function, of graph dataframe
plot(hypernet, vertex.size=1)

# interactive plot function!
tkplot(hypernet)
## choose layouts: Fruchterman-Reingold or Kamada-Kawai = all force-directed graphs
#spring forces proportional to the graph theoretic distances = similar attract, differences separated


##### Basic analysis
# Inspect vertex degree and use the degree for vertex size
hist(degree(hypernet))

# retrieve vertex ids:
# vertex names are usually stored in a vertex attribute named name in igraph. use V(g)$name to retrieve the names of all the vertices.
match("Seymour, Lesley K.", V(hypernet)$name)
match("Nelson, Sarah C.", V(hypernet)$name)

## Get vertex ids of ALL actors:
V(hypernet)$name

## Number of paths between one vertex () and another vertex ()
edge.disjoint.paths(hypernet, 1, 1918)


# I can now test the various centrality measures to have a feel of the whole network centralization
# betweenness centrality
betweenness(hypernet)
#closeness centrality
closeness(hypernet)
#local clustering coefficient
transitivity(hypernet, type="average")
#eigenvector centrality
evcent(hypernet)$vector

#summary for each centrality
summary(betweenness(hypernet))
summary(closeness(hypernet))


## SORTING function highest , and view top 10 
## top ten BETWEENNESS
sort(betweenness(hypernet), decreasing = TRUE)[1:10]
## top ten EIGENVECTOR
sort(evcent(hypernet)$vector, decreasing = TRUE)[1:10]
## top ten DEGREE
sort(degree(hypernet), decreasing = TRUE)[1:10]

## more types of analysis we can compare
graph.density(hypernet)
triangles(hypernet)
diameter(hypernet)


############ SCALE-FREE WHOLE NETWORK

## average path length  ; should be LOW 
average.path.length(hypernet, directed=FALSE, unconnected=TRUE)

## clustering coefficient; should be HIGH 
transitivity(hypernet, type="average")


deg = degree(hypernet, mode = "all")
hist(deg)

degreedist = degree.distribution(hypernet, mode = "all", cumulative = FALSE)

# Plot degree distribution
plot(degreedist)


## Fit POWER LAW onto hypernet
fit_power_law = function(hypernet) {
  # calculate degree
  d = degree(hypernet, mode = "all")
  dd = degree.distribution(hypernet, mode = "all", cumulative = FALSE)
  degree = 1:max(d)
  probability = dd[-1]
  # delete blank values
  nonzero.position = which(probability != 0)
  probability = probability[nonzero.position]
  degree = degree[nonzero.position]
  reg = lm(log(probability) ~ log(degree))
  cozf = coef(reg)
  power.law.fit = function(x) exp(cozf[[1]] + cozf[[2]] * log(x))
  alpha = -cozf[[2]]
  R.square = summary(reg)$r.squared
  print(paste("Alpha =", round(alpha, 3)))
  print(paste("R square =", round(R.square, 3)))
  # plot
  plot(probability ~ degree, log = "xy", xlab = "Degree (log)", ylab = "Probability (log)", 
       col = 1, main = "Chroms1 Degree Distribution")
  curve(power.law.fit, col = "red", add = T, n = length(d))
}

fit = fit_power_law(hypernet)
## Alpha = 0.582
## R^2 = 0.431

# Fit power-law

fit = power.law.fit(degreedist, implementation = 'plfit')
fit

mean(degree(hypernet))

### Restrict to only certain vertices in subgraph
V(hypernet)[name == "Gregory, Simon G."]

sub1 = V(hypernet)$name == c("Gregory, Simon G.","Barlow, Karen F.","Mclay, Kirsten E.","Kaul, R.","Swarbreck, D.","Dunham, A.","Scott, C.E.","How K.L.","Woodfine, K.","Spencer, C.C.A.","Jones, M.C.","Gillson, C.","Searl S.","Zhou, Y.","Kokocinski, F.","McDonald L.","Evans, R.","Phillips, K.","Atkinson, A.","Cooper, R.","Jones, C.","Hall, R.E.","Andrews, T.D.","Lloyd, C.","Ainscough, R.","Almeida, J.P.","Ambrose, K.D.","Anderson, F.","Andrew, R.W.","Ashwell, R.I.S.","Aubin, K.","Babbage, A.K.","Bagguley, C.L.","Bailey, J.","Beasley, H.","Bethel, G.","Bird, C.P.","Bray-Allen, S.","Brown, J.Y.","Brown, A.J.","Buckley, D.","Burton, J.","Bye, J.","Carder, C.","Chapman, J.C.","Clark, S.Y.","Clarke, G.","Clee, C.","Cobley, V.","Collier, R.E.","Corby, N.","Covill G.J.","Davies, J.","Deadman, R.","Dunn, M.","Earthrowl, M.","Ellington, A.G.","Errington, H.","Frankish, A.","Frankland, J.","French, L.","Garner, P.","Garnett, J.","Gay, L.","Ghori, M.R.J.","Gibson, R.","Gilby, L.M.","Gillett, W.","Glithero, R.J.","Grafham, D.V.","Griffiths, C.","Griffiths-Jones, S.","Grocock, R.","Hammond, S.","Harrison, E.S.I.","Hart, E.","Haugen, E.","Heath, P.D.","Holmes, S.","Holt, K.","Howden, P.J.","Hunt, A.R.","Hunt, S.E.","Hunter, G.","Isherwood, J.","James, R.","Johnson, C.","Johnson, D.","Joy, A.","Kay, M.","Kershaw, J.K.","Kibukaw M.","Kimberley, A.M.","King, A.","Knights, A.J.","Lad, H.","Laird, G.","Lawlor, S.","Leongamornlert, D.A.","Lloyd, D.M.","Loveland, J.","Lovell, J.","Lush, M.J.","Lyne, R.","Martin, S.","Mashreghi-Mohammadi, M.","Matthews, L.","Matthews, N.S.W.","McLaren, S.","Milne, S.","Mistry, S.","Moore, M.J.F.","Nickerson, T.","O'Dell, C.N.","Oliver, K.","Palmeiri, A.","Palmer, S.A.","Parker, A.","Patel, D.","Pearce, A.V.","Peck, A.I.","Pelan, S.","Phelps, K.","Phillimore, B.J.","Plumb, R.","Rajan, J.","Raymond, C.","Rouse, G.","Saenphimmachak, C.","Sehra, H.K.","Sheridan, E.","Shownkeen, R.","Sims, S.","Skuce, C.D.","Smith, M.","Steward, C.","Subramanian, S.","Sycamore, N.","Tracey, A.","Tromans, A.","Van Helmond, Z.","Wall, M.","Wallis, J.M.","White, S.","Whitehead, S.L.","Wilkinson, J.E.","Willey, D.L.","Williams, H.","Wilming, L.","Wray, P.W.","Wu, Z.","Coulson, A.","Vaudin, M.","Sulston, J.E.","Durbin, R.","Hubbard, T.","Wooster, R.","Dunham, I.","Carter, N.P.","McVean, G.","Ross, M.T.","Harrow, J.","Olson, M.V.","Beck, S.","Rogers, J.","Bentley, D.R.")


## Induced Subgraphs (means: containing only the specified vertices and all the edges among them.)
## create subset of ids 
subid = as.numeric(V(hypernet)[(c("Gregory, Simon G.","Barlow, Karen F.","McLay, Kristen E.","Kaul, Rajinder K.","Swarbreck, David","Dunham, Andrew","Scott, Carol E.","Howe, Kevin Lee","Woodfine, Kathryn","Spencer, Chris CA","Jones, Matthew C.","Gillson, Christopher","Searle, Stephen M.J.","Zhou, Yang","Kokocinski, Felix","McDonald L.","Evans, Richard S.","Phillips, Kimbly","Atkinson, Alexandre L.","Cooper, Rachel A.","Jones, Carol A.","Hall, Rebekah E.","Andrews, T. Daniel","Lloyd, Christine R.","Ainscough, Rachael","Almeida, Jeff P.","Ambrose, Kerrie D.","Anderson, F.","Andrew, R.W.","Ashwell, Robert I S","Aubin, K.","Babbage, Anne K.","Bagguley, Claire L.","Bailey, Jonathon","Beasley, Helen","Bethel, Graeme","Bird, Christine P.","Bray-Allen, Sarah P.","Brown, Jacqueline Y.","Brown, Alex J.","Buckley, Danielle G.","Burton, John L.","Bye, Jacqueline","Carder, Carol","Chapman, Joanne C.","Clark, Sue Y.","Clarke, Graham D.Paul","Clee, Chris M.","Cobley, Victoria E.","Collier, R. E.","Corby, Nicole R.","Coville, G. J.","Davies, Joy R.","Deadman, Rebecca","Dunn, Matthew","Earthrowl, Mark E.","Ellington, A. G.","Errington, H.","Frankish, Adam G.","Frankland, John A.","French, Lisa","Garner, Patrick","Garnett, Jane","Gay, L.","Ghori, Jilur","Gibson, Richard S.","Gilby, L.M.","Gillett, Will G.","Glithero, Rebecca J.","Grafham, Darren V.","Griffiths, Coline","Griffiths-Jones, Sam","Grocock, Russell James","Hammond, Sian Austin","Harrison, Elliot S I","Hart, Elizabeth A.","Haugen, Eric","Heath, Paul D.","Holmes, Sharon","Holt, Karen","Howden, Philip J.","Hunt, Adrienne R.","Hunt, Sarah E.","Hunter, Giselle","Isherwood, Judith","James, Rosalina","Johnson, Christopher M.","Johnson, David C.","Joy, Ann A.","Kay, Mike","Kershaw, Joanne K.","Kibukawa, Miho","Kimberley, Andrew M.","King, Andrew","Knights, Andrew J.","Lad, H.","Laird, Gavin K.","Lawlor, Stephanie","Leongamornlert, Daniel A.","Lloyd, David M.L.","Loveland, Jane E.","Lovell, Jamieson D.","Lush, Michael J.","Lyne, Rachel","Martin, Sancha L.","Mashreghi-Mohammadi, Maryam","Matthews, Lucy H.","Matthews, Nicholas S.W.","McLaren, S.","Milne, Sarah A.","Mistry, Shailesh L.","Moore, Melissa J.F.","Nickerson, T.","O'Dell, Christopher","Oliver, Karen","Palmeiri, Anthony","Palmer, Sophie A.","Parker, Anne","Patel, Dina","Pearce, Alex V.","Peck, A. I.","Pelan, Sarah E.","Phelps, Karen A.","Phillimore, Benjamin J.C.T.","Plumb, Robert W.","Rajan, Jeena Joseph Sahaya","Raymond, Christopher K.","Rouse, Greg","Saenphimmachak, Channakhone","Sehra, Harminder K.","Sheridan, Elizabeth","Shownkeen, Ratna","Sims, Sarah K.","Skuce, Carl D.","Smith, Michelle L.","Steward, Charles A.","Subramanian, Sandhya","Sycamore, Neil","Tracey, Alan","Tromans, Anthony C.","Van Helmond, Z.","Wall, Melanie","Wallis, Justine M.","White, Simon S.","Whitehead, Siobhan L.","Wilkinson, Jane","Willey, David L.","Williams, Hawys","Wilming, Laurens G.","Wray, Paul W.","Wu, Zaining","Coulson, Alan R.","Vaudin, Mark D.","Sulston, John E.","Durbin, Richard M.","Hubbard, Tim J.P.","Wooster, Richard W.","Dunham, Ian","Carter, Nigel P.","McVean, Gil","Ross, Mark T.","Harrow, Jennifer L.","Olson, Maynard V.","Beck, Stephan R.","Rogers, Jane","Bentley, David R."))])

subg <- induced_subgraph(hypernet, subid)
plot(subg, vertex.label=subg$name, layout=layout_with_fr, vertex.label.color="black", vertex.label.cex=0.8, vertex.label.dist=2, edge.curved=0.2)




## SORTING function highest , and view top 10 
## top ten BETWEENNESS
sort(betweenness(subg), decreasing = TRUE)[1:10]
## top ten EIGENVECTOR
sort(evcent(subg)$vector, decreasing = TRUE)[1:10]
## top ten DEGREE
sort(degree(subg), decreasing = TRUE)[1:10]

ecount(subg)
vcount(subg)

save(subg, file="subgraph.csv")




## Community Detection 
# Grivan-Newman algorithm
# 1st we calculate the edge betweenness, merges, etc...
ebc <- edge.betweenness.community(hypernet, directed=F)

plot(ebc, hypernet)


# Now we have the merges/splits and we need to calculate the modularity
# for each merge for this we'll use a function that for each edge
# removed will create a second graph, check for its membership and use
# that membership to calculate the modularity
mods <- sapply(0:ecount(hypernet), function(i){
  g2 <- delete.edges(hypernet, ebc$removed.edges[seq(length=i)])
  cl <- clusters(g2)$membership
  # March 13, 2014 - compute modularity on the original graph g 
  modularity(hypernet,cl)
})

# we can now plot all modularities
plot(mods, pch=20)




######### GRAPH MODELS
## ERGMS??? 

## generate random graphs according to the G(n,m) Erdos-Renyi model. n = number of vertices , m = number of edges
randomgraph = sample_gnm(n = vcount(hypernet), m = ecount(hypernet), directed=FALSE)
degree_distribution(randomgraph)

diameter(randomgraph)
graph.density(randomgraph)
triangles(randomgraph)
summary(betweenness(randomgraph))

plot(degree_distribution(randomgraph))


## generate scale-free graphs according to the Barabasi-Albert model. n = number of vertices , m = number of edges added on each time step
pa = sample_pa(vcount(hypernet), m=2, directed = FALSE)

plot(degree_distribution(pa))

diameter(pa)
graph.density(pa)
triangles(pa)
summary(betweenness(pa))


##  PLOTS & visualization
# I can also plot the degree distribution to see how spread out are the differences in the number of connections each node has
plot(degree.distribution(hypernet), xlab="node degree")


# Basic visualization with plot.igraph, a step from simple plot function
plot.igraph(hypernet,layout=layout.fruchterman.reingold, vertex.color='red',vertex.size=degree(hypernet), vertex.label=V(hypernet)$name)

# or, if we want to remove the vertex labels --- vertex.label=NA
plot.igraph(hypernet,vertex.label=NA, layout=layout.fruchterman.reingold, vertex.color='red',vertex.size=degree(lesmis_gdf), vertex.label=lesmis)

## try different layouts  layout=layout.fruchterman.reingold
plot.igraph(hypernet, layout=layout.circle, vertex.color='red',vertex.size=degree(hypernet), vertex.label=V(hypernet)$name)


# One nice thing for 'plot' is that I can plot centrality measures against each other to see correlations among various metrics, and identify overlapping central actors.

## Why plot eigenvector with betweenness 
## High betweenness and low eigenvector centrality may be an important gatekeeper to a central actor.
# Low betweeness and high eigenvector centrality may have unique access to central actors.

# plot eigenvector with betweenness 
plot(evcent(hypernet)$vector, betweenness(hypernet), xlab='eigenvector', ylab='betweenness', xlim=c(0,1), cex=1.5)

textxy(evcent(hypernet)$vector, betweenness(hypernet), hypernet$From, pos=2, cex=0)
