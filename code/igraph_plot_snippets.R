
plot_TF <- graph_from_edgelist(as.matrix(find.effects_TF[,.(geneA, geneB)]),directed=TRUE)
test <- find.effects_TF[,.(geneA, geneB)]

#test.plot <- test[test$geneA %in% test$geneA[duplicated(test$geneA)],]
test.plot <- test[, if (.N>1) .SD, by=geneA]


# l <- layout_with_fr(test.plot) 
# plot(net.bg, layout=l)

plot_TF <- graph_from_data_frame(find.effects_TF[1:100],directed=TRUE)


#V(plot_TF)$label <- NA

# E(allcases.graph)$color[E(allcases.graph)$case == 1] <- 'red'    
# E(allcases.graph)$color[E(allcases.graph)$case == 2] <- 'blue'  
# E(allcases.graph)$color[E(allcases.graph)$case == 4] <- 'red'  
V(plot_TF)$size <- 5

plot(plot_TF, layout=layout_as_tree, main="T & F", edge.arrow.size=.5, vertex.label=NA, edge.curved=.1)
plot(plot_TF, layout=layout_with_lgl, main="T & F", edge.arrow.size=.4, vertex.label=NA, edge.curved=.1)
plot(plot_TF, layout=layout_on_sphere, main="T & F", edge.arrow.size=.4, vertex.label=NA, edge.curved=.1)

toplot <- graph_from_edgelist(as.matrix(test))
V(toplot)$size <- 5

plot(toplot, layout=layout_with_fr, vertex.label=NA, edge.arrow.size=.4)




E(toplot)[V(toplot)[genes.A=="YAL003W"] %--% V(g)[genes.B=="YAL033W"]]
E(toplot) [ 1:3 %--% 2:6 ]



genes <- unique(find.effects_TF$geneA)
gene <- genes[1]
res.table <- data.table(matrix(ncol = 6, nrow=0))
colnames(res.table) <- colnames(find.effects_TF)
# for (gene in genes){
#   temp <- rbindlist(find.effects_TF[])
# }
temp <- rbind(res.table, find.effects_TF[geneA==gene | geneB==gene])
if (gene == temp$geneA){
  gene <- temp$geneB
  else if (gene == temp$geneB){
    gene <- temp$geneA
  }
}


########



library(GGally)
library(network)
library(sna)
library(ggplot2)


grph <- graph_from_data_frame(find.effects_TF[,.(geneA, geneB)],directed=TRUE)

ggnet2(simplify(grph))





#####

effects <- find.effects

# Define cases
# one gene affects the other but we don't know the action of the other
case1 <- effects[(is.na(effects$`A->B`) & effects$`B->A` == TRUE) | (is.na(effects$`B->A`) & effects$`A->B` == TRUE)]
# one gene affects the other but is not affected by the other gene
case2 <- effects[(effects$`A->B` == TRUE & effects$`B->A` == FALSE) | (effects$`A->B` == FALSE & effects$`B->A` == TRUE)]
# we don't know if the genes affect each other
case3 <- effects[is.na(effects$`A->B`) & is.na(effects$`B->A`)]
# a gene affects the other but it's also affected by that other gene
case4 <- effects[(effects$`A->B` == TRUE & effects$`B->A` == TRUE)]
# no gene affects the other
case5 <- effects[(effects$`A->B` == FALSE & effects$`B->A` == FALSE)]




case1.A_to_B <- case1[case1$`A->B`==TRUE & is.na(case1$`B->A`)]
case1.A_to_B <- case1.A_to_B[,c("geneA","geneB", "A->B")]

case1.B_to_A <- case1[case1$`B->A`==TRUE & is.na(case1$`A->B`)]
case1.B_to_A <- case1.B_to_A[,c("geneB", "geneA", "B->A")]

if (nrow(case1.A_to_B) > 0 & nrow(case1.B_to_A) > 0){
  case1.links <- rbind(case1.A_to_B, case1.B_to_A, use.names=FALSE)
  
} else if (nrow(case1.A_to_B) == 0 & nrow(case1.B_to_A) > 0){
  case1.links <- case1.B_to_A
  
} else if (nrow(case1.A_to_B) > 0 & nrow(case1.B_to_A) == 0){
  case1.links <- case1.A_to_B
  
} else if (nrow(case1.A_to_B) == 0 & nrow(case1.B_to_A) == 0){
  case1.links <- case1.A_to_B
}

colnames(case1.links) <- c("from", "to", "type")
case1.links$case <- 1
case1.links <- unique(case1.links)

case2.A_to_B <- case2[case2$`A->B`==TRUE & case2$`B->A`==FALSE]
case2.A_to_B <- case2.A_to_B[,c("geneA","geneB", "A->B")]

case2.B_to_A <- case2[case2$`B->A`==TRUE & case2$`A->B`==FALSE]
case2.B_to_A <- case2.B_to_A[,c("geneB", "geneA", "B->A")]

if (nrow(case2.A_to_B) > 0 & nrow(case2.B_to_A) > 0){
  case2.links <- rbind(case2.A_to_B, case2.B_to_A, use.names=FALSE)
  
} else if (nrow(case2.A_to_B) == 0 & nrow(case2.B_to_A) > 0){
  case2.links <- case2.B_to_A
  
} else if (nrow(case2.A_to_B) > 0 & nrow(case2.B_to_A) == 0){
  case2.links <- case2.A_to_B
  
} else if (nrow(case2.A_to_B) == 0 & nrow(case2.B_to_A) == 0){
  case2.links <- case2.A_to_B
}

colnames(case2.links) <- c("from", "to", "type")
case2.links$case <- 2
case2.links <- unique(case2.links)

case3.links <- case3[,c("geneA", "geneB")]
case3.links$type <- NA
case3.links <- unique(case3.links)
colnames(case3.links) <- c("from", "to", "type")
case3.links$case <- 3

case4.A_to_B <- case4[case4$`A->B`==TRUE & case4$`A->B`==TRUE]
case4.A_to_B <- case4.A_to_B[,c("geneA","geneB", "A->B")]

case4.B_to_A <- case4[case4$`B->A`==TRUE & case4$`B->A`==TRUE]
case4.B_to_A <- case4.B_to_A[,c("geneB", "geneA", "B->A")]

if (nrow(case4.A_to_B) > 0 & nrow(case4.B_to_A) > 0){
  case4.links <- rbind(case4.A_to_B, case4.B_to_A, use.names=FALSE)
  
} else if (nrow(case4.A_to_B) == 0 & nrow(case4.B_to_A) > 0){
  case4.links <- case4.B_to_A
  
} else if (nrow(case4.A_to_B) > 0 & nrow(case4.B_to_A) == 0){
  case4.links <- case4.A_to_B
  
} else if (nrow(case4.A_to_B) == 0 & nrow(case4.B_to_A) == 0){
  case4.links <- case4.A_to_B
}

colnames(case4.links) <- c("from", "to", "type")
case4.links$case <- 4
case4.links <- unique(case4.links)


case5.A_to_B <- case5[case5$`A->B`==FALSE & case5$`A->B`==FALSE]
case5.A_to_B <- case5.A_to_B[,c("geneA","geneB", "A->B")]

case5.B_to_A <- case5[case5$`B->A`==FALSE & case5$`B->A`==FALSE]
case5.B_to_A <- case5.B_to_A[,c("geneB", "geneA", "B->A")]

if (nrow(case5.A_to_B) > 0 & nrow(case5.B_to_A) > 0){
  case5.links <- rbind(case5.A_to_B, case5.B_to_A, use.names=FALSE)
  
} else if (nrow(case5.A_to_B) == 0 & nrow(case5.B_to_A) > 0){
  case5.links <- case5.B_to_A
  
} else if (nrow(case5.A_to_B) > 0 & nrow(case5.B_to_A) == 0){
  case5.links <- case5.A_to_B
  
} else if (nrow(case5.A_to_B) == 0 & nrow(case5.B_to_A) == 0){
  case5.links <- case5.A_to_B
}

colnames(case5.links) <- c("from", "to", "type")
case5.links$case <- 5
case5.links <- unique(case5.links)



allcases.links <- do.call("rbind", list(case1.links, case2.links, case3.links, case4.links, case5.links))




#### Case2 plot
allcases.graph <- graph_from_data_frame(allcases.links[allcases.links$case == 2],directed=TRUE)

E(allcases.graph)$color <- as.factor(E(allcases.graph)$case)

V(allcases.graph)$label <- NA

E(allcases.graph)$color[E(allcases.graph)$case == 1] <- 'red'    
E(allcases.graph)$color[E(allcases.graph)$case == 2] <- 'blue'  
E(allcases.graph)$color[E(allcases.graph)$case == 4] <- 'red'  

# V(allcases.graph)$label.cex = 1.5

# plot(allcases.graph, layout = layout_with_fr,vertex.label.degree=0)
# plot(allcases.graph, layout = layout_with_fr, vertex.label.dist=2)
plot(allcases.graph, layout = layout_as_tree, vertex.label.dist=2, main="Case2 - TRUE and FALSE")
plot(allcases.graph, layout=layout_with_graphopt, vertex.label.dist=2, main="Case2 - TRUE and FALSE")








# gr.meta2 <- graph_from_edgelist(as.matrix(lc.edges[, .(node1, node2)]), directed=T)
# plot(gr.meta2, vertex.label=NA, vertex.color=lc.edges$cluster, layout=layout.fruchterman.reingold)
