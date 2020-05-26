# 1. reshuffle gene names and h2
chi2.causal.2

h.cutoff <- seq(0.5,0.95,0.05)
h2.shuffle <- chi2.causal.2[,.(h, causal)]
h2.shuffled <- data.table(h2.shuffled)
#LOOP
res <- list()
O.U_causal_not_hot <- data.table(h.cutoff)
O.U_causal_hot <- data.table(h.cutoff)
O.U_not_causal_not_hot <- data.table(h.cutoff)
O.U_not_causal_hot <- data.table(h.cutoff)
for (i in 1:1000){
  
  h2.shuffled$h <- sample(h2.shuffle$h)
  
  for (hcut in h.cutoff){
    higher <- h2.shuffled[h>=hcut, .N, by="causal"][order(causal)]
    lower <- h2.shuffled[h<hcut, .N, by="causal"][order(causal)]
    
    
    for (ca in unique(h2.shuffled$causal)){
      if (!ca %in% higher$causal){
        toadd <- data.table(causal=ca, N=0)
        higher <- rbind(higher, toadd)
      }
      if (!ca %in% higher$causal){
        toadd <- data.table(causal=ca, N=0)
        lower <- rbind(lower, toadd)
      }
    }
    
    table.chi2 <- data.table(causal=unique(h2.shuffled$causal),higher=higher$N, lower=lower$N)
    table.chi2 <- transpose(table.chi2, keep.names = "cond", make.names = "causal")
    
    # O.U_causal_not_hot[,paste0("V", i)] <- table.chi2[cond=="higher"]$`causal not in hotspot`/table.chi2[cond=="lower"]$`causal not in hotspot`
    
    col.name <- paste0("iter_", i)

    O.U_causal_not_hot[h.cutoff==hcut, (col.name) := table.chi2[cond=="higher"]$`causal not in hotspot`/table.chi2[cond=="lower"]$`causal not in hotspot`]
    O.U_causal_hot[h.cutoff==hcut, (col.name) :=table.chi2[cond=="higher"]$`causal in hotspot`/table.chi2[cond=="lower"]$`causal in hotspot`]
    O.U_not_causal_not_hot[h.cutoff==hcut, (col.name) :=table.chi2[cond=="higher"]$`not causal not in hotspot`/table.chi2[cond=="lower"]$`not causal not in hotspot`]
    O.U_not_causal_hot[h.cutoff==hcut, (col.name) := table.chi2[cond=="higher"]$`not causal in hotspot`/table.chi2[cond=="lower"]$`not causal in hotspot`]
  }
}

O.U_causal_not_hot.melt <- melt(O.U_causal_not_hot, id.vars=1)
O.U_causal_not_hot.melt$group <- "causal_not_hot"

O.U_causal_hot.melt <- melt(O.U_causal_hot, id.vars=1)
O.U_causal_hot.melt$group <- "causal_hot"

O.U_not_causal_not_hot.melt <- melt(O.U_not_causal_not_hot, id.vars=1)
O.U_not_causal_not_hot.melt$group <- "not_causal_not_hot"

O.U_not_causal_hot.melt <- melt(O.U_not_causal_hot, id.vars=1)
O.U_not_causal_hot.melt$group <- "not_causal_hot"

rb1 <- rbind(O.U_causal_not_hot.melt, O.U_causal_hot.melt)
rb2 <- rbind(rb1, O.U_not_causal_not_hot.melt)
final_res <- rbind(rb2, O.U_not_causal_hot.melt)

library(stringr)
realdata.melt <- melt(chi2.res[,.(h.cutoff, `O/U.causal.nothot`, `O/U.causal.hot`, `O/U.notcausal.nothot`, `O/U.notcausal.hot`)], id.vars=1)
realdata.melt$variable <- str_replace(realdata.melt$variable, "O/U.", "")
cat <- unique(realdata.melt$variable)

ggsave(filename = "mixed_heritability.pdf", plot=heri, device="pdf", limitsize=F, )

# heri <- final_res[order(group)] %>% 
heri <- ggplot(data=final_res[order(group)], aes(x=factor(h.cutoff),y=as.numeric(value)))+
  geom_boxplot(aes(color=group), outlier.size = 0.5) +
  geom_line(data = realdata.melt[order(variable)], aes(x = as.factor(h.cutoff), y=value, group = variable, color=variable))+
  geom_point(data = realdata.melt[order(variable)], aes(x = as.factor(h.cutoff), y=value, group = variable, color=variable))
plot(heri)

ggplot(data=realdata.melt[order(variable)], aes(x=as.factor(h.cutoff), y=value, group=variable)) +
  geom_line()
