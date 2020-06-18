# find gene pairs in cis

test <- gene.location.order[chr.A==chr.B][,("diff"):= start.A-start.B]
# 1Mb == 1,000,000bp
any(T %in% (abs(test$diff) > 1000000))

test[abs(diff) > 1000000]


cispairs <- test[abs(diff) < 1000000]
