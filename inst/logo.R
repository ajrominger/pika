x <- sad(model = 'fish', par = 0.01)
r <- sad2Rank(x, 100)
foo <- 0.8

pdf('inst/rad.pdf', width = 3, height = 3)
par(mar = c(1, 1, 0, 0))
plot(r, log = 'y', type = 'n', ylim = c(foo, max(r)), axes = FALSE)
polygon(c(1, 1:length(r), length(r)), c(foo, r, foo), 
        col = hsv(0.58, 0.4), lwd = 2)
socorro::logAxis(2, col.axis = 'transparent', lwd = 2)
axis(1, col.axis = 'transparent', lwd = 2)
dev.off()



x <- exp(seq(log(0.01), log(100), length.out = 50))
y <- log(x + 1) * 10000 - 50

pdf('inst/sar.pdf', width = 3, height = 3)
par(mar = c(1, 1, 0, 0))
plot(x, y, log = 'xy', type = 'n', axes = FALSE)
polygon(c(x, max(x)), c(y, min(y)), 
        col = hsv(0.7, 0.4), lwd = 2)
socorro::logAxis(1:2, col.axis = 'transparent', lwd = 2)
dev.off()
