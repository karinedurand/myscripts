
p = seq(0,1, length=100)
plot(p, dbeta(p, 1, 5), ylab="density", type ="l", col=4)
lines(p, dbeta(p, 1, 200), type ="l", col=3)
lines(p, dbeta(p, 5, 5), col=2) 
lines(p, dbeta(p, 1, 1), col=1) 
legend(0.8,3, c("Be(1,5)","Be(1,200)","Be(5,5)", "Be(1,1)"),lty=c(1,1,1,1),col=c(4,3,2,1))
