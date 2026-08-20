plot.glinternet.cv = function(x, ...){
    #plot cv curve
    bestIndex = which.min(x$cvErr)
        barwidth = 0.25
        xi = 1:length(x$lambda)
        y = x$cvErr
        delta = x$cvErrStd
        limits = c(min(y-delta), max(y+delta))
        if (diff(limits) == 0) limits = limits + c(-1,1) * max(1,abs(limits[1])) * .04
        plot(xi, y, type="n", xlab="Lambda index", ylab="CV error", xaxt="n", ylim=limits)
        segments(xi-barwidth, y+delta, xi+barwidth, y+delta,col="grey")
        segments(xi-barwidth, y-delta, xi+barwidth, y-delta,col="grey")
        segments(xi, y+delta,xi,y-delta,col="grey")
        lines(xi,y,lwd=2)

        abline(v=bestIndex, lty=3)
        ticks = if (length(xi) == 1) 1 else unique(c(1,seq(2,length(xi),2),length(xi)))
        axis(1, at=ticks, las=1)
}
