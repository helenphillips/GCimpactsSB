setwd("C:/Users/hp39wasi/WORK/GCimpactsSB")
## I got these numbers by hand
## These are the number of papers screened

gcd <- data.frame(gcd = c("LUI", "Fragmentation", "Climate Change", "Pollution", "Nutrient \nEnrichment", "Invasive \nspecies"),
                  total = c(502, 28, 113, 194, 218, 44))

# jpeg(filename = "NumberofArticles.jpeg", quality = 100, res = 300, width = 2000, height = 1000)

par(mar=c(6, 4, 1, 1))
par(bg=NA)


b <- barplot(gcd$total, space=0.1, ylab = "Number of articles")
end_point = 0 + nrow(gcd) + nrow(gcd)-1 #this is the line which does the trick (together with barplot "space = 1" parameter)

text(b[,1], par("usr")[3]-15, 
     srt = 60, adj= c(1, 0.2), xpd = TRUE,
     labels = gcd$gcd, cex=1)

# dev.off()



b <- barplot(t$Freq, space=0.1, ylab = "Number of papers")
#axis(1, at = b[,1], labels = t$t, space=1)
text(b[,1], par("usr")[3]-15, 
     srt = 60, adj= c(1, 0.2), xpd = TRUE,
     labels = gcd$gcd, cex=1)
dev.off()


## NUMBER OF EXTRACTED DATA

# got all these numbers by hand

dat <- data.frame(gcd = c("LUI", "Fragmentation", "Climate \nChange", "Pollution", "Nutrient \nEnrichment", "Invasive \nspecies"),
                  articles = c(92, 22, 95, 57, 122, 43), 
                  cases = c(295, 96, 493, 343, 417, 154),
                  todo = c(313/2,4,16,125/2,68/2,0))

par(mar = c(6, 4, 1, 7))

b <- barplot(dat$articles, xaxs = "i", ylab = "Number of articles (white bars)", xlim = c(0, 7), col="white")
par(new=TRUE)
plot(b,dat$cases,xaxs = "i", xlim=c(0,7),pch = 19,col="red",axes=FALSE,ylim=c(0,500),ann=FALSE)
axis(4,at=seq(0,500,50), las = 2, pos = 7.3)
mtext("Number of cases (red dots)", line =  4.5, side = 4)


text(b[,1], par("usr")[3]-15, 
     srt = 60, adj= c(1, 0.5), xpd = TRUE,
     labels = dat$gcd, cex=1)



stackedDat <- matrix(c(dat$articles, dat$todo), nrow = 2, byrow = TRUE)
colnames(stackedDat) <- dat$gcd


jpeg(filename = "NumberofArticles.jpeg", quality = 100, res = 300, width = 3500, height = 2000)

par(mar = c(6, 6, 1, 7))


b <- barplot(stackedDat, 
        col=colors()[c(23,89)] , 
        border="white", 
        space=0.04, 
        xlim = c(0,7),
        font.axis=2, 
        xlab="Global Change Driver",
        ylab="Number of articles \n(brown = data extracted, \ngreen = estimated still to extract)")
par(new=TRUE)
tt <- plot(b,dat$cases,pch = 19,col="black" ,axes=FALSE,ylim=c(0,500),ann=FALSE,  xlim=c(0,7))
axis(4,at=seq(0,500,50), las = 2, pos = 6.3)
mtext("Number of cases (black dots)", line =  -2, side = 4)
dev.off()
