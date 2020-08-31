#barplot LCBD
windows()
layout(matrix(c(1,2,3,4),nrow= 4, ncol= 1, byrow=TRUE))
par(mar= c(0,10,8,7))
par(las= 1)
barplot(LCBD_barplot[,4],
        xaxs = "i",
        yaxs="i",
        ylim = c(0, 585), xlim = c(-0.5, 14),
        space= 0,
        col= "black",
        axisnames= FALSE,
        axes= FALSE,
        ylab= "Altitude")
axis(side= 2, at= c(0, seq(100, 500, by= 100), 585))

dim_names<- c("Taxonomic LCBD", "Phylogenetic LCBD", "Functional LCBD")
par(mar= c(4,10,1,12.5))
for(i in 1:3){
    barplot(LCBD_barplot[,i],
            names.arg= rownames(LCBD_barplot),
            xaxs = "i",
            yaxs="i",
            ylim = c(0, 0.40), xlim = c(-0.5, 14),
            col= "white",
            axisnames= FALSE,
            axes= FALSE,
            ylab= dim_names[i]
    )
    axis(side = 1, 
         at = c(-0.5,0.7,1.9,3.2,4.3,5.5,6.7,7.9,9.1,10.3,11.5,12.7,13.5), 
         labels= c("","LIvi", "LVaca", "LBri", "UIvi", "LDou", "UVaca",
                   "MBri", "MDou", "SMaria", "UBri", "UDou", ""))
    par(las= 1)
    axis(side = 2, at= , ylab= dim_names[i])   
}
amb
mean_altitude_sort<-sort(sapply(levels(amb$porcao), function(x) {mean(amb$ALTITUDE[amb$porcao == x])}))
mean_altitude<-sapply(levels(amb$porcao), function(x) {mean(amb$ALTITUDE[amb$porcao == x])})
match(mean_altitude_sort,mean_altitude)
LCBD_barplot <-cbind(LCBD_barplot[match(mean_altitude_sort,mean_altitude),c(1,2,3)], mean_altitude_sort)

#para cada barplot separado
windows()
barplot(LCBD_barplot[,3],
        names.arg= rownames(LCBD_barplot),
        xaxs = "i",
        yaxs="i",
        ylim = c(0, 0.20), xlim = c(-0.5, 14),
        col= "white",
        axisnames= FALSE,
        axes= FALSE,
        ylab= dim_names[3]
)
axis(side = 1, 
     at = c(-0.5,0.7,1.9,3.2,4.3,5.5,6.7,7.9,9.1,10.3,11.5,12.7,13.5), 
     labels= c("","LIvi", "LVaca", "LBri", "UIvi", "LDou", "UVaca",
               "MBri", "MDou", "SMaria", "UBri", "UDou", ""))
par(las= 1)
axis(side = 2, at= , ylab= dim_names[3]) 
max(LCBD_barplot[,3])
