library(grid)

plotsirD = function(input,outfile = NULL)
{
  if(!is.null(outfile))png(outfile,height=8,width=15,units='cm',res=300,pointsize=10)
  #pushViewport(plotViewport(c(1,1,1,1)))
  pushViewport(viewport(layout=grid.layout(nrow=1,ncol=2)))
  xlm=c(0,max(input$t))
  ylm1=c(0,100)
  ylm2=c(0,max(input$I)*1.05)/1000
  pushViewport(viewport(layout.pos.col=1,layout.pos.row=1))
  pushViewport(plotViewport(c(4,4,1,1),xscale=xlm,yscale=ylm1,clip=TRUE))
  grid.lines(input$t,100*input$R/(input$S[1]+input$I[1]+input$R[1]),gp=gpar(col='darkorange',lwd=2),default.units = 'native')
  popViewport()
  pushViewport(plotViewport(c(4,4,1,1),xscale=xlm,yscale=ylm1))
  grid.lines(c(0,1,1,0,0),c(0,0,1,1,0))
  tick_locations = 0:floor(xlm[2]/30)
  grid.xaxis(at=30*tick_locations,label=tick_locations)
  grid.yaxis(at=seq(0,100,20),label=paste0(seq(0,100,20),'%'))
  grid.text('Time (mo)',y=unit(-3,'lines'))
  grid.text('Recovered',x=unit(-3.25,'lines'),rot=90)
  popViewport()
  popViewport()
  pushViewport(viewport(layout.pos.col=2,layout.pos.row=1))
  pushViewport(plotViewport(c(4,4,1,1),xscale=xlm,yscale=ylm2,clip=TRUE))
  grid.lines(input$t,input$I/1000,gp=gpar(col='darkorange',lwd=2),default.units = 'native')
  popViewport()
  pushViewport(plotViewport(c(4,4,1,1),xscale=xlm,yscale=ylm2,clip=FALSE))
  grid.lines(c(0,1,1,0,0),c(0,0,1,1,0))
  
  tick_locations = 0:floor(xlm[2]/30)
  grid.xaxis(at=30*tick_locations,label=tick_locations)
  yticks = pretty(ylm2)
  yticks = yticks[yticks<ylm2[2]]
  grid.yaxis(at=yticks,label=paste0(yticks,'k'))
  grid.text('Time (mo)',y=unit(-3,'lines'))
  grid.text('Infected',x=unit(-3.25,'lines'),rot=90)
  popViewport()
  popViewport()
  #popViewport()
  if(!is.null(outfile))dev.off()
  
}


plotsirS = function(input,outfile = NULL)
{
  if(!is.null(outfile))png(outfile,height=8,width=15,units='cm',res=300,pointsize=10)
  #pushViewport(plotViewport(c(1,1,1,1)))
  pushViewport(viewport(layout=grid.layout(nrow=1,ncol=2)))
  xlm=c(0,max(input$t))
  ylm1=c(0,100)
  ylm2=c(0,max(input$I)*1.05)/1000
  pushViewport(viewport(layout.pos.col=1,layout.pos.row=1))
  pushViewport(plotViewport(c(4,4,1,1),xscale=xlm,yscale=ylm1,clip=TRUE))
  grid.lines(input$t,100*input$R/(input$S[1]+input$I[1]+input$R[1]),gp=gpar(col='#003171',lwd=2),default.units = 'native')
  popViewport()
  pushViewport(plotViewport(c(4,4,1,1),xscale=xlm,yscale=ylm1))
  grid.lines(c(0,1,1,0,0),c(0,0,1,1,0))
  tick_locations = 0:floor(xlm[2]/30)
  grid.xaxis(at=30*tick_locations,label=tick_locations)
  grid.yaxis(at=seq(0,100,20),label=paste0(seq(0,100,20),'%'))
  grid.text('Time (mo)',y=unit(-3,'lines'))
  grid.text('Recovered (%)',x=unit(-3.5,'lines'),rot=90)
  popViewport()
  popViewport()
  pushViewport(viewport(layout.pos.col=2,layout.pos.row=1))
  pushViewport(plotViewport(c(4,4,1,1),xscale=xlm,yscale=ylm2,clip=TRUE))
  grid.lines(input$t,input$I/1000,gp=gpar(col='#003171',lwd=2),default.units = 'native')
  popViewport()
  pushViewport(plotViewport(c(4,4,1,1),xscale=xlm,yscale=ylm2,clip=FALSE))
  grid.lines(c(0,1,1,0,0),c(0,0,1,1,0))
  
  tick_locations = 0:floor(xlm[2]/30)
  grid.xaxis(at=30*tick_locations,label=tick_locations)
  yticks = pretty(ylm2)
  yticks = yticks[yticks<ylm2[2]]
  grid.yaxis(at=yticks,label=paste0(yticks,'k'))
  
  grid.text('Time (mo)',y=unit(-3,'lines'))
  grid.text('Infected',x=unit(-3,'lines'),rot=90)
  popViewport()
  popViewport()
  #popViewport()
  if(!is.null(outfile))dev.off()
  
}


