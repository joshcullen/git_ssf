get.cumtime.8dir=function(coord1,mean.time,window1){
  #calculate 8 directions

  direction=c('up','do','le','ri','ur','ul','lr','ll')
  ndirection=length(direction)
  fim=matrix(NA,ndirection*window1,3)
  oo=1
  for (i in 1:ndirection){
    if (direction[i]=='up'){
      final.coord=coord1['y']+window1
      seq1.x=rep(coord1['x'],window1)
      seq1.y=(coord1['y']+1):final.coord
    }
    if (direction[i]=='do'){
      final.coord=coord1['y']-20
      seq1.x=rep(coord1['x'],window1)
      seq1.y=(coord1['y']-1):final.coord
    }    
    if (direction[i]=='ri'){
      final.coord=coord1['x']+20
      seq1.x=(coord1['x']+1):final.coord
      seq1.y=rep(coord1['y'],window1)
    }
    if (direction[i]=='le'){
      final.coord=coord1['x']-20
      seq1.x=(coord1['x']-1):final.coord
      seq1.y=rep(coord1['y'],window1)      
    }
    if (direction[i]=='ur'){
      final.coord.y=coord1['y']+20
      final.coord.x=coord1['x']+20
      seq1.y=(coord1['y']+1):final.coord.y
      seq1.x=(coord1['x']+1):final.coord.x
    }
    if (direction[i]=='ul'){
      final.coord.y=coord1['y']+20
      final.coord.x=coord1['x']-20
      seq1.y=(coord1['y']+1):final.coord.y
      seq1.x=(coord1['x']-1):final.coord.x
    }
    if (direction[i]=='lr'){
      final.coord.y=coord1['y']-20
      final.coord.x=coord1['x']+20
      seq1.y=(coord1['y']-1):final.coord.y
      seq1.x=(coord1['x']+1):final.coord.x
    }
    if (direction[i]=='ll'){
      final.coord.y=coord1['y']-20
      final.coord.x=coord1['x']-20
      seq1.y=(coord1['y']-1):final.coord.y
      seq1.x=(coord1['x']-1):final.coord.x
    }
    tmp=cbind(seq1.y,seq1.x); colnames(tmp)=c('y','x')
    tempo=cumsum(mean.time[tmp])
    
    #adjust time traveling in diagonal to account for greater distance
    #3 pixels up and down are not the same as 3 pixels in the diagonal
    if (direction[i]%in%c('ll','lr','ur','ul')) tempo=tempo*sqrt(2)
    fim[oo:(oo+window1-1),]=cbind(tempo,seq1.x,seq1.y)
    oo=oo+window1
  }
  colnames(fim)=c('cum.time','x','y')
  fim
}