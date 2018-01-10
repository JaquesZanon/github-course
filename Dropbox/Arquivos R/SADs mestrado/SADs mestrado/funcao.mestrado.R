
ttt<-read.table('clipboard',header=T)


library(sads)
?sads
dados<-read.table('teste.txt', header=T)
geral<-read.table('clipboard', header=T,row.names=1)


teta<-sad(geral,'zsm')

### função utilizada no mestrado
sad<-function(dados,type='all'){
  
  ### tirando os zeros
  resu<-list()
  for(i in 1:ncol(dados)){
    col<-dados[,i]
    resu[[i]]<-col[col!=0]
    resu[[i]]<-as.integer(resu[[i]])
  }#tira os zeros resulta em uma lista com cada vetor de abundancia.
  
  ### ajustando os modelos
  if(type=='all'){  
  mat<-matrix(NA,ncol(dados),3)
  mat2<-matrix(NA,ncol(dados),4)
  rownames(mat)<-colnames(dados)
  rownames(mat2)<-colnames(dados)
  colnames(mat)<-c('alpha','theta','prop')
  colnames(mat2)<-c('meanlog','sdlog','mu-poisson','sig-poisson')
  for(i in 1:ncol(dados)){
    
      ls<-fitsad(resu[[i]],'ls')
      mat[i,1]<-coef(ls)# coef log-series
    
      mzsm<-fitsad(resu[[i]],"mzsm")
      mat[i,2]<-coef(mzsm)# coef zero-sum neutral theory
    
      geom<-fitsad(resu[[i]],"geom")
      mat[i,3]<-coef(geom)# coef geometric series
    
    
  }# analisa as fitdsad logseries, mzsm, e a geo
  for(i in 1:ncol(dados)){
      mat2[i,1]<-(coef(fitsad(resu[[i]],"lnorm",trunc=0.5,trueLL=F))[1])
      mat2[i,2]<-(coef(fitsad(resu[[i]],"lnorm",trunc=0.5,trueLL=F))[2])#coef lognormal
      mat2[i,3]<-(coef(fitsad(resu[[i]],"poilog"))[1])
      mat2[i,4]<-(coef(fitsad(resu[[i]],"poilog"))[2])    
  }# analisa lognoraml e poisson lognormal
  
  ### agora os AICs
  AIC<-list()
  
  for(i in 1:ncol(dados)){
    ls<-fitsad(resu[[i]],'ls')
    mzsm<-fitsad(resu[[i]],"mzsm")
    geom<-fitsad(resu[[i]],"geom")
    lognorm<-fitsad(resu[[i]],"lnorm",trunc=0.5,trueLL=F)
    poilog<-fitsad(resu[[i]],"poilog")
    AIC[[i]]<-AICtab(ls,mzsm,geom,lognorm,poilog)
    
   
  }
  return(list(ajuse1=mat,ajuste2=mat2,deltaAIC=AIC))
  
}
  if(type=='lognormal'){
    mat2<-matrix(NA,ncol(dados),4)
    rownames(mat2)<-colnames(dados)
    colnames(mat2)<-c('meanlog','sdlog','mu-poisson','sig-poisson')
    for(i in 1:ncol(dados)){
      mat2[i,1]<-(coef(fitsad(resu[[i]],"lnorm",trunc=0.5,trueLL=F))[1])
      mat2[i,2]<-(coef(fitsad(resu[[i]],"lnorm",trunc=0.5,trueLL=F))[2])#coef lognormal
      mat2[i,3]<-(coef(fitsad(resu[[i]],"poilog"))[1])
      mat2[i,4]<-(coef(fitsad(resu[[i]],"poilog"))[2])
    }
    return(list(ajuste=mat2))
      
  }
    
  
  if(type=='logseries'){
    mat<-matrix(NA,ncol(dados),1)
    rownames(mat)<-colnames(dados)
    colnames(mat)<-c('alpha')
    
    for(i in 1:ncol(dados)){
      ls<-fitsad(resu[[i]],'ls')
      mat[i,1]<-coef(ls)# coef log-series
    }
    return(list(ajuste=mat))
  }
  if(type=='zsm'){
    mat<-matrix(NA,ncol(dados),1)
    rownames(mat)<-colnames(dados)
    colnames(mat)<-c('theta')
    for(i in 1:ncol(dados)){
      zsm<-fitsad(resu[[i]],'mzsm')
      mat[i,1]<-coef(zsm)# coef log-series
    }
    
    
    return(list(ajuste=mat))
  }  
  if(type=='geometric'){
    mat<-matrix(NA,ncol(dados),1)
    rownames(mat)<-colnames(dados)
    colnames(mat)<-c('geom')
    for(i in 1:ncol(dados)){
      geom<-fitsad(resu[[i]],'geom')
      mat[i,1]<-coef(geom)# coef log-series
    }
    
    return(list(ajuste=mat))
  }  
    
}


sad(ttt)

