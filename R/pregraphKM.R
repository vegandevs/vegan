`pregraphKM` <-
    function(matrice)
{
    `row.col.number` <- function(mat,number){
        nr<-nrow(mat)
        nc<-ncol(mat)
        mod<-number %% nr
        div<-number/nr
        ##First column
        if(mod!=0 & div>1){
            nr.f<-mod
            nc.f<-trunc(div)+1
        }else{
            if(div<=1){
                nc.f<-1
                nr.f<-number
            }else{
                if(mod==0){
                    nc.f<-div
                    nr.f<-nr
                }
            }
        }
        list(nr=nr.f,nc=nc.f)
    }

    ## Beginning of the function

    for(k in 1:(ncol(matrice)-1)){
        i=1
        j=1
        tmp<-0
        while(j <= max(matrice[,k])){
            if(i==1){
                mat<-table(matrice[,k],matrice[,k+1])
                number<-which.max(mat)
                tmp<-row.col.number(mat,number)
                tmp0<-tmp
                ## Change les indices 
                if(tmp$nr!=tmp$nc){
                    find.nc<-which(matrice[,k+1]==tmp$nc)
                    find.nr<-which(matrice[,k+1]==tmp$nr)
                    matrice[find.nc,k+1]<-tmp$nr
                    matrice[find.nr,k+1]<-tmp$nc
                }else{}
                i=2
            }else{
                mat<-table(matrice[,k],matrice[,k+1])
                mat[tmp0$nr,]<-0
                mat[,tmp0$nr]<-0
                number<-which.max(mat)
                tmp<-row.col.number(mat,number)
                ## Change les indices 
                if(tmp$nr!=tmp$nc){
                    find.nc<-which(matrice[,k+1]==tmp$nc)
                    find.nr<-which(matrice[,k+1]==tmp$nr)
                    matrice[find.nc,k+1]<-tmp$nr
                    matrice[find.nr,k+1]<-tmp$nc
                }else{}
                tmp0$nr<-c(tmp0$nr,tmp$nr)
            }
            j=j+1
        }
    }
    matrice
}
