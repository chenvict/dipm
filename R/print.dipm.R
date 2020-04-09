
prin <- function(mess,a) {

    cat(sprintf("\n%s\n",mess))
    print(a)
}



prinm <- function(mess) {

    cat(sprintf("\n%s\n",mess))
}



findrows_bin <- function(data,covar,coval) {
#
#  This function accepts as arguments a BINARY 
#  covariate and one of its values and returns a 
#  vector containing "TRUE"/"FALSE" values, where 
#  "TRUE" denotes that the covariate in a 
#  particular row of data equals the specified 
#  value.
#

      ifrow= ( data[,covar] == coval )
      return(ifrow)
}



findrows_ord <- function(data,covar,coval,sign) {
#
#  This function accepts as arguments an ORDINAL 
#  covariate and its values and returns a vector 
#  containing "TRUE"/"FALSE" values, where "TRUE" 
#  denotes that the covariate in a particular row 
#  of data equals the specified values.
#

    dat=data[,covar]

    if( sign == 1 ) {  # less than or equal to
        ifrow=( dat <= coval )
    } 

    if( sign == 2 ) {  # greater than
        ifrow=( dat > coval )
    }

    return(ifrow)
}



findrows_nom <- function(data,covar,coval,ncat) {
#
#  This function accepts as arguments a NOMINAL 
#  covariate and its values and returns a vector 
#  containing "TRUE"/"FALSE" values, where "TRUE" 
#  denotes that the covariate in a particular row 
#  of data equals any of the specified values.
#

#
#   ... extract nominal categories from "coval" input
#
      vals2=get_nomvals(coval,ncat)
      vals2=sub("\\{","",vals2)
      vals2=sub("\\}","",vals2)
      vals2=strsplit(vals2,split=",")
      vals=unlist(vals2)

      lll=length(vals)
      ifrow=rep(FALSE,nrow(data))

      for (i in 1:lll) {

          ifrow2= ( data[,covar] == vals[i] )
          ifrow=ifrow | ifrow2

      }

      return(ifrow)
}



findrows_node <- function(node,tree,data,ncat) {
#
#  This function accepts as arguments the integer
#  index of a node starting from 1, the tree, and
#  the dataset containing all of the candidate
#  split variables and returns a vector containing 
#  "TRUE"/"FALSE" values, where "TRUE" denotes that 
#  a particular row of data is in the node.
#

#
#   ... get indices of relevant nodes in tree
#
    index=node

    if ( node != 1 ) {

        for ( i in 1:node ) {

            j=length(index)
            current_node=index[j]

            if ( current_node == 1 ) break()

            index=c(index,tree$"parent"[current_node])
        }

        index=sort(index)
    }

#
#   ... get the covariates and values
#
      covars=tree$"splitvar"[index]
      covals=tree$"splitval"[index]
      covar_types=tree$"type"[index]
      signs=tree$"sign"[index]

      ncat0=0
      ncat=c(ncat0,ncat[covars[-1]])

#
#   ... find the rows that have the desired values
#
      ifnode=rep(TRUE,nrow(data))

      lll=length(covars)

      if ( lll > 0 ) {

      for ( i in 1:length(covars) ) {

          if ( covar_types[i] == 0 ) {  # all rows are in root node
              next()
          }

          if ( covar_types[i] == 1 ) {  # binary
              ifnode2=findrows_bin(data,covars[i],covals[i])
          }

          if ( covar_types[i] == 2 ) {  # ordinal
              ifnode2=findrows_ord(data,covars[i],covals[i],signs[i])
          }

          if ( covar_types[i] == 3 ) {  # nominal
              ifnode2=findrows_nom(data,covars[i],covals[i],ncat[i])
          }

          ifnode=ifnode & ifnode2
      }
      }

      return(ifnode)
}



get_nomvals <- function(value,ncat) {
#
#  This function returns the categories of a nominal
#  split.
#
    icat=ncat
    cats0=NULL

    for ( i in 1:ncat ) {

        value=as.numeric(value)
        value=value/2

        if ( floor(value) == value ) {
            icat=icat-1
            next
        }

        cats0=c(cats0,icat)

        if ( value == 0.5 ) break

        icat=icat-1
        value=floor(value)
    }

    cats0=sort(cats0)
    cats=paste("{",cats0[1],sep="")

    if ( length(cats0) > 1 ) {
    
        for ( i in 2:length(cats0) ) {
            cats=paste(cats,",",cats0[i],sep="")
        }
    }

    cats=paste(cats,"}",sep="")

    return(cats)
}



print.dipm <- function(tree_txt,X,types,ncat,method,ntree,print) {
#
#  This function accepts as an argument a re-formatted
#  "tree_txt" object from a tree made in our C code and 
#  prints out the tree in text form to the screen in R.
#

#    for survival methods, best treatment determined by
#    mean survival
    if ( method %in% c(11,12,13,20,21,25) ) {

        for ( i in 1:nrow(tree_txt) ) {

#            get vector denoting which records are in the current node
            ifnode=findrows_node(i,tree_txt,X,ncat)

#            get the kaplan-meier curves by treatment group for 
#            subjects in the current node
            survobj=with(X[ifnode,],Surv(Y[ifnode],C[ifnode] == 1))
            km=survfit(survobj ~ treatment[ifnode],data=X[ifnode,])
            val=survival:::survmean(km,rmean="individual")[[1]]

#            get the mean survival values for each treatment group
            km_means=val[,"*rmean"]

#            get the corresponding treatment group values
            km_trts=rownames(val)
            km_trts=sub(".*=","",km_trts)

#            the best treatment group has the highest mean value
            bestval=km_trts[which.max(km_means)]
            tree_txt$"besttrt"[i]=as.integer(bestval)
        }
    }

#    exclude unwanted columns from tree output
    tree_txt=tree_txt[,-which(colnames(tree_txt) == "parent")]
    tree_txt=tree_txt[,-which(colnames(tree_txt) == "sign")]
    tree_txt=tree_txt[,-which(colnames(tree_txt) == "ntrt0")]
    tree_txt=tree_txt[,-which(colnames(tree_txt) == "ntrt1")]
    tree_txt=tree_txt[,-which(colnames(tree_txt) == "r0")]
    tree_txt=tree_txt[,-which(colnames(tree_txt) == "r1")]
    tree_txt=tree_txt[,-which(colnames(tree_txt) == "p0")]
    tree_txt=tree_txt[,-which(colnames(tree_txt) == "p1")]

#    make splitvars the variable that splits the current node
#    into its subsequent child nodes instead of the variable
#    that split the current node
    for ( i in 1:nrow(tree_txt) ) {

        lchild=tree_txt$lchild[i]  # current left child

        if ( lchild != 0 ) {

            tree_txt$splitvar[i]=tree_txt$splitvar[lchild]
            tree_txt$type[i]=tree_txt$type[lchild]
            tree_txt$splitval[i]=tree_txt$splitval[lchild]

        } else if ( lchild == 0 ) {

            tree_txt$splitvar[i]=NA
            tree_txt$type[i]=0
            tree_txt$splitval[i]=NA
        }
    }

#    make covariate types more readable
    if ( nrow(tree_txt) > 0 ) {

#        make covariate types more readable
        for ( i in 1:nrow(tree_txt) ) {
            if ( tree_txt$type[i] == 1 ) tree_txt$type[i]="bin"
            else if ( tree_txt$type[i] == 2 ) tree_txt$type[i]="ord"
            else if ( tree_txt$type[i] == 3 ) tree_txt$type[i]="nom"
        }

#        make "sign" values more readable
###        for ( i in 2:nrow(tree_txt) ) {
###            if ( tree_txt$sign[i] == 0 ) tree_txt$sign[i]="="
###            else if ( tree_txt$sign[i] == 1 ) tree_txt$sign[i]="LE"
###            else if ( tree_txt$sign[i] == 2 ) tree_txt$sign[i]="GT"
###        }
    }

#    make split values of nominal variables more readable
    ifnominal=any(tree_txt$type == "nom")
    if ( ifnominal == TRUE ) {

        index=which(types == 3)
        nom_colnames=colnames(types)[index]
        nomX=data.frame(X[,colnames(types)[index]])
        nomX=data.frame(lapply(nomX,
                               function(x) as.numeric(levels(x)[x])))
        colnames(nomX)=nom_colnames
        ncat=apply(nomX,2,max)
        temp_covar=colnames(X)[tree_txt$splitvar]

        for ( i in 1:nrow(tree_txt) ) {

            if ( is.na(temp_covar[i]) ) next

            for ( j in 1:ncol(nomX) ) {

                if ( temp_covar[i] == colnames(nomX)[j] ) {
                    newval=get_nomvals(tree_txt$splitval[i],ncat[j])
                    tree_txt$splitval[i]=newval
                }
            }
        }
    }

#    for terminal nodes, set type, lchild, and rchild to NA
    for ( i in 1:nrow(tree_txt) ) {

        lchild=tree_txt$lchild[i]  # current left child

        if ( lchild == 0 ) {
            tree_txt$type[i]=NA
            tree_txt$lchild[i]=NA
            tree_txt$rchild[i]=NA
        }
    }

#    print the tree
    if ( print == TRUE ) {

#        print header
        if ( method == -1 ) {

            prinm(paste("PM Tree ",
                        "(Continuous Y, ",
                        "2 treatments",
                        "):",
                        sep=""))

        } else if ( method == 6 ) {

            prinm(paste("DIPM Tree ",
                        "(Continuous Y, ",
                        "2 treatments, ",
                        "no mtry, ",
                        "ntree=",ntree,
                        "):",
                        sep=""))

        } else if ( method == 7 ) {

            prinm(paste("DIPM Tree ",
                        "(Continuous Y, ",
                        "2 treatments, ",
                        "yes mtry, ",
                        "ntree=",ntree,
                        "):",
                        sep=""))

        } else if ( method == 11 ) {

            prinm(paste("PM Tree ",
                        "(Survival Y, ",
                        "2 treatments",
                        "):",
                        sep=""))

        } else if ( method == 12 ) {

            prinm(paste("DIPM Tree ",
                        "(Survival Y, ",
                        "2 treatments, ",
                        "no mtry, ",
                        "ntree=",ntree,
                        "):",
                        sep=""))

        } else if ( method == 13 ) {

            prinm(paste("DIPM Tree ",
                        "(Survival Y, ",
                        "2 treatments, ",
                        "yes mtry, ",
                        "ntree=",ntree,
                        "):",
                        sep=""))

        } else if ( method == 20 ) {

            prinm(paste("DIPM Tree ",
                        "(Survival Y, ",
                        "2+ treatments, ",
                        "no mtry, ",
                        "ntree=",ntree,
                        "):",
                        sep=""))

        } else if ( method == 21 ) {

            prinm(paste("DIPM Tree ",
                        "(Survival Y, ",
                        "2+ treatments, ",
                        "yes mtry, ",
                        "ntree=",ntree,
                        "):",
                        sep=""))

        } else if ( method == 22 ) {

            prinm(paste("DIPM Tree ",
                        "(Continuous Y, ",
                        "2+ treatments, ",
                        "no mtry, ",
                        "ntree=",ntree,
                        "):",
                        sep=""))

        } else if ( method == 23 ) {

            prinm(paste("DIPM Tree ",
                        "(Continuous Y, ",
                        "2+ treatments, ",
                        "yes mtry, ",
                        "ntree=",ntree,
                        "):",
                        sep=""))

        } else if ( method == 24 ) {

            prinm(paste("PM Tree ",
                        "(Continuous Y, ",
                        "2+ treatments",
                        "):",
                        sep=""))

        } else if ( method == 25 ) {

            prinm(paste("PM Tree ",
                        "(Survival Y, ",
                        "2+ treatments",
                        "):",
                        sep=""))
        }

        print(tree_txt,row.names=FALSE)
    }

#    add split variable names to tree output
    tree_txt=data.frame(tree_txt,splitvar_name=NA)
    tree_txt=tree_txt[,c(1,2,10,3:9)]

    for ( i in 1:nrow(tree_txt) ) {

        var_i=tree_txt$splitvar[i]

        if ( is.na(var_i) ) next()

        tree_txt$splitvar_name[i]=colnames(X)[var_i]
    }

    return(tree_txt)
}