coxph_R_to_C <- function(X) {
#
#  This function accepts as an argument a list 
#  containing: 
#      the vector of outcomes Y, 
#      the vector of censoring indicators C, 
#      the vector of treatment assignments treatment,
#      and the vector of split indicators split.
#
#  This function returns the squared z-statistic
#  of the split by treatment interaction term
#  in a Cox proportional hazards model.
#

    Y=X[,1]
    C=X[,2]
    treatment=X[,3]
    split=X[,4]

#    if a treatment group only has censored observations, return 0
    if ( length( table(C[which(treatment==0)]) ) == 1 ||
         length( table(C[which(treatment==1)]) ) == 1 ) {

        return(0)
    }

#    save "coxph" warnings in x
    options(warn=1) 

    x=capture.output({
        fit0=coxph(Surv(Y,C)~treatment+split+treatment*split)
    },type="message")

#    if there is a warning, return 0
    if ( length(x) > 0 ) {  # if there is at least 1 warning
        return(0)
    }

#    get z-statistic of split by treatment interaction term
    fit=summary(fit0)
    val=fit$"coefficients"["treatment:split","z"]

    return(val^2)
}