############################################################
# Family utilities
############################################################


initialize_family <- function(family){


  fam <- family$family


  if(!fam %in%
     c(
       "gaussian",
       "binomial",
       "poisson"
     )
  ){

    stop(
      "Family not supported"
    )
  }


  family
}




inverse_link <- function(eta, family){


  switch(
    family$family,


    gaussian =
      eta,


    binomial =
      1/(1+exp(-eta)),


    poisson =
      exp(eta)
  )

}




working_weights <- function(mu, family){


  switch(
    family$family,


    gaussian =
      rep(1,length(mu)),


    binomial =
      mu*(1-mu),


    poisson =
      mu
  )
}
