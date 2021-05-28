
###
## Function which contains details of parameter names and units
## NOTE: that this function is not in use
###

# This function is based on an original Matlab function development by A. A. Bloom (UoE, now at the Jet Propulsion Laboratory).
# Translation to R and subsequent modifications by T. L Smallman (t.l.smallman@ed.ac.uk, UoE).

parameter_details<- function (modelname,parameter_type,ctessel_pft) {
    if (modelname == "DALEC_CDEA") {
        output=list(npars=23)
        output$parameter_symbols=c("p1","p2","p3","p4","p5","p6","p7","p8","p9","p10","p11","p12","p13","p14","p15","p16","p17","i1","i2","i3","i4","i5","i6")
        output$parameter_names=c("D_{rate}","F_{Rau}","F_{Fol}","F_{Roots}","L","TOR_{Wood}","TOR_{Roots}","MR_{litter}"
                                ,"MR_{SOM}","Temp_{par}","C_{eff}","B_{day}","F_{Lab}","R_{L}","F_{day}","R_{F}","LMA"
                                ,"C_{LABILE}","C_{FOLIAR}","C_{ROOT}","C_{WOOD}","C_{LITTER}","C_{SOM}")
        output$parameter_units=c("day^{-1}","fraction","fraction","fraction","years","day^{-1}","day^{-1}","day^{-1}","day^{-1}"
                                ,"(unitless)","(unitless)","doy","fraction","days","doy","days","gC m^{-2}","gC m^{-2}","gC m^{-2}"
                                ,"gC m^{-2}","gC m^{-2}","gC m^{-2}","gC m^{-2}")
    } else {
        print("model name not found in 'parameter_details'")
    }
    return(output)
} # end function parameter_details
