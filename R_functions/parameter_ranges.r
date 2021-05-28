
###
## Function stores parameter information for DALECv2 / DALEC_CDEA
## NOTE: this function is not currently in use
###

# This function is based on an original Matlab function development by A. A. Bloom (UoE, now at the Jet Propulsion Laboratory).
# Translation to R and subsequent modifications by T. L Smallman (t.l.smallman@ed.ac.uk, UoE).

parameter_ranges<- function(modelname,parameter_type,ctessel_pft){

    if (modelname == "DALEC_CDEA") {

        # decomposition rate (day-1)
        parmin=0.000001
        parmax=0.01
        parlog=1

        # fraction of GPP respired
        parmin=append(parmin,0.3)
        parmax=append(parmax,0.7)
        parlog=append(parlog,0)

        #Fraction of (1-fgpp) to foliage*/
        parmin=append(parmin,0.01)
        parmax=append(parmax,0.5)
        parlog=append(parlog,1)

        # Fraction of (1-fgpp) to roots*/
        parmin=append(parmin,0.01)
        parmax=append(parmax,1)
        parlog=append(parlog,1)

        #Leaf Lifespan*/
        #Wright et al. 2004*/
        parmin=append(parmin,1.001)
        parmax=append(parmax,8)
        parlog=append(parlog,1)

        #TOR wood*/
        parmin=append(parmin,0.00001)
        parmax=append(parmax,0.005)
        parlog=append(parlog,1)

        #TOR roots*/
        parmin=append(parmin,0.0001)
        parmax=append(parmax,0.01)
        parlog=append(parlog,1)

        #TOR litter*/
        parmin=append(parmin,0.0003)
        parmax=append(parmax,0.01)
        parlog=append(parlog,1)

        #TOR SOM*/
        parmin=append(parmin,0.00001)
        parmax=append(parmax,0.003)
        parlog=append(parlog,1)

        #Temp factor* = Q10 = 1.2-2.2*/
        parmin=append(parmin,0.018)
        parmax=append(parmax,0.08)
        parlog=append(parlog,0)

        #Canopy Efficiency*/
        parmin=append(parmin,1)
        parmax=append(parmax,80)
        parlog=append(parlog,0)

        #Bday*/
        parmin=append(parmin,1)
        parmax=append(parmax,365.25)
        parlog=append(parlog,0)

        #Fraction to Clab*/
        parmin=append(parmin,0.01)
        parmax=append(parmax,0.5)
        parlog=append(parlog,1)

        #PERIODS must be > 3.45*dt

        #Clab Release period*/
        parmin=append(parmin,10)
        parmax=append(parmax,100)
        parlog=append(parlog,1)

        #Fday*/
        parmin=append(parmin,1)
        parmax=append(parmax,365.25)
        parlog=append(parlog,0)

        #Leaf fall period*/
        parmin=append(parmin,20)
        parmax=append(parmax,200)
        parlog=append(parlog,1)

        #LMA*/
        #Kattge et al. 2011*/
        parmin=append(parmin,10)
        parmax=append(parmax,400)
        parlog=append(parlog,1)

        #INITIAL VALUES DECLARED HERE*/

        #C labile*/
        inimin=10
        inimax=1000.0
        inilog=1

        #C foliar*/
        inimin=append(inimin,10)
        inimax=append(inimax,1000.0)
        inilog=append(inilog,1)

        #C roots*/
        inimin=append(inimin,10)
        inimax=append(inimax,1000.0)
        inilog=append(inilog,1)

        #C_agb*/
        inimin=append(inimin,10.0)
        inimax=append(inimax,50000.0)
        inilog=append(inilog,1)

        #C litter*/
        inimin=append(inimin,10)
        inimax=append(inimax,1000.0)
        inilog=append(inilog,1)

        #C_som*/
        inimin=append(inimin,10.0)
        inimax=append(inimax,1000000.0)
        inilog=append(inilog,1)

    }
    # create output list
    output=list(allparsmin=c(parmin,inimin),allparsmax=c(parmax,inimax))

    # normalizing each parameter distribution
    output$logparstdev=(log(output$allparsmax./output$allparsmin))/4

    return(output)

} # end function parameter_ranges
