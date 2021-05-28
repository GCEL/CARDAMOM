
###
## Function to extract met drivers from binary file
###

# This function is based on an original Matlab function development by A. A. Bloom (UoE, now at the Jet Propulsion Laboratory).
# Translation to R and subsequent modifications by T. L Smallman (t.l.smallman@ed.ac.uk, UoE).
# Exceptions are given within specific functions.

read_binary_file_format<- function(infile) {

      #print("Beginning read of binary input files...")
      # Open and read the DALEC binary driver file
      # open this chains binary file into R, instructing 'r' to read and 'b' for binary
      bob=file(infile,'rb') ; nos_var=1e6
      bd=readBin(bob, double(),nos_var)
      # keep reading until we have read all that can be read
      set1=NA
      while (length(set1) > 0) {
	      set1=readBin(bob, double(),nos_var)
	      bd=append(bd,set1)
      }
      # now close this chain
      close(bob)

      # begin preparing to disagregate the different sections of the file
      k=0
      # extract static data (100 places)
      static=bd[(k+1):(k+100)]
      k=k+100
      # extract priors
      pr=bd[(k+1):(k+100)]
      k=k+100
      # extract prior uncertainties (50 places)
      pru=bd[(k+1):(k+100)]
      k=k+100
      # other priors (50 places)
      opr=bd[(k+1):(k+100)]
      k=k+100
      # other prior uncertainties (50 places)
      opru=bd[(k+1):(k+100)]
      k=k+100

      # store prior information
      md=list(parpriors=pr,parpriorunc=pru,otherpriors=opr,otherpriorunc=opru)

      # id code (not currently used)
      md$id=static[1]
      # latitude
      md$lat=static[2]
      # number of days in simulation
      md$nodays=static[3]
      # number of met fields
      md$nomet=static[4]
      # number of observation streams
      md$noobs=static[5]
      # implement reality checks (1=yes,0=no)
      md$RC=static[6]
      # ctessel pft
      md$ctessel_pft=static[7]
      # yield if forest
      md$yield=static[8]
      # age if forest
      md$age=static[9]
      # 10 will be nos_pars
      ###
      # start searching EDCs from anywhere (1) or from prescribed starting point (0)
      md$rc_random_search= static[11] == 1
      # Top sand %
      md$top_sand=static[12]
      # bot sand %
      md$bot_sand=static[13]
      # Top clay %
      md$top_clay=static[14]
      # Bot clay %
      md$bot_clay=static[15]


      # extract temporal data (met and obs)
      tempdata=bd[(k+1):(k+((md$nomet+md$noobs)*md$nodays))]
      # restructure correctly
      tempdata=array(tempdata,dim=c(md$nomet+md$noobs,md$nodays))

      # pass all met data
      md$met=t(tempdata[1:md$nomet,1:md$nodays])
      # pass all observations (if any)
      # current defaults are:
      # 1st column = GPP
      # 2nd column = LAI
      # 3rd column = NEE
      md$obs=t(tempdata[(md$nomet+1):(md$nomet+md$noobs),1:md$nodays])

      # pass back out information
      return(md)

} # end of function
## Use byte compile
read_binary_file_format<-cmpfun(read_binary_file_format)
