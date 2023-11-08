"""
17 April 2018 - JFE
The object has been renamed to CARDAMOM_C and only works with the 
no longer supported C version of CARDAMOM

15 February 2018 - JFE
The argument path2projects has been added to point the CARDAMOM object
to where the project files should be saved to allow more flexibility

25 September 2015 - JFE
The wrapper now writes the MODEL ID as first value in the data files to 
allow the framework to choose which model to use

25 September 2014 - JFE
Data downloaded from the cluster will be ordered in sub-folders with run number

25 July 2014 - JFE
Added one method to copy a project to a newly named one <copy_project>.
Added one method to reperform all the steps of the setup <resetup>

22 May 2014 - JFE
Added some methods to write the input data locally, upload on the cluster

15 May 2014 - JFE
This file defines the class of object CARDAMOM. It stores the parameters
of a CARDAMOM project setup and provides methods to compile code, write input,
execute the MCMC, read output, etc...
"""

#load modules
import os, cPickle, time, struct         
import datetime as dt       
import numpy as np
from netCDF4 import Dataset
import shutil, socket

class CARDAMOM_C(object):
    
    def __init__(self, **kwargs):
        """
        Defines the CARDAMOM object 
        """
       
        #set the paths for the current project
        self.set_paths()
 
        if 'path2projects' in kwargs:
            self.paths['projects'] = kwargs['path2projects']

        # check if project details have been defined
        if 'project_name' in kwargs:
            self.project_name = kwargs['project_name']
        else:
            self.project_name = raw_input('Enter project name: ')
      
        # load if project exists
        
        if self.project_name in os.listdir(self.paths["projects"]): 
            if self.project_name+".pData" in os.listdir('/'.join([self.paths["projects"],self.project_name])):
                print "Project \"%s\" found" % self.project_name          
                self.load_project()
            else:
                print "Project directory found but file %s.pData could not be found" % self.project_name
                self.new_project() 
        else:
            os.mkdir('/'.join([self.paths["projects"],self.project_name]))
            self.new_project(**kwargs)

    def new_project(self,**kwargs):
        """
        This method creates the variables needed by CARDAMOM but not yet set to definite values:
        - nPoints.............. the number of points to run CARDAMOM on
        - latitude............. a sequence of length nPoints with the latitude 
        - longitude............ a sequence of length nPoints with the longitude
        - drivers.............. an array of dimensions [nPoints, nTsteps, nDrivers]
        - observations......... an array of dimensions [nPoints, nTsteps, nObs]

        nPoints, nTsteps and nObs will be defined automatically from the shape of 

        """
        
        # This is initialized here but will contain a dictionary will all the project details
        self.details = {}
        if 'project_type' in kwargs:
            self.project_type = kwargs['project_type']
        else:
            projtype=raw_input("Enter project number of type: \n1\tDALEC_GSI\n2\tDALEC_CDEA_LU_FIRES\n3\tDALEC_GSI_newalloc\n4\tDALEC_CDEA_LU_FIRES_LIU\n5\tDALEC_GSI_LIU\n6\tDALEC_CDEA_LU_FIRES_ET\n7\tDALEC_CDEA_CWD_LU_FIRES\n8\tDALEC_GSI_GLEAM\n9\tDALEC_CDEA_GLEAM\n10\tDALEC_CDEA_HBV\n")
            if projtype != "":       
                self.project_type=projtype
                if self.project_type == "1":
                    self.project_type = "DALEC_GSI"
                elif self.project_type == "2":
                    self.project_type = "DALEC_CDEA_LU_FIRES"
                elif self.project_type == "3":
                    self.project_type = "DALEC_GSI_newalloc"
                elif self.project_type == "4":
                    self.project_type = "DALEC_CDEA_LU_FIRES_LIU"
                elif self.project_type == "5":
                    self.project_type = "DALEC_GSI_LIU"
                elif self.project_type == "6":
                    self.project_type = "DALEC_CDEA_LU_FIRES_ET"
                elif self.project_type == "7":
                    self.project_type = "DALEC_CDEA_CWD_LU_FIRES"
                elif self.project_type == "8":
                    self.project_type = "DALEC_GSI_GLEAM"
                elif self.project_type == "9":
                    self.project_type = "DALEC_CDEA_GLEAM"
                elif self.project_type == "10":
                    self.project_type = "DALEC_CDEA_LU_FIRES_HBV"


            else:
                self.project_type="CARDAMOM_LU"

        self.save_project()


    def save_project(self):
        """
        Saves the object attributes in a file named <project_name>.pData
        """   
    
        f = open('/'.join([self.paths["projects"],self.project_name,self.project_name+".pData"]),"wb")
        cPickle.dump(self.__dict__,f,2)
        f.close()

    def copy_project(self,newproject):
        """
        Copies the project into a <newproject>
        """
        
        keep_name = self.project_name
        self.project_name = newproject
        if newproject not in os.listdir(self.paths["projects"]):
            os.mkdir(self.paths["projects"]+"/"+newproject)
        f = open('/'.join([self.paths["projects"],newproject,newproject+".pData"]),"wb")
        cPickle.dump(self.__dict__,f,2)
        f.close()
        self.project_name = keep_name


    def load_project(self):
        """
        Loads the object attributes stored in a cPickle file named <project_name>.pData
        """

        f = open('/'.join([self.paths["projects"],self.project_name,self.project_name+".pData"]),"rb")
        tmp_dict = cPickle.load(f)
        f.close()

   #     if "paths" in tmp_dict:
   #         tmp_dict.pop("paths",None)

        self.__dict__.update(tmp_dict)

    def set_paths(self):
        """
        Loads paths to the different directories required by CARDAMOM
        """

        hostname = socket.gethostname()
        if "paths" in dir(self):
            print "Paths found in project"
        else:
            print "Paths need to be defined"
            if "default_paths_%s.pData" % (hostname) in os.listdir(os.getcwd()):
                usedefault = raw_input("Use default paths <y/n>? ")
                if usedefault == 'y':
                    f = open("default_paths_%s.pData" % (hostname),"r")
                    self.paths = cPickle.load(f)
                    f.close()
                else:
                    self.define_paths()
            else:
                print "No default paths found... defining paths now"
                self.define_paths()

    def define_paths(self):
        """
        Define the paths used by the project and asked whether to use them
        as default ones
        """

        hostname = socket.gethostname()
        self.paths = {}
        
        keepcurrent = raw_input("Keep current directory (%s) as root one <y/n>? " % os.getcwd())
        
        if keepcurrent == 'y':
            self.paths["CARDAMOM"] = os.getcwd()+"/"
        else:
            self.paths["CARDAMOM"] = raw_input("Enter CARDAMOM root directory: ")

        self.paths["projects"] = raw_input("Where will this project be saved? (provide full path or leave blank to use subdirectory \"projects\"): ")
        if self.paths["projects"] == "":
            self.paths["projects"] = self.paths["CARDAMOM"]+"/projects/"

        self.paths["library"] = raw_input("Where are source codes kept? (provide full path or leave blank to use subdirectory \"library\"): ")
        if self.paths["library"] == "":
            self.paths["library"] = self.paths["CARDAMOM"]+"/library/"

        usecluster = raw_input("Will you run this project on a cluster <y/n>? ")
        if usecluster == "y":
            self.paths["cluster_username"] = raw_input("Enter your username (leave blank for default): ")
            if self.paths["cluster_username"] == "":
                self.paths["cluster_username"] = "jexbraya"

            self.paths["cluster_address"] = raw_input("Enter address of cluster (leave blank for eddie): ")
            if self.paths["cluster_address"] == "":
                self.paths["cluster_address"] = "eddie.ecdf.ed.ac.uk"

            self.paths["cluster_directory"] = raw_input("Enter cluster working directory (full path or leave blank for default): ")
            if self.paths["cluster_directory"] == "":
                self.paths["cluster_directory"] = "/exports/work/scratch/jexbraya/"
            
        savedefault = raw_input("Save current paths as default ones for this machine <y/n>? ")
        
        if savedefault == 'y':
            f = open("default_paths_%s.pData" % hostname,"w")
            cPickle.dump(self.paths,f)
            f.close()


    def setup(self,latitude,longitude,drivers,observations,parprior,parpriorunc,otherprior,otherpriorunc):
        """
        This method is a wrapper to the function used to setup the project
        """
        # pixel data
     
        self.details["latitude"] = latitude
        self.details["longitude"] = longitude
        self.details["drivers"] = drivers
        self.details["observations"] = observations

        # MCMC specific values
        self.details["parprior"] = parprior
        self.details["parpriorunc"] = parpriorunc
        self.details["otherprior"] = otherprior
        self.details["otherpriorunc"] = otherpriorunc

        # technical options
        EDCs = raw_input("Use EDCs <y/n>? ")
        if EDCs == 'y':
            self.details["EDCs"] = True

        warning = 0
        #get the number of time steps
        if drivers.ndim == 2:
            self.details["no_pts"] = 1
            self.details["met_fields"] = drivers.shape[1]
            self.details["obs_fields"] = observations.shape[1]
            if drivers.shape[0] != observations.shape[0]:
                print " /!\ Warning: dimensions of drivers and observations arrays do not agree on the number of time steps /!\ "           
                warning += 1
            else:
                self.details["tsteps"] = drivers.shape[0]
                 
        else:
            self.details["no_pts"] = drivers.shape[0]
            self.details["met_fields"] = drivers.shape[2]
            self.details["obs_fields"] = observations.shape[2]
            if drivers.shape[1] != observations.shape[1]:
                print " /!\ Warning: dimensions of drivers and observations arrays do not agree on the number of time steps /!\ "
                warning += 1
            else:
                self.details["tsteps"] = drivers.shape[1]

            if drivers.shape[0] != latitude.shape[0]:
                print " /!\ Warning: dimensions of drivers and latitude arrays do not agree on the number of sites /!\ "  
                warning += 1
    
        #assign the number of met fields
        if warning == 0:
            print "Saving project before writing binary files...   ",
            self.save_project()
            print "OK"

            #write input data      
            self.createInput()
            #copy the source code
            self.backup_source()

            #create a directory for the executable
            if "exec" not in os.listdir(self.paths["projects"]+self.project_name):
                os.mkdir(self.paths["projects"]+self.project_name+"/exec")

            compile_code = raw_input("Compile local code <y/n>? ")
            if compile_code == "y":
                print "Compiling local code"
                self.backup_source()
                self.compile_local_code()
                

            cluster = raw_input("Send project files on the cluster <y/n>? ")
            
            if cluster == "y":
                self.details["cluster"] = True
                #copy data on cluster
                self.send_to_cluster()

            cluster = raw_input("Compile code on the cluster <y/n>? ")
            if cluster == "y":
                self.compile_cluster()

        else:
            print "Too many warning raised... please check input data according to previous messages" 




    def createInput(self):
        """
        This method writes the input files in the local directory
        """

        # Defines where the input files will be written        
        path2project = '/'.join([self.paths["projects"],self.project_name])
        path2data = path2project+"/data/"

        if "data" not in os.listdir(path2project):
            print "Directory \"%s\" not found... creating" % path2data
            os.mkdir(path2data)

        print "Now creating input data in \"%s\" for project \"%s\" with type \"%s\"" % (path2data,self.project_name,self.project_type)



        if self.project_type == "DALEC_GSI":
            modelid = 1
        elif self.project_type == "DALEC_CDEA_LU_FIRES":
            modelid = 2
        elif self.project_type == "DALEC_GSI_newalloc":
            modelid = 3 
        elif self.project_type == "DALEC_CDEA_LU_FIRES_LIU":
            modelid = 4 
        elif self.project_type == "DALEC_GSI_LIU":
            modelid = 5 
        elif self.project_type == "DALEC_CDEA_LU_FIRES_ET":
            modelid = 6
        elif self.project_type == "DALEC_CDEA_CWD_LU_FIRES":
            modelid = 7
        elif self.project_type == "DALEC_GSI_GLEAM":
            modelid = 8
        elif self.project_type == "DALEC_CDEA_GLEAM":
            modelid = 9
        elif self.project_type == "DALEC_CDEA_LU_FIRES_HBV":
            modelid = 10

        for ii in xrange(self.details["no_pts"]):
            #create an empty array to store data to be written
            towrite=np.zeros(300+self.details["tsteps"]*(self.details["met_fields"]+self.details["obs_fields"]),dtype="d")-9999.
                    
            #provide fixed values


            towrite[0]=modelid                       # pixel number           
            towrite[2]=self.details["tsteps"]        # no of time steps
            towrite[3]=self.details["met_fields"]    # no of met fields
            towrite[4]=self.details["obs_fields"]    # no of obs fields
            towrite[5]=self.details["EDCs"]          # use EDCs? 1: Yes / 0: No
           # towrite[6]=self.details["unc"]           # is uncertainty of obs provided
            #provide priors and met data

            #assign dummies to make code easier to read
            parprior = self.details["parprior"]
            parpriorunc = self.details["parpriorunc"]
            otherprior = self.details["otherprior"]
            otherpriorunc = self.details["otherpriorunc"]

            if self.details["no_pts"] == 1:           
                towrite[1]=self.details["latitude"]  # pixel latitude

                towrite[100:100+len(parprior)]=parprior
                towrite[150:150+len(parpriorunc)]=parpriorunc
                towrite[200:200+len(otherprior)]=otherprior
                towrite[250:250+len(otherpriorunc)]=otherpriorunc

                metobs=np.hstack([self.details["drivers"],self.details["observations"]])

            else:
                towrite[1]=self.details["latitude"][ii]  # pixel latitude

                towrite[100:100+len(parprior[ii])]=parprior[ii]
                towrite[150:150+len(parpriorunc[ii])]=parpriorunc[ii]
                towrite[200:200+len(otherprior[ii])]=otherprior[ii]
                towrite[250:250+len(otherpriorunc[ii])]=otherpriorunc[ii]

                metobs=np.hstack([self.details["drivers"][ii],self.details["observations"][ii]])

            towrite[300:]=metobs.flatten()

            #create binary data
            f=file(path2data+"%s_%05i.bin" % (self.project_name,ii+1),'wb')
            f.write(struct.pack(len(towrite)*'d',*towrite))
            f.close()

        print "Written a total of %i/%i input file" % (ii+1,self.details["no_pts"])

    def backup_source(self):
        """
        This method copies the source code in a project sub-directory
        """
        if "src" not in os.listdir(self.paths["projects"]+self.project_name):
            os.mkdir(self.paths["projects"]+self.project_name+"/src")
        os.system("cp -r %s/* %s/%s/src" % (self.paths["library"],self.paths["projects"],self.project_name))

        recompile_local = raw_input("Recompile local code <y/n>? ")
        if recompile_local=="y":
            self.compile_local_code()

    def update_source_cluster(self):
        """
        This method backs up the source code in a project sub-directory and sends it to 
        the cluster
        """
        self.backup_source()
        dest = self.paths["cluster_username"]+"@"+self.paths["cluster_address"]+":"+self.paths["cluster_directory"]+"/"+self.project_name
        os.system("scp -r %s/%s/src %s" % (self.paths["projects"],self.project_name,dest))

        recompile_cluster = raw_input("Recompile on cluster <y/n>? ")
        if recompile_cluster == "y":
            self.compile_cluster()

    def compile_local_code(self):
        """
        This method compiles the code and saves a backup
        """

        path2source=self.paths["projects"]+self.project_name+"/src/"
        path2include="%s/models/%s/likelihood/MODEL_LIKELIHOOD.c" % (path2source,self.project_type)
        path2exe = self.paths["projects"]+self.project_name+"/exec/"

        if "exec" not in os.listdir(self.paths["projects"]+self.project_name):
            os.mkdir(self.paths["projects"]+self.project_name+"/exec")
        #compile directly in good directory
        os.system("gcc -O3 %s/general/cardamom_main.c --include %s -o %s.exe -lm" % (path2source,path2include,path2exe+self.project_name))
                 

        
    def send_to_cluster(self):
        """
        This method sends the whole project to the cluster
        """ 

        dest=self.paths["cluster_username"]+"@"+self.paths["cluster_address"]+":"+self.paths["cluster_directory"]+"/"+self.project_name
        print "Copying binary files to remote destination \"%s\" " % dest
        os.system("ssh %s@%s mkdir %s/%s" % (self.paths["cluster_username"],self.paths["cluster_address"],self.paths["cluster_directory"],self.project_name))
 
        os.system("scp -r %s %s" % (self.paths["projects"]+self.project_name+"/data", dest))
        os.system("scp -r %s %s" % (self.paths["projects"]+self.project_name+"/src", dest))
        os.system("scp -r %s %s" % (self.paths["projects"]+self.project_name+"/exec", dest))

        print "Successfully copied the data on the cluster"

    def compile_cluster(self):
        """
        This method compiles the code on the cluster
        """
        
        path2cluster = self.paths["cluster_address"]
        path2source = self.paths["cluster_directory"]+self.project_name+"/src/"
        path2include = "%s/models/%s/likelihood/MODEL_LIKELIHOOD.c" % (path2source,self.project_type)
        path2exe = self.paths["cluster_directory"]+self.project_name+"/exec/"

        os.system("ssh %s 'gcc -O3 %s/general/cardamom_main.c --include %s -o %s.exe -lm'" % (path2cluster,path2source,path2include,path2exe+self.project_name))

    def resetup(self):
        """
        Redo the setup with attributes of the object
        """

        latitude      =  self.details["latitude"]
        longitude     =  self.details["longitude"]
        drivers       =  self.details["drivers"]
        observations  =  self.details["observations"]

        # MCMC specific values
        parprior      =  self.details["parprior"]
        parpriorunc   =  self.details["parpriorunc"]
        otherprior    =  self.details["otherprior"]
        otherpriorunc =  self.details["otherpriorunc"]

        print "Starting setup for project "+self.project_name

        self.setup(latitude,longitude,drivers,observations,parprior,parpriorunc,otherprior,otherpriorunc)


    def download_cluster(self, **kwargs):
        """
        This method downloads the results from the cluster
        """
        
        print "Preparing to download data from \"%s\"" % self.paths["cluster_address"]

        cluster_details = "%s@%s" % (self.paths["cluster_username"],self.paths["cluster_address"])

        src = "%s/%s/output/*" % (self.paths["cluster_directory"],self.project_name)
        
        #create the output folder 
        if "output" not in os.listdir("%s/%s" % (self.paths["projects"],self.project_name)):
            os.mkdir("%s/%s/output" % (self.paths["projects"],self.project_name))

        runlist = os.listdir("%s/%s/output/" % (self.paths["projects"],self.project_name))
        if "runid" in kwargs:
            runid = kwargs["runid"]
        else:
            if len(runlist) == 0:
                runid = 1
            else:
                runlist.sort()
                runid = int(runlist[-1].split("_")[-1])+1


        dst = "%s/%s/output/run_%03i/" % (self.paths["projects"],self.project_name,runid)
        if "run_%03i" % runid not in runlist:
            os.mkdir(dst)
            print "scp -r %s:%s %s" % (cluster_details,src,dst)
            os.system("scp -r %s:%s %s" % (cluster_details,src,dst))
        else:
            if len(os.listdir(dst)) != 0.:
                notempty = raw_input("Destination folder \"%s\" not empty... continue <y/n>?" % dst)
                if notempty == "y":
                    print "Downloading data in \"%s\"" % dst
                    os.system("scp -r %s:%s %s" % (cluster_details,src,dst))
            else:
                os.system("scp -r %s:%s %s" % (cluster_details,src,dst))
      

if __name__ == "__main__":

    # First associate a variable to the CARDAMOM object that holds / will hold the data
    FORMA=Dataset("../data/CARDAMOM_GLOBAL/FORMA_0.5_ERA-Interim_361x720.nc")
    d0 = dt.datetime(2005,12,19)
    dates=[]
    for dd in FORMA.variables["time"][1:-1]:
        dates.append(d0+dt.timedelta(float(dd)))

    simlength=len(dates)  

    degradation = FORMA.variables["FORMA"][:].data
    degradation[degradation==-9999.]=0.
    
    #test will use the pixel in the AMAZON with the most deforestation since 1 Jan 2006
    #we start from time step 1 rather than 0, :360 is to limit to south-western hemisphere
    selec=np.where(np.nansum(degradation[1:,180:,:360],0)==np.nanmax(np.nansum(degradation[1:,180:,:360],0)))
    
    latitude = FORMA.variables["latitude"][selec[0]+180][0]  
    longitude= FORMA.variables["longitude"][selec[1]][0]

    # met drivers are:
    # day of simulation, Tmn (C), Tmx (C), radation (MJ d-1), CO2 (ppmv), doy
    # newly added for LU: removal (fraction of AGB), timestep length
    met=np.zeros([simlength,8])-9999.

    #observations are:
    # nee; LAI; GPP
    obs=np.zeros([simlength,3])-9999.

    #get removal time series
    print "Extracting FORMA data...   ",     
    removal=FORMA.variables["FORMA"][2:,selec[0]+180,selec[1]].data #FORMA dates correspond to end of time step
    print "OK"

    # get time step length
    tsteps=[];dayofsim=[];dayofyear=[]
    # get met drivers
    tmn=[];tmx=[];rad=[]
    #ref day 0 for day of simulation
    dayref0=dt.datetime(2005,12,31)

    #add a last date 
    dates.append(dates[-1]+dt.timedelta(days=8))
    for ii in xrange(len(dates)-1):
        
        tsteps.append((dates[ii+1]-dates[ii]).days)    
        dayref=dt.datetime(dates[ii].year-1,12,31)
        dayofyear.append((dates[ii]-dayref).days+(tsteps[-1]-1)/2.)
        dayofsim.append((dates[ii]-dayref0).days+(tsteps[-1]-1)/2.)

        #get first and last day of timestep
        firstday=dates[ii];lastday=dates[ii+1]-dt.timedelta(days=1)
        if lastday.year==firstday.year:
            print "Getting ERA-Interim data for timestep starting on %02i/%02i/%4i" % (firstday.day,firstday.month,firstday.year)
            erai=Dataset("../data/ERA-Interim/0.5deg/ERA-Interim_%4i.nc" % firstday.year)
            
            #get index of first time step
            first_step = ((firstday-dayref).days-1)*2
            last_step  = first_step+(lastday-firstday).days*2
    
            #get coordinates of points
            latid=erai.variables["latitude"]==latitude

            #longitude starts from GMT in ERA-Interim data            
            lonid=erai.variables["longitude"]==longitude+360
          
            #extract the fields - mean daily rad in MJ m-2 d-1 (transformed from MJ m-2 per 12 hours)
            rad.append(1e-6*erai.variables['ssrd'][first_step:last_step,latid,lonid].sum()/(0.5*(last_step-first_step)))

            #for each day in the era data, get the min and max             
            tmpmn,tmpmx=[],[]
            for dd in range(tsteps[-1]-1):
                tmpmn.append(erai.variables['mn2t'][first_step+dd*2:first_step+(dd+1)*2,latid,lonid].min(0))            
                tmpmx.append(erai.variables['mx2t'][first_step+dd*2:first_step+(dd+1)*2,latid,lonid].max(0))
            tmn.append(np.mean(tmpmn)-273.15);tmx.append(np.mean(tmpmx)-273.15)

            erai.close()

    # add data to drivers' array
    met[:,0] = np.array(dayofsim); 
    met[:,1] = np.array(tmn)
    met[:,2] = np.array(tmx)   
    met[:,3] = np.array(rad);
    met[:,4] = np.zeros(simlength)+385.
    met[:,5] = np.array(dayofyear)
    met[:,6] = removal[:,0]
    met[:,7] = np.array(tsteps)

    #remove last date from list
    dates=dates[:-1]

    #get HWSD
    hwsd=Dataset("../data/CARDAMOM_GLOBAL/HWSD_OC_0.5deg_ERA-Interim_361x720.nc")
    latid=hwsd.variables["latitude"]==latitude
    lonid=hwsd.variables["longitude"]==longitude
    som=hwsd.variables["orgC"][:].data[latid,lonid][0]
    hwsd.close()

    #get woody biomass
    saatchi=Dataset("../data/CARDAMOM_GLOBAL/Saatchi_0.5deg_ERA-Interim_361x720.nc")
    latid=saatchi.variables["latitude"]==latitude
    lonid=saatchi.variables["longitude"]==longitude
    agb=saatchi.variables["AGB"][:].data[latid,lonid][0]*.5 # 0.5 for biomass to C
    saatchi.close()

    #get MODIS LAI observations
    

    path2modis="../data/MODIS/LAI/0.5deg/"
    lai=np.zeros(simlength)
    
    for ii,dat in enumerate(dates): 
        print "Getting MODIS LAI data on %02i/%02i/%4i" % (dat.day,dat.month,dat.year)
        #leap years with day after 28/2
        if dat.year %4 == 0 and dat.month >=3:
            tmpdat=dat+dt.timedelta(days=1)
            dateascii = "%4i-%02i-%02i" % (tmpdat.year,tmpdat.month,tmpdat.day)
        #leap year
        else:
            dateascii = "%4i-%02i-%02i" % (dat.year,dat.month,dat.day)
        #print ii, dateascii,
        modisfile = "MOD15A2_E_LAI_%s_rgb_720x360.CSV" % dateascii
        if modisfile in os.listdir(path2modis):
            #print ' ok'
            lai[ii] = np.loadtxt(path2modis+"MOD15A2_E_LAI_%s_rgb_720x360.CSV" % dateascii,delimiter=',')[selec[0]+180,selec[1]]
    print "OK"

    lai[lai==99999]=-9999.
    obs[:,1]=lai  

    #define priors

    parprior=np.zeros(50)-9999.;parpriorunc=np.zeros(50)-9999.
    #commonly used priors (from AAB code)
    parprior[1]=0.5;parpriorunc[1]=1.2
    parprior[9]=0.03;parpriorunc[9]=1.15 # temp_rate (Mahecha 2010)
    parprior[16]=70;parpriorunc[16]=2 #LMA - Kattge 2011
    #parprior[10]=20;parpriorunc[10]=1.5 #Ameriflux GPP & LAI assimilated data
    parprior[22]=som;parpriorunc[22]=1.5 
    parprior[20]=max(agb,100.);parpriorunc[20]=1.5 #saatchi agb for woody

    #
    parprior[4]=max(lai) / (max(lai)-min(lai[lai>0.]))
    parpriorunc[4]=1.2;

    #leaf onset & leaf fall dates priors
    if np.abs(latitude)>45:
        parprior[11]=120+365.25;parpriorunc[11]=1.1
        parprior[14]=270+365.25;parpriorunc[14]=1.1

    otherprior=np.zeros(50)-9999.;otherprior[0]=1
    otherpriorunc=np.zeros(50)-9999.;otherpriorunc[0]=2

### now that data is loaded - setup CARDAMOM

    
    






