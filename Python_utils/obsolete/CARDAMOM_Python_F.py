"""
19 September 2018 - JFE
Added:
- a method to read output files and return them directly as a pandas
DataFrame
- a method to rerun a pixel with parameter sets stored in an output file

17 April 2018 - JFE
This file is an adaptation of CARDAMOM_Python_C.py and contains the
CARDAMOM_F class to store CARDAMOM projects running with the Fortran
code
"""

#load modules
import os, cPickle, time, struct
import datetime as dt
import numpy as np
from netCDF4 import Dataset
import shutil, socket
import sys
import pandas as pd
from itertools import combinations
import xarray as xr
import re
import subprocess
import glob

class CARDAMOM_F(object):

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
                print("Project \"%s\" found" % self.project_name)
                self.load_project()
            else:
                print("Project directory found but file %s.pData could not be found" % self.project_name)
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
        types = ['ACM','DALEC_CDEA','DALEC_GSI_BUCKET','AT_DALEC', \
                'AT_DALEC_CROP','DALEC_CDEA_FR','DALEC_GSI_FR',
                'DALEC_GSI_FR_DBio','DALEC_GSI_MFOL_FR','DALEC_GSI_FR_LABILE', \
                'DALECN_GSI_FR','DALEC_GSI_DFOL_FR','DALEC_GSI_DFOL_FROOT_FR', \
                'DALEC_GSI_DFOL_LABILE_FR','DALECN_GSI_DFOL_LABILE_FR', \
                'DALECN_GSI_DFOL_LABILE_FROOT_FR','DALEC_GSI_DFOL_CWD_FR', \
                'DALECN_GSI_BUCKET','DALEC_CDEA_LU_FIRES']

        if 'model' in kwargs:
            self.model = kwargs['model']
            if self.model not in types:
                print("Warning - Unrecognized model version")
        else:
            print("Available model types")
            for mm,modname in enumerate(types):
                print('%02i - %s' % (mm,modname))
            projtype=raw_input("Choose model type from list above\n> ")
            if projtype != "":
                self.model = types[int(projtype)]
                self.modelid = int(projtype)
            else:
                print("Unknown project type")

        npools = (2, 6, 8, 6, 8, 6, 6, 10, 7, 7, 10, 6, 9, 9, 12, 12, 7, 8, 6)
        npars  = (20, 23, 40, 22, 34, 23, 33, 53, 36, 37, 49, 36, 46, 44, 57,
                    61, 38, 48, 23)
        nfluxes= (3, 16, 21, 16, 16, 18, 18, 28, 18, 20, 21, 18, 19, 19, 21,
                    21, 24, 25, 28)

        self.nfluxes = nfluxes[self.modelid]
        self.npars   = npars[self.modelid]
        self.npools  = npools[self.modelid]

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
            print("Paths found in project")
        else:
            print("Paths need to be defined")
            if "default_paths_%s.pData" % (hostname) in os.listdir(os.getcwd()):
                usedefault = raw_input("Use default paths <y/n>? ")
                if usedefault == 'y':
                    f = open("default_paths_%s.pData" % (hostname),"r")
                    self.paths = cPickle.load(f)
                    f.close()
                else:
                    self.define_paths()
            else:
                print("No default paths found... defining paths now")
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

        usehpc = raw_input("Will you run this project on a HPC <y/n>? ")
        if usehpc == "y":
            self.paths["hpc_username"] = raw_input("Enter your username (leave blank for default): ")
            if self.paths["hpc_username"] == "":
                self.paths["hpc_username"] = "jexbraya"

            self.paths["hpc_address"] = raw_input("Enter address of HPC (leave blank for eddie): ")
            if self.paths["hpc_address"] == "":
                self.paths["hpc_address"] = "eddie3.ecdf.ed.ac.uk"

            self.paths["hpc_directory"] = raw_input("Enter HPC working directory (full path or leave blank for default): ")
            if self.paths["hpc_directory"] == "":
                self.paths["hpc_directory"] = "/exports/csce/eddie/geos/groups/gcel/"

        savedefault = raw_input("Save current paths as default ones for this machine <y/n>? ")

        if savedefault == 'y':
            f = open("default_paths_%s.pData" % hostname,"w")
            cPickle.dump(self.paths,f)
            f.close()


    def setup(self, lat, lon, drivers, obs, obsunc, parprior, parpriorunc, \
    otherprior, otherpriorunc, EDCs = True, PFT = False):
        """
        This method is a wrapper to the function used to setup the project.
        Arguments are:
        - lat: a 1D array with pixels latitude
        - lon: a 1D array with pixels longitude
        - drivers: a 3D array [pixel, timestep, variable] containing drivers
        - obs: a 3D array containing time resolved observational constraints
            (-9999. for missing values)
        - obsunc: a 3D array containing the uncertainty in obs
        - parprior: a 3D array containing parameters prior values (-9999. for
            missing values)
        - parpriorunc: a 3D array containing the uncertainty around parameter
            priors (-9999. for missing values)
        - otherprior: a 3D array containing other priors (eg total vegetation C)
            (-9999. for missing values)
        - otherpriorunc: a 3D array containing the uncertainty of these other
            priors
        """

        # first store data
        self.lat            = lat
        self.lon            = lon
        self.drivers        = drivers
        self.obs            = obs
        self.obsunc         = obsunc

        # MCMC specific values
        self.parprior       = parprior
        self.parpriorunc    = parpriorunc
        self.otherprior     = otherprior
        self.otherpriorunc  = otherpriorunc

        # a couple of options
        self.edcs           = int(EDCs)
        self.pft            = int(PFT)

        #get the number of pixels, time steps, drivers and obs
        self.npts           = drivers.shape[0]
        self.nsteps         = drivers.shape[1]
        self.ndrivers       = drivers.shape[2]
        self.nobs           = obs.shape[2]

        print("Project data succesfully loaded, now saving...   ",)
        self.save_project()
        print("DONE")

    def createInput(self):
        """
        This method writes the input files in the local directory
        """

        # Defines where the input files will be written
        path2project = '/'.join([self.paths["projects"],self.project_name])
        path2data = path2project+"/data/"

        if "data" not in os.listdir(path2project):
            print("Directory \"%s\" not found... creating" % path2data)
            os.mkdir(path2data)

        print("Now creating input data in \"%s\" for project \"%s\" with type \"%s\"" % (path2data,self.project_name,self.model))

        for ii in xrange(self.npts):
            # create an empty array to store data to be written
            towrite=np.zeros(500+self.nsteps*(self.ndrivers+self.nobs*2),dtype="d")-9999.

            # provide fixed values in first 100
            towrite[0] = self.modelid
            towrite[1] = self.lat[ii]       # pixel number
            towrite[2] = self.nsteps        # no of time steps
            towrite[3] = self.ndrivers      # no of met fields
            towrite[4] = self.nobs*2        # no of obs fields * 2 with UC
            towrite[5] = self.edcs          # use EDCs? 1: Yes / 0: No
            towrite[6] = self.pft           # crop site?

            # now get priors and uncertainty - 100 for each following Luke
            towrite[100:100+self.parprior.shape[1]]=self.parprior[ii]
            towrite[200:200+self.parpriorunc.shape[1]]=self.parpriorunc[ii]
            towrite[300:300+self.otherprior.shape[1]]=self.otherprior[ii]
            towrite[400:400+self.otherpriorunc.shape[1]]=self.otherpriorunc[ii]

            #loop over time steps to extract drivers, obs and uncertainty
            metobs = np.zeros([self.nsteps,self.ndrivers+2*self.nobs])
            for nn in xrange(self.nsteps):
                # store the met drivers of the corresponding time step
                metobs[nn,:self.ndrivers] = self.drivers[ii,nn]
                # store observations and corresponding uncertainty
                for oo in xrange(self.nobs):
                    metobs[nn,self.ndrivers+(oo*2)]   = self.obs[ii,nn,oo]
                    metobs[nn,self.ndrivers+(oo*2)+1] = self.obsunc[ii,nn,oo]
            towrite[500:]=metobs.flatten()

            #create binary data
            f=file(path2data+"%s_%05i.bin" % (self.project_name,ii+1),'wb')
            f.write(struct.pack(len(towrite)*'d',*towrite))
            f.close()

        print("Written input file for %i sites" % (self.npts))

    def backup_source(self):
        """
        This method copies the source code in a project sub-directory.
        """

        if "src" not in os.listdir(self.paths["projects"]+self.project_name):
            os.mkdir(self.paths["projects"]+self.project_name+"/src")
        os.system("cp -r %s/* %s/%s/src/" % (self.paths["library"],self.paths["projects"],self.project_name))

        print('Copied source code from repository %s to local project %s/%s/src/' % (self.paths['library'],self.paths["projects"],self.project_name))


    def update_source_hpc(self):
        """
        This method backs up the source code in a project sub-directory and sends it to
        the HPC
        """
        self.backup_source()
        dest = self.paths["hpc_username"]+"@"+self.paths["hpc_address"]+":"+self.paths["hpc_directory"]+"/"+self.project_name
        os.system("scp -r %s/%s/src %s" % (self.paths["projects"],self.project_name,dest))

        recompile_hpc = raw_input("Recompile on HPC <y/n>? ")
        if recompile_hpc == "y":
            self.compile_hpc()

    def compile_local(self, compiler ='ifort', flags ='-O2'):
        """
        This method compiles the code locally, using ifort
        with optimization flags by default.
        """

        # create folder to hold exec
        if "exec" not in os.listdir(self.paths["projects"]+self.project_name):
            os.mkdir(self.paths["projects"]+self.project_name+"/exec")

        #define path to library, model version and path to exec
        path2lib = '%s/%s/src/' % (self.paths["projects"],self.project_name)
        model = self.model
        #executable bears the name of the project
        path2exe = '%s/%s/exec/%s.exe' % (self.paths["projects"],self.project_name,self.project_name)

        #set compiler and options in the command
        cmd = '%s %s' % (compiler, flags)

        # copied file order from Luke's
        cmd += ' %s/misc/math_functions.f90 %s/misc/oksofar.f90' % (path2lib,path2lib)  # helpful functions
        cmd += ' %s/model/%s/src/%s.f90' % (path2lib,model,model)                       # the model itself
        if model+'_CROP.f90' in os.listdir('%s/model/%s/src/' % (path2lib,model)):
            cmd += ' %s/model/%s/src/%s_CROP.f90' % (path2lib,model,model)              # the crop model if it exists
        cmd += ' %s/general/cardamom_structures.f90' % (path2lib)                       # structure definition
        cmd += ' %s/method/MHMCMC/MCMC_FUN/MHMCMC_STRUCTURES.f90' % (path2lib)          # MCMC specific structures
        cmd += ' %s/model/%s/src/%s_PARS.f90' % (path2lib,model,model)                  # the file holding boundary values of parameters
        cmd += ' %s/general/cardamom_io.f90' % (path2lib)                               # the file with the IO functions
        cmd += ' %s/method/MHMCMC/MCMC_FUN/MHMCMC.f90' % (path2lib)                     # the actual MCMC function
        cmd += ' %s/model/%s/likelihood/MODEL_LIKELIHOOD.f90' % (path2lib,model)        # the likelihood files / includes EDCs
        cmd += ' %s/general/cardamom_main.f90' % (path2lib)                             # the main file
        cmd += ' -o %s' % path2exe

        print(cmd)
        os.system(cmd)

    def compile_f2py(self, fcompiler = 'intelem', opt ='-O2'):
        """
        This method compiles the f2py version of the source code that should exist in the
        same folder as the code
        """
        path2lib = '%s/%s/src/' % (self.paths["projects"],self.project_name)
        model = self.model
        #first check that there is an f2py version - by convention it should be name <model>_f2py.f90
        if model+'_f2py.f90' in os.listdir('%s/model/%s/src/' % (path2lib,model)):
            path2src = '%s/model/%s/src/%s_f2py.f90' % (path2lib,model,model)

            cmd = 'f2py -c -m f2py_model --fcompiler="%s" --opt="%s" %s' % (fcompiler,opt,path2src)

            print(cmd)
            # compile and copy
            os.system(cmd)
            os.system('mv f2py_model.so %s/%s/exec/' % (self.paths["projects"],self.project_name) )
        else:
            print('No f2py source file found')


    def send_to_hpc(self):
        """
        This method sends the necessary parts of the project to the hpc
        """

        dest=self.paths["hpc_directory"]+"/"+self.paths['hpc_username']+'/'+self.project_name
        print("Copying binary files to remote destination \"%s\" " % dest)
        #create the directory
        os.system("ssh %s@%s mkdir %s" % (self.paths["hpc_username"],self.paths["hpc_address"],dest))
        #create directory for the output
        os.system("ssh %s@%s mkdir %s/output" % (self.paths["hpc_username"],self.paths["hpc_address"],dest))
        #create the directory to store log files
        os.system("ssh %s@%s mkdir %s/log" % (self.paths["hpc_username"],self.paths["hpc_address"],dest))
        #copy data, source and exec
        os.system("scp -r %s %s@%s:%s" % (self.paths["projects"]+self.project_name+"/data/", self.paths["hpc_username"],self.paths["hpc_address"],dest))
        os.system("scp -r %s %s@%s:%s" % (self.paths["projects"]+self.project_name+"/src/", self.paths["hpc_username"],self.paths["hpc_address"],dest))
        os.system("scp -r %s %s@%s:%s" % (self.paths["projects"]+self.project_name+"/exec/", self.paths["hpc_username"],self.paths["hpc_address"],dest))

        print("Successfully copied the data on the hpc")

    def compile_hpc(self,compiler ='ifort', flags ='-O2'):
        """
        This method compiles the code on the hpc
        """

        path2lib = '%s/%s/src/' % (self.paths["projects"],self.project_name)
        path2lib_hpc = '%s/%s/%s/src/' % (self.paths['hpc_directory'],self.paths['hpc_username'],self.project_name)
        model = self.model
        #executable bears the name of the project
        path2exe_hpc = '%s/%s/%s/exec/%s.exe' % (self.paths["hpc_directory"],self.paths['hpc_username'],self.project_name,self.project_name)

        #set compiler and options in the command
        cmd = '%s %s' % (compiler, flags)

        # copied file order from Luke's
        cmd += ' %s/misc/math_functions.f90 %s/misc/oksofar.f90' % (path2lib_hpc,path2lib_hpc)  # helpful functions
        cmd += ' %s/model/%s/src/%s.f90' % (path2lib_hpc,model,model)                           # the model itself
        if model+'_CROP.f90' in os.listdir('%s/model/%s/src/' % (path2lib,model)):
            cmd += ' %s/model/%s/src/%s_CROP.f90' % (path2lib_hpc,model,model)                  # the crop model if it exists
        cmd += ' %s/general/cardamom_structures.f90' % (path2lib_hpc)                           # structure definition
        cmd += ' %s/method/MHMCMC/MCMC_FUN/MHMCMC_STRUCTURES.f90' % (path2lib_hpc)              # MCMC specific structures
        cmd += ' %s/model/%s/src/%s_PARS.f90' % (path2lib_hpc,model,model)                      # the file holding boundary values of parameters
        cmd += ' %s/general/cardamom_io.f90' % (path2lib_hpc)                                   # the file with the IO functions
        cmd += ' %s/method/MHMCMC/MCMC_FUN/MHMCMC.f90' % (path2lib_hpc)                         # the actual MCMC function
        cmd += ' %s/model/%s/likelihood/MODEL_LIKELIHOOD.f90' % (path2lib_hpc,model)            # the likelihood files / includes EDCs
        cmd += ' %s/general/cardamom_main.f90' % (path2lib_hpc)                                 # the main file
        cmd += ' -o %s' % path2exe_hpc

        if compiler == 'ifort':
            cmd = 'module load intel && '+cmd

        print(cmd)
        os.system("ssh %s@%s '%s'" % (self.paths['hpc_username'],self.paths["hpc_address"],cmd))

    def resetup(self):
        """
        Redo the setup with attributes of the object
        """

        lat             =   self.lat
        lon             =   self.lon
        drivers         =   self.drivers
        obs             =   self.obs
        obsunc          =   self.obsunc

        # MCMC specific values
        parprior        =   self.parprior
        parpriorunc     =   self.parpriorunc
        otherprior      =   self.otherprior
        otherpriorunc   =   self.otherpriorunc

        edcs            =   self.edcs
        pft             =   self.pft

        print("Starting setup for project "+self.project_name)

        self.setup(lat,lon,drivers,obs,obsunc,parprior,parpriorunc,otherprior,otherpriorunc,edcs,pft)

    def submit_hpc(self, nchains=3, nsims = 1000000, nout = 1000, runtime = '06:00:00'):
        """
        This method submits the jobs on the cluster, with string of options.
        """

        #basic logic is to create a submittable script, upload to eddie and run it
        #this way the script can get the SGE_TASK_ID env variable to run on
        #arrays of jobs

        #get different paths
        path2proj = "%s/%s/%s/" % (self.paths["hpc_directory"],self.paths["hpc_username"],self.project_name)
        path2log  = path2proj+'log/'
        path2exec = path2proj+'exec/%s.exe' % (self.project_name)
        hpc_username = self.paths['hpc_username']
        hpc_address = self.paths['hpc_address']

        #define header
        header = '#!/usr/bin/python\n'
        header+= '#$ -l h_rt=%s\n' % (runtime)
        header+= '#$ -o %s\n' % path2log
        header+= '#$ -e %s\n\n' % path2log
        header+= 'import os\n'
        header+= 'taskid=int(os.environ["SGE_TASK_ID"])\n\n' #get the task id of element in array

        for ch in range(1,1+nchains):
            print("Writing script for chain ", ch)
            script = header
            script+= 'os.system("rm %s/output/%s_%%05i_%i_*" %% (taskid))\n\n' % (path2proj,self.project_name,ch)
            script+= 'os.system("%s/exec/%s.exe ' % (path2proj,self.project_name)
            script+= '%s/data/%s_%%05i.bin ' % (path2proj,self.project_name)
            script+= '%s/output/%s_%%05i_%i_ %i 0 %i" %% (taskid,taskid))\n\n' % (path2proj,self.project_name,ch,nsims,nout)

            #write the script
            scriptname = 'script_%s_%i.sub' % (self.project_name,ch)
            file(scriptname,'w').write(script)

            #upload the script in the home directory of the hpc
            os.system('scp %s %s@%s:' % (scriptname,hpc_username,hpc_address))
            #submit the script
            os.system('ssh %s@%s qsub -t %i-%i -N %s_%i %s' % (hpc_username,hpc_address,1,min(self.npts,10000),self.project_name,ch,scriptname))
            #submit an additional job if more than 10000 pts
            if self.npts > 10000:
                os.system('ssh %s@%s qsub -t %i-%i -N %s_%i_ %s' % (hpc_username,hpc_address,10001,self.npts,self.project_name,ch,scriptname))


    def download_hpc(self, **kwargs):
        """
        This method downloads the results from the HPC
        """

        print("Preparing to download data from \"%s\"" % self.paths["hpc_address"])

        hpc_details = "%s@%s" % (self.paths["hpc_username"],self.paths["hpc_address"])

        src = "%s/%s/%s/output/*" % (self.paths["hpc_directory"],self.paths["hpc_username"],self.project_name)

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
            print("scp -r %s:%s %s" % (hpc_details,src,dst))
            os.system("scp -r %s:%s %s" % (hpc_details,src,dst))
        else:
            if len(os.listdir(dst)) != 0.:
                notempty = raw_input("Destination folder \"%s\" not empty... continue <y/n>?" % dst)
                if notempty == "y":
                    print("Downloading data in \"%s\"" % dst)
                    os.system("scp -r %s:%s %s" % (hpc_details,src,dst))
            else:
                os.system("scp -r %s:%s %s" % (hpc_details,src,dst))

    def read_output(self, run=1, pixel=1, chain=1, burnin = 0.5):
        """
        This method reads output files for a specified pixel and run and returns
        it in a pandas DataFrame. Output files have to be stored under
        self.
        The burnin option indicate the fraction of the output file to be
        disregarded as the MCMC burnin. If burnin < 1., the method disregards a
        fraction of the file, else it ignores the number of parameter sets specified.
        """

        #first define the path to the file
        path = self.paths['projects'] + self.project_name + '/output/run_%03i/' % (run)
        path += self.project_name+'_%05i_%i_PARS' % (pixel,chain)

        #get the number of parameters from the object
        npars = self.npars
        #JFE checking if the are output file exists
        if len(glob.glob(path)) > 0:
            fpars=file(path,'rb')
            content=fpars.read()

            #if file is not empty, extract the parameter sets
            if content != '':
                parsets=np.array(struct.unpack((len(content)/8)*'d',content)).reshape([len(content)/(8*(npars+1)),npars+1])

                if burnin < 1:
                    keep = int(parsets.shape[0]*burnin)
                elif burnin > 1:
                    keep = burnin

                parsets = parsets[keep:]

            # if file is empty, return empty array of parameters
            else:
                parsets = np.zeros([500,npars+1])-9999.

        #no output file, return empty array of parameters
        else:
            parsets = np.zeros([500,npars+1])-9999.

        #return the parameter sets as a pandas dataframe
        return pd.DataFrame(parsets)

    def rerun_pixel(self, run=1, pixel=1, chains=1, thresh=1.2, burnin = 0.5):
        """
        This method reruns the output for a defined pixel using parameter sets
        stored in the corresponding file. The burnin option indicates the fraction
        of the parameter sets that should be disregarded as burn in of the MCMC
        """

        # first reads the parameters - either a single file or get_converged_parameters
        if type(chains) == int:
            chains=[chains]

        parsets = self.get_converged_parameters(run=run,pixel=pixel,chains=chains,
                                                burnin=burnin,thresh=thresh,output_pars=True)
        #read_output and get_converged_parameters outpustrip()t pandas, now transform
        #to array while removing the last column
        parsets = parsets.values[:,:-1]
        nparsets= parsets.shape[0]
        #add path to f2py compiled module to sys
        path2f2py = "%s/%s/exec/" % (self.paths["projects"],self.project_name)
        if path2f2py not in sys.path:
            sys.path.append(path2f2py)
        # import the f2py module
        if 'f2py_model.so' in os.listdir(path2f2py):
            from f2py_model import carbon_model
    #        print('Loaded f2py_model.so')
        else:
            print('Module f2py_model.so not found, use method compile_f2py')
            return None

        #### now entering the actual running phase
        # first get the pixel's general info -
        # /!\ pixels are indexed from 0, but files start at 1
        pixid   = pixel-1
        #create / load invariant parameters
        start   = 1
        finish  = self.nsteps
        met     = self.drivers[pixid].T
        #create the deltat time series
        deltat  = np.append([self.drivers[pixid,0,0]],np.diff(self.drivers[pixid,:,0]))
        lat     = self.lat[pixid]
        #now create arrays that will store the output
        lai     = np.zeros([nparsets,self.nsteps])-9999.
        nee     = np.zeros([nparsets,self.nsteps])-9999.
        fluxes  = np.zeros([nparsets,self.nsteps,self.nfluxes])-9999.
        pools   = np.zeros([nparsets,self.nsteps+1,self.npools])-9999.
        gpp     = np.zeros([nparsets,self.nsteps])-9999.

        #now loop over the parameter sets
        for pp,pars in enumerate(parsets):
            #do not do anything if all parameters are -9999. => missing data
            if pars.size != (pars == -9999.).sum():
                lai[pp],nee[pp],fluxes[pp],pools[pp],gpp[pp] = carbon_model(start,finish,
                                                                        met,pars,deltat,
                                                                        lat,lai[pp],
                                                                        nee[pp],fluxes[pp],
                                                                        pools[pp],gpp[pp])

        return lai,nee,fluxes,pools,gpp

    def GR(self,chains):
        """
        This method returns the Gelman-Rubin convergence criterion between
        parameter distributions stored in chains
        """
        n = len(chains[0]) # number of draws per chain
        m = len(chains)    # number of chains

        # calculate the within chain variance W: mean of the variance of each chain
        W = 0.
        for cc in range(m):

            sj2 = np.var(chains[cc],ddof=1)
            W += sj2
        W = W / m

        # calculate the between chain variance B: variance of chains means
        theta = 0.
        for cc in range(m):
            theta += np.mean(chains[cc])
        theta = theta*1./m

        B = 0.
        for cc in range(m):
            B += (np.mean(chains[cc])-theta)**2
        B = B*n/(m+1)

        # estimate variance
        vartheta = (1.-1./n)*W+B/n

        return vartheta/W

    def get_converged_parameters(self, run=1, pixel=1, chains=[1,2,3], burnin = 0.5,
                            thresh = 1.2, output_pars = True):
        """
        This method checks whether convergence has been reached between specified
        chains for a specified run, pixel according to a threshold for the
        Gelman-Rubin convergence metric.
        """

        if type(chains) == int:
            chains = [chains]

        # define the number of chains
        nchains = len(chains)

        #define the possible combinations of 2+ chains
        chain_combinations = []
        for n in range(nchains,1,-1):
            chain_combinations += list(combinations(chains,n))

        conv = np.zeros(len(chain_combinations),dtype='bool')
        #loop over the possible combinations until criterion met
        for cc,combi in enumerate(chain_combinations):
            pars = []
            for ch in combi:
                pars.append(self.read_output(run=run,pixel=pixel,chain=ch,burnin=burnin).values[:,-1])
            #JFE concat dataframes to test whether data's missing
            pars1d = np.concatenate(pars)
            if pars1d.size != (pars1d==-9999.).sum():
                crit = self.GR(pars)
                if crit < thresh:
                    conv[cc] = True
            else:
                conv[cc] = True

        zipped_res = zip(chain_combinations,conv)
        if not output_pars:
            #output the results of the convergence test
            return zipped_res
        else:
            #start a counter to check whether one of the chain combinations has converged
            conv = 0
            #loop over inverted combinations list, to be sure the converged combination with max no of chains is last
            for combi in zipped_res[::-1]:
                #if combination as converged get the data
                if combi[1] == True:
                    conv+=1
                    for ch,chainid in enumerate(combi[0]):
                        if ch == 0:
                            pars = self.read_output(run=run,pixel=pixel,chain=chainid,burnin=burnin)
                        else:
                            pars = pd.concat([pars,self.read_output(run=run,pixel=pixel,chain=chainid,burnin=burnin)])
            # if no combination has converged, return empty parameter sets
            if conv == 0:
                pars = pd.DataFrame(np.zeros([500,self.npars])-9999)
                pars = pd.concat([self.read_output(run=run,pixel=pixel,chain=ch,burnin=burnin) for ch in chains],axis='rows')
            pars.index = np.arange(pars.shape[0])
            return pars


    def get_parameter_maps(self,run=1,chains=[1,2,3], burnin = 0.5,thresh=1.2,
                        percentiles=[2.5,5.,25.,50.,75.,95.,97.5],save=True):
        """
        This method returns a map of the parameters as xarray.
        """

        #first define the grid by extracting the lat and lon res
        dlon = np.abs(self.lon.max()-self.lon); reslon = dlon[dlon!=0.].min()
        dlat = np.abs(self.lat.max()-self.lat); reslat = dlat[dlat!=0.].min()
        #define the lat/lon of the grid...
        longrid = np.arange(self.lon.min(),self.lon.max()+reslon,reslon)
        latgrid = np.arange(self.lat.max(),self.lat.min()-reslat,-reslat)
        #create the xarray.Dataset that will hold all the parameters
        #specify data_vars
        dummy = np.zeros([len(percentiles),latgrid.size,longrid.size])-9999.
        data_vars = {}
        for ii in range(self.npars):
            data_vars['p%02i' % (ii)] = (['percentile','lat','lon'],dummy.copy(),{'_FillValue':-9999.})
        data_vars['likelihood'] = (['percentile','lat','lon'],dummy.copy(),{'_FillValue':-9999.})
        coords = {'percentile': (['percentile'],percentiles,{'units':'%'}),
                     'lat': (['lat'],latgrid,{'units':'degrees_north'}),
                     'lon': (['lon'],longrid,{'units':'degrees_east'})}

        self.parameter_maps = xr.Dataset(data_vars=data_vars,coords=coords)
        #iterate over parameters
        for pp in np.arange(self.npts):
            print '\rGetting parameter sets for pixel %05i / %05i' % (pp+1,self.npts),
            latid = np.where(latgrid==self.lat[pp])[0][0]
            lonid = np.where(longrid==self.lon[pp])[0][0]
            pars = self.get_converged_parameters(run=run,pixel=pp+1,chains=chains,
                                            thresh=thresh, output_pars=True)
            pct = np.percentile(pars,percentiles,axis=0)
            for ii in range(self.npars):
                self.parameter_maps['p%02i' % (ii)][:,latid,lonid] = pct[:,ii]
            self.parameter_maps['likelihood'][:,latid,lonid] = pct[:,-1]

        print 'Done extracting parameters'
        #now check whether it has to be save
        if save:
            dst = self.paths["projects"]+self.project_name+'/post/'+self.project_name+'_run_%03i_parameters.nc' % run
            print 'Saving parameter maps as netCDF at '+dst
            if 'post' not in os.listdir(self.paths["projects"]+self.project_name):
                os.mkdir(self.paths["projects"]+self.project_name+'/post/')

            self.parameter_maps.to_netcdf(dst,'w')

    def rerun_all_pixels(self, run=1, chains=[1,2,3], thresh=1.2, burnin = 0.5,
                            percentiles=[2.5,5.,25.,50.,75.,95.,97.5],
                            keep_fluxes=['0','2','0-2','12+13','2+12+13','2+12+13-0','16','2+12+13+16-0'],
                            keep_pools=['0','1','2','3','4','5','0+1','0+1+2+3','4+5','0+1+2+3+4+5']):
        """
        This methods runs all pixels and stores the output in two large arrays.
        If other_fluxes is defined, the method returns specified fluxes with possible
        combinations.
        """

        #create a function to understand the content of keep_fluxes
        def compound(method,array):
            meth = method.replace(' ','')
            indexes = re.split('[+-]',meth)
            if len(indexes) == 1:
                return array[:,:,int(meth)]
            else:
                operators = re.findall('[+-]',meth)
                operators = [np.add if ops == '+' else np.subtract for ops in operators]
                #print operators
                #if only missing vals in first flux/pool => return -9999.
                if (array[:,:,int(indexes[0])] == -9999.).sum() == (array[:,:,int(indexes[0])]).size:
                    tmp = np.zeros(array.shape[:2])-9999.
                else:
                    for oo, ops in enumerate(operators):
                        if oo == 0:
                            tmp = operators[0](array[:,:,int(indexes[0])],array[:,:,int(indexes[1])])
                        else:
                            tmp = operators[oo](tmp,array[:,:,int(indexes[oo+1])])
                return tmp

        #create arrays to store output
        if keep_fluxes == 'all':
            all_fluxes = np.zeros([self.npts,len(percentiles),self.nsteps,self.nfluxes])-9999.
        else:
            all_fluxes = np.zeros([self.npts,len(percentiles),self.nsteps,len(keep_fluxes)])-9999.

        if keep_pools == 'all':
            all_pools   = np.zeros([self.npts,len(percentiles),self.nsteps+1,self.npools])-9999.
        else:
            all_pools   = np.zeros([self.npts,len(percentiles),self.nsteps+1,len(keep_pools)])-9999.

        #loop over pixels rerun the model and store output
        for pixel in range(1,self.npts+1):
            print '\rRerunning pixel %05i / %05i' % (pixel,self.npts),
            lai,nee,fluxes,pools,gpp = self.rerun_pixel(run=run,pixel=pixel,chains=chains,thresh=thresh,burnin=burnin)
            if keep_fluxes == 'all':
                all_fluxes[pixel-1] = np.percentile(fluxes,percentiles,axis=0)
            else:
                for ii,method in enumerate(keep_fluxes):
                    all_fluxes[pixel-1,:,:,ii] = np.percentile(compound(method,fluxes),percentiles,axis=0)

            if keep_pools == 'all':
                all_pools[pixel-1] = np.percentile(fluxes,percentiles,axis=0)
            else:
                for ii, method in enumerate(keep_pools):
                    all_pools[pixel-1,:,:,ii] = np.percentile(compound(method,pools),percentiles,axis=0)

        #save the rerun fluxes and pools
        self.fluxes_rerun = all_fluxes
        self.pools_rerun = all_pools
        #save the pctiles and run id
        self.percentiles_rerun = percentiles
        self.rerun_id = run

    def save_fluxes_to_netcdf(self,start='2000-01-01',freq='M',fluxnames=['gpp','ra','npp','rh','reco','nee','fire','nbe'],dest=None):
        """
        This method outputs the fluxes array created with rerun_all_pixels into
        a netcdf file.
        """
        #first define the grid by extracting the lat and lon res
        dlon = np.abs(self.lon.max()-self.lon); reslon = dlon[dlon!=0.].min()
        dlat = np.abs(self.lat.max()-self.lat); reslat = dlat[dlat!=0.].min()
        #define the lat/lon of the grid...
        longrid = np.arange(self.lon.min(),self.lon.max()+reslon,reslon)
        latgrid = np.arange(self.lat.max(),self.lat.min()-reslat,-reslat)
        #define the date_range
        dates = pd.date_range(start=start,freq=freq,periods=self.nsteps)

        dummy = np.zeros([len(dates),len(self.percentiles_rerun),latgrid.size,longrid.size])-9999.

        attrs = {'_FillValue':-9999.,'units':'g C m-2 d-1'}
        data_vars = {}
        for ii,fluxname in enumerate(fluxnames):
            data_vars[fluxname] = (['time','percentile','lat','lon'],dummy.copy(),attrs)

        coords = {'time': (['time'],dates),
                  'percentile': (['percentile'],self.percentiles_rerun,{'units':'%'}),
                  'lat': (['lat'],latgrid,{'units':'degrees_north'}),
                  'lon': (['lon'],longrid,{'units':'degrees_east'})}

        self.fluxes_maps = xr.Dataset(data_vars=data_vars,coords=coords)

        #loop over pixels to save in a map
        for pp in np.arange(self.npts):
            #get the pixel indexes on the lat/lon grid
            latid = np.where(latgrid==self.lat[pp])[0][0]
            lonid = np.where(longrid==self.lon[pp])[0][0]
            #loop over the fluxes
            for ii in range(self.fluxes_rerun.shape[-1]):
                fluxname = fluxnames[ii]
                #transpose as the output array saved percentile dim before time dim
                self.fluxes_maps[fluxname][:,:,latid,lonid] = self.fluxes_rerun[pp,:,:,ii].T

        if dest == None:
            dst = self.paths["projects"]+self.project_name+'/post/'+self.project_name+'_run_%03i_fluxes.nc' % self.rerun_id
        else:
            dst = dest

        self.fluxes_maps.to_netcdf(dst,'w')

    def save_pools_to_netcdf(self,start='2000-01-01',freq='M',poolnames=['Clab','Cleaf','Croot','Cwood','Clitter','Csom','Cfol','Cveg','Cdom','Ctot'],dest=None):
        """
        This method outputs the pools array created with rerun_all_pixels into
        a netcdf file.
        """
        #first define the grid by extracting the lat and lon res
        dlon = np.abs(self.lon.max()-self.lon); reslon = dlon[dlon!=0.].min()
        dlat = np.abs(self.lat.max()-self.lat); reslat = dlat[dlat!=0.].min()
        #define the lat/lon of the grid...
        longrid = np.arange(self.lon.min(),self.lon.max()+reslon,reslon)
        latgrid = np.arange(self.lat.max(),self.lat.min()-reslat,-reslat)
        #define the date_range
        dates = pd.date_range(start=start,freq=freq,periods=self.nsteps)

        dummy = np.zeros([len(dates),len(self.percentiles_rerun),latgrid.size,longrid.size])-9999.

        attrs = {'_FillValue':-9999.,'units':'g C m-2'}
        data_vars = {}
        for ii,poolname in enumerate(poolnames):
            data_vars[poolname] = (['time','percentile','lat','lon'],dummy.copy(),attrs)

        coords = {'time': (['time'],dates),
                  'percentile': (['percentile'],self.percentiles_rerun,{'units':'%'}),
                  'lat': (['lat'],latgrid,{'units':'degrees_north'}),
                  'lon': (['lon'],longrid,{'units':'degrees_east'})}

        self.pools_maps = xr.Dataset(data_vars=data_vars,coords=coords)

        #loop over pixels to save in a map
        for pp in np.arange(self.npts):
            #get the pixel indexes on the lat/lon grid
            latid = np.where(latgrid==self.lat[pp])[0][0]
            lonid = np.where(longrid==self.lon[pp])[0][0]
            #loop over the pools
            for ii in range(self.pools_rerun.shape[-1]):
                poolname = poolnames[ii]
                #transpose as the output array saved percentile dim before time dim
                #for pools, remove the first time step
                self.pools_maps[poolname][:,:,latid,lonid] = self.pools_rerun[pp,:,1:,ii].T

        if dest == None:
            dst = self.paths["projects"]+self.project_name+'/post/'+self.project_name+'_run_%03i_pools.nc' % self.rerun_id
        else:
            dst = dest

        self.pools_maps.to_netcdf(dst,'w')

if __name__ == "__main__":


    # example project using data in drivers.csv file
    data                = pd.read_csv('drivers.csv',parse_dates = True, index_col = 'date')

    # get drivers and reshape as 3D array - only 1st year
    drivers             = np.zeros([data.shape[0],8])
    drivers[:,:-2]      = data.get_values()[:,:-1]
    drivers             = np.expand_dims(drivers,0)

    #get observations and reshape as 3D array
    obs                 = np.zeros([drivers.shape[1],17])-9999.
    obs[:,1]            = data.LAI.get_values()
    obs                 = np.expand_dims(obs,0)

    #define observations uncertainty and reshape as 3D array
    obsunc              = np.zeros(obs.shape)-9999.
    obsunc[obs!=-9999.] = .5 #uncertainty as fraction of obs

    # define lat / lon
    lat                 = np.array([-12.75])
    lon                 = np.array([131.25])

    # define parpriors for allocation to autotrophic respiration
    # and canopy efficiency
    parprior            = np.zeros([1,100])-9999.
    parprior[:,1]       = 0.5
    parprior[:,10]      = 17.5

    # define their uncertainty
    parpriorunc         = np.zeros([1,100])-9999.
    parpriorunc[:,1]    = 1.2
    parpriorunc[:,10]   = 1.2

    #define otherprior and their uncertainty
    otherprior          = np.zeros([1,100])-9999.
    otherpriorunc       = np.zeros([1,100])-9999.

    #create / load project
    prj = CARDAMOM_F(project_name='fortran_test_cdea_lu_fires')

    # setup / store data in object
    prj.setup(lat,lon,drivers,obs,obsunc,parprior,parpriorunc,otherprior,otherpriorunc)

    # write binary files
    prj.createInput()

    # backup source and compile
    prj.backup_source()
    prj.compile_local(flags='-O2')
