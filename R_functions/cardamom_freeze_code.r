
###
## Function to freeze CARAMOM
###

# This function is based on an original Matlab function development by A. A. Bloom (UoE, now at the Jet Propulsion Laboratory).
# Translation to R and subsequent modifications by T. L Smallman (t.l.smallman@ed.ac.uk, UoE).

cardamom_freeze_code<- function (PROJECT) {

      #function CARDAMOM_FREEZE_CODE(PROJECT)
      #PROJECT = PROJECT structure provided with CARDAMOM
      #if PROJECT is empty, then your local copy will still be made
      #in the CARDAMOM_LOCAL directory.

      print('***********CADRAMOM_FREEZE_CODE***************')
      print('"Freezing" ALL CARDAMOM code :)')
      print('This function ensures that your current       ')
      print('CARDAMOM code is saved in its present state.  ')
      print('If you wish to retrieve any information from  ')
      print('the frozen CARDAMOM code, make sure to unfreeze')
      print('code by using CARDAMOM_UNFREEZE_CODE')
      print('----------------------------------------------')
      print('NOTE: CADRAMOM_FREEZE_CODE will be run each   ')
      print('time you create a project. It is up to you to:')
      print('(a) update code on cluster machines')
      print('(b) initiate new projects if any permanent')
      print('changes are made to the code')
      print('**********************************************')

      # step 0. delete the old local backup copy
      if (file.exists('./CARDAMOM_LOCAL/CARDAMOM_RECENT.zip') == TRUE) {
          system('rm CARDAMOM_LOCAL/CARDAMOM_RECENT.zip')
      }

      #step 1. zip cardamom
      print('Compressing CARDAMOM folder ...')
      system(paste("zip -r -q CARDAMOM.zip ",PROJECT$localpath,"*"))

      #step 2. move to PROJECT.exepath
      print(paste('Copying CARDAMOM.zip folder to ',PROJECT$exepath, ' ...',sep=""))
      system(paste('cp CARDAMOM.zip',PROJECT$exepath))

      print('Storing local copy: CARDAMOM_LOCAL/CARDAMOM_RECENT.zip ...')

      #step 4. store a local copy
      system(paste('mv CARDAMOM.zip CARDAMOM_LOCAL/CARDAMOM_RECENT.zip'))

      print('CARDAMOM successfully backed up!!!')
      print('**********************************************')

} # end function
