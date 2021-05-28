
###
## Function automates accessing the UoE HPC service Eddie
## and issues commands to Eddie regarding CARDAMOM
###

# This function is by T. L Smallman (t.l.smallman@ed.ac.uk, UoE).

ecdf_execute<- function(command,cluster_address) {

    # create fifo  connection for ssh
    if (file.exists("sshOut") == TRUE ){ system('rm sshOut')}
    system('mkfifo sshOut')
    #Then, connect to the pipe and the fifo from within R:
    # The &> redirects both stdout and stderr to the fifo
    sshIn <- pipe( paste("ssh -A ",username,"@",cluster_address," &> sshOut",sep=""), open = 'w+')
    sshOut <- fifo( 'sshOut', 'r+' )
    # write commands to list
    for (i in seq(1, length(command))) {writeLines( command[i], sshIn )}
    # 'flush' pass arguements and then read remote systems return
    flush( sshIn )
    readLines( sshOut )
    # clean up
    close(sshIn)
    system('rm sshOut')

} # end function ecdf_execute
