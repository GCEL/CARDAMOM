
# Docker

Docker allows for software to be packaged together with the operating system environment necessary for it to run.
From [Wikipedia](https://en.wikipedia.org/wiki/Docker_(software))

> Docker is a tool that is used to automate the deployment of applications in lightweight containers so that applications can work efficiently in different environments in isolation.

## Running the CARDAMOM Docker container

You will need a laptop or desktop computer with Windows, MacOS or Linux that you have administrative privelages on.

Follow these steps:

  1. Create an account at [docker hub](http://hub.docker.com).
  2. Download and install Docker Desktop for your operating system.
     More information available at [the Docker documentation](https://docs.docker.com/get-docker/).
  3. Open Docker Desktop and login with your Docker Hub credentials.
  4. In the search bar, look for cardamomtest/cardamom. Click 'Pull' next to the result.
  5. After the image is downloaded, select 'Images' on left panel and then 'Local'
     in tab at the top of the centre panel.
  6. Play under actions
  1. Click on optional settings
  1. Under volumes, Host path browse to shared directory on the host computer and set 'Container path' to '/home/cardamom/host'
  1. Click Run
  1. In the center panel, go the 'Terminal' tab and follow the steps in [README_GIT_GUIDANCE.md]().
     
