# GeneidTRAINerDocker
Docker container containing the perl script that we use to train the ab initio program geneid


To start the docker container GeneidTRAINerDocker type the following:

docker build -t GeneidTRAINerDocker .

In order to run it:

docker run -v $(pwd)/input:/input -v $(pwd)/output:/output GeneidTRAINerDocker
