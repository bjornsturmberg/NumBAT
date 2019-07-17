Using NumBAT with docker
========================

Start by obtaining the docker image and source code. Pull the docker image:

~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
docker pull morblockdock/numbat
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The docker image can be run simply with:

~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
docker run -it morblockdock/numbat bash
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

from which you can navigate to the example files and run them as described
elsewhere in the docs.

 

To have local access to the source files and simulation results the process is a
bit more involved.

If on linux, you must either run the following commands in sudo, or first grant yourself permission (otherwise you'll get the
`docker: Got permission denied while trying to connect to the Docker daemon socket` error).

~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
su - ${USER}
sudo usermod -aG docker ${USER}
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Now, create a directory `NumBAT-local`

~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
mkdir NumBAT-local
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

We now copy the source code into the NumBAT-local folder from the docker image. From the present working directory (ie without moving into NumBAT-local) execute the following:

~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
docker run --rm --entrypoint tar morblockdock/numbat cC /home/NumBAT/ . | tar xvC ./NumBAT-local
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

This command runs the docker image, copies the source code and then removes the
running docker container. 

Now we run the image as before, but in addition we mount the local directory,
allowing us to modify the files and obtain the results directly:

~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
docker run -v $(pwd)/NumBAT-local/:/home/NumBAT/ -it morblockdock/numbat bash
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

This runs the docker image, mounting the local directory to the path
/home/NumBAT/ in the docker dontainer, and runs the bash shell).

Changes made in the NumBAT source code now occur for the container also,
allowing for local development and easy access of results.

Finally, scripts can be run directly in a container by passing the path, for
example:

~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
docker run -v $(pwd)/NumBAT-local/:/home/NumBAT/ -it morblockdock/numbat python3 /home/NumBAT/tutorials/simo-tut_01-first_calc.py
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Runs the `simo-tut_01-first_calc.py` script directly from the local commandline
(without having to enter the docker container manually).

 

If you wish to build your own docker image, the local dockerfile
is provided as an example. From the higher directory you can build a docker
image using:

~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
docker build -t "numbat" NumBAT
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The dockerfile is a direct mapping of the setup script, in the syntax
currently preferred in the docker image build process
