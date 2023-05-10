# PATH-SURVEYORS Docker Setup

Within this repository we have provided Dockerfiles and setup scripts for users to deploy the R shiny apps locally. Each of the folders represents a different R Shiny application from the PATH-SURVEYORS suite of tools that includes the Dockerfile, app files, and shell scripts users can follow to build, run, and save the Docker images of each application.

## Build

When the Dockerfile and app folder are in the users desired directory, the image can be built with the `step1_build_the_docker.sh` script, this will run the build function in the current directory.

```
docker build -t [app name] .
```

## Run

When the image build is completed, the app can then be run with the `step2_run_the_docker_image.sh` script

```
docker run -it --rm -p 3838:3838 [app name]
```

### Use

The app will then be accessible through a local instance. In your browser you can type http://127.0.0.1:3838 and the application will load.

## Save

The image that was generated can then be saved for easy access in the future. Users can save the image by running the `step4_saving_the_docker_image.sh` script.

```
docker save -o [app save name].tar [app name]
```

### Load app from save

In the future the app can then be loaded with the `step5_load_and_run_the_tar_image.sh` script.

```
docker load -i [app save name].tar
docker run -it --rm -p 3838:3838 [app name]:latest
```

#### Use

The application can then be accessed locally through typing the address http://127.0.0.1:3838 in your browser.
