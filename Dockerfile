####this is the one to do read and write 
#    used in combination with compose.yml file.
#      (we don't have to do in combination, but we could)
#	https://www.robertmylesmcdonnell.com/content/posts/docker/
# but we do convert not to use yaml (docker compose)
#	but only use docker to add named volumes
FROM rocker/rstudio:3.6.3
#create a folder in container and entering it.
#something it is not working well at the root directory.
WORKDIR /home/rstudio/IgSeqR3
#note about directory, we use the bind volume to the container volum main
#we do this binding with the docker run command switches -w /main -v "$(pwd):/main" 
# At running time (not compiling/building time)
#  docker run \
#     -w /main -v "$(pwd):/main" \
#     -it "readwrite"
# #building:
# docker build -t readwrite .
#RUN touch little_script.R
#Run this in order to plot x11 figures

RUN sudo apt-get update && apt-get install -y --no-install-recommends libtcl8.6 tcl8.6-dev libtk8.6 tk8.6-dev

RUN sudo apt-get update && apt-get install -y --no-install-recommends libxt6 \
    && apt-get install -y libz-dev cmake sudo\
    && apt-get install -y libpng-dev libfontconfig1-dev \
    && apt-get install -y libharfbuzz-dev libfribidi-dev \
    && apt-get install -y libfreetype6-dev libtiff5-dev libjpeg-dev \
    && apt-get install -y libxml2-dev libbz2-dev libssl-dev \
    && apt-get install -y --no-install-recommends libpng-dev libproj-dev pandoc python3 libcurl4-openssl-dev libglu1-mesa-dev libgl1-mesa-dev zlib1g-dev libicu-dev imagemagick libmagick++-dev gsfonts 
    
    

RUN R -e "install.packages('renv', repos = c(CRAN = 'https://cloud.r-project.org'))"

# approach two
# to set up the 
COPY renv.lock renv.lock
#COPY renv.lock /home/renv.lock
RUN  mkdir -p renv
COPY .Rprofile .Rprofile
COPY renv/activate.R renv/activate.R
COPY renv/settings.dcf renv/settings.dcf

#COPY renv/settings.dcf renv/settings.dcf
# Note for above, we don't have renv/settings.json somehow (?),
# it is not necessary!!! we could just need .Rprofile and activate.R and renv.lock
#

#change owner ship to rstudio, RECURSIVELY with -R
# and then restor the renv library as the user "rstudio"
RUN chown -R rstudio . \
    && sudo -u rstudio R -e "renv::restore()"
#RUN chown rstudio /home/renv.lock

RUN sudo -u rstudio echo "setwd(\"/home/rstudio/IgSeqR3/\")" > /home/rstudio/.Rprofile

RUN sudo -u rstudio echo "renv::load()" >> /home/rstudio/.Rprofile 
#    && echo ".rs.restartR()" >>/home/rstudio/.Rprofile 

COPY README.md README.md
COPY Readme.txt Readme.txt
COPY environment.yml environment.yml
COPY R_code/ ./R_code/

#this is for rstudio. so we need to change the ownership
#RUN chown rstudio /app/KerasR_test

#add current directory to show in the print.
#RUN echo "print('hi, every-body!$(pwd)')" >> little_script.R
#RUN echo "created an R file. now running.....!!\n"
#CMD ["Rscript", "readwrite.R"]
