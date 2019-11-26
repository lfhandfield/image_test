#This is a sample Image 
FROM ubuntu 

RUN apt-get update 
RUN apt-get -y install build-essential libsdl2-2.0-0 libsdl2-net-2.0-0 libglu1-mesa-dev mesa-common-dev r-base
RUN ls
RUN make
RUN echo Image created
CMD ["helloworld"]
