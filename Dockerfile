#This is a sample Image 
FROM ubuntu 

RUN apt-get update 
RUN apt-get -y install build-essential
RUN apt-get -y install libglu1-mesa-dev mesa-common-dev
RUN apt-get -y install libsdl2-2.0-0 libsdl2-net-2.0-0
RUN apt-get -y install libsdl2-dev libsdl2-net-dev
WORKDIR /opt/display
COPY src /src
COPY Makefile /Makefile
RUN ls
RUN make
ENV PATH="/opt/display/:${PATH}"
RUN echo Image created
CMD ["helloworld"]
