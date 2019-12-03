FROM ubuntu 

ENV DEBIAN_FRONTEND=noninteractive
RUN apt-get update 
RUN apt-get -y install build-essential
RUN apt-get -y install libglu1-mesa-dev mesa-common-dev libsdl2-2.0-0 libsdl2-dev glew-utils libglew-dev libsdl2-net-2.0-0 libsdl2-net-dev 
RUN apt-get -y install libc6-dbg gdb valgrind
WORKDIR /opt/display
COPY src /opt/display/src
COPY Images /opt/display/Images
COPY Makefile /opt/display/Makefile
RUN ls
RUN make
ENV PATH="/opt/display/:${PATH}"
RUN echo Image created
ENV PATH="/opt/display/:${PATH}"
RUN ./helloworld
CMD ["helloworld"]
