#This is a sample Image 
FROM ubuntu 

RUN apt-get update 
RUN apt-get install –y libsdl2-2.0-0 libsdl2-net-2.0-0
CMD [“echo”,”Image created”] 
