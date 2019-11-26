Helloworld: main.o
  g++ -std=c++11 main.o -l SDL2 -l GL -l SDL2_net

Helloworld: src/main.cpp src/main.h
  g++ -std=c++11 src/main.cpp
  
  
