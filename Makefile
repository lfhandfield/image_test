helloworld: main.o display.o
	g++ -std=c++11 main.o display.o -l SDL2 -l GL -l SDL2_net -o helloworld

display.o: src/display.cpp src/display.h src/display.hpp
	g++ -std=c++11 -o display.o src/display.cpp

main.o: src/main.cpp src/main.h
	g++ -std=c++11 -o main.o -c src/main.cpp 
