helloworld: main.o bastructs.o display.o
	g++ -std=c++11 main.o bastructs.o display.o -l SDL2 -l GL -l SDL2_net -o helloworld

display.o: src/Display.cpp src/Display.h src/Display.hpp
	g++ -std=c++11 -o display.o src/Display.cpp
	
bastructs.o: src/bastructs.cpp src/bastructs.hpp src/bastructs.h
	g++ -std=c++11 -o bastructs.o src/bastructs.cpp

main.o: src/main.cpp src/main.h
	g++ -std=c++11 -o main.o -c src/main.cpp 
