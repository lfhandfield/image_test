helloworld: main.o
	g++ -std=c++11 main.o -l SDL2 -l GL -l SDL2_net -o helloworld

main.o: src/main.cpp src/main.h
	g++ -std=c++11 -o main.o -c src/main.cpp 
