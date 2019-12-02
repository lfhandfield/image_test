helloworld: main.o bastructs.o primitive.o display.o glad.o
	g++ -std=c++11 main.o bastructs.o primitive.o display.o glad.o -l SDL2 -l GL -l SDL2_net -lpthread -o helloworld

display.o: src/Display.cpp src/Display.h src/Display.hpp
	g++ -std=c++11 -o display.o src/Display.cpp
	
bastructs.o: src/bastructs.cpp src/bastructs.hpp src/bastructs.h
	g++ -std=c++11 -o bastructs.o src/bastructs.cpp
	
primitive.o: src/primitive.cpp src/primitive.hpp src/primitive.h
	g++ -std=c++11 -o primitive.o src/primitive.cpp

glad.o: src/glad.c src/glad.h src/khrplatform.h
	g++ -std=c++11 -o ./glad.o -c src/glad.c

main.o: src/main.cpp src/main.h
	g++ -std=c++11 -o main.o -c src/main.cpp 
