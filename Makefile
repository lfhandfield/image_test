helloworld: main.o core.o cstructs.o bastructs.o primitive.o display.o glad.o
	g++ -std=c++11 main.o core.o cstructs.o bastructs.o primitive.o display.o glad.o -l SDL2 -l GL -l  SDL2_net -ldl -lpthread -o helloworld

display.o: src/Display.cpp src/Display.h src/Display.hpp
	g++ -std=c++11 -o display.o -c src/Display.cpp
	
bastructs.o: src/bastructs.cpp src/bastructs.hpp src/bastructs.h
	g++ -std=c++11 -o bastructs.o -c src/bastructs.cpp

cstructs.o: src/cstructs.cpp src/cstructs.h
	g++ -std=c++11 -o cstructs.o -c src/cstructs.cpp

core.o: src/core.cpp src/core.h	src/core.hpp
	g++ ${COMPILER_FLAGS} -o core.o -c src/core.cpp
	
primitive.o: src/primitive.cpp src/primitive.hpp src/primitive.h
	g++ -std=c++11 -o primitive.o -c src/primitive.cpp

glad.o: src/glad.c src/glad.h src/khrplatform.h
	g++ -std=c++11 -o ./glad.o -c src/glad.c

main.o: src/main.cpp src/main.h
	g++ -std=c++11 -o main.o -c src/main.cpp 

