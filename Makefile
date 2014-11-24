all:
	g++ -O3 -lopencv_core -lopencv_highgui feature.cpp deepmaxent.cpp densitySample.cpp densityGrid.cpp main.cpp data.cpp plot.cpp

debug:
	g++ -O0 -g -lopencv_core -lopencv_highgui feature.cpp deepmaxent.cpp densitySample.cpp densityGrid.cpp main.cpp data.cpp plot.cpp
