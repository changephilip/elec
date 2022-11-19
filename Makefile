CXX = g++
CXXFLAGS= -std=gnu++17 -shared -fPIC -fopenmp
external.so: external.cpp
	$(CXX)  $(CXXFLAGS) $^ -o $@

clean:
	rm *.so
