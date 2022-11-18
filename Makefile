CXX = g++
CXXFLAGS= -std=gnu++17 -shared -fPIC -fopenmp
dielec.so: dielec.cpp
	$(CXX)  $(CXXFLAGS) $^ -o $@

clean:
	rm *.so
