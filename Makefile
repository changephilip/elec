CXX = g++
CXXFLAGS= -std=gnu++17 -shared -fPIC
dielec.so: dielec.cpp
	$(CXX)  $(CXXFLAGS) $^ -o $@
