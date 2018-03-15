CXX = clang++-3.8
CXXFLAGS = -std=c++14 -Wall -Wextra -pedantic -g -march=native  -I /usr/include/eigen3/Eigen -fsanitize=undefined

SRC = source/
build_dir = build/

vpath %.cpp $(SRC)
vpath %.h $(SRC)

all: main
	./build/a
	gnuplot graph
	gnuplot error

main: $(addprefix $(build_dir), main.o)
	$(CXX) $^ $(CXXFLAGS) -o $(build_dir)a

$(build_dir)%.o: %.cpp %.h
	$(CCX) $(CXXFLAGS) -c $< -o $@

$(build_dir)%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

clean:
	$(RM) -rf $(build_dir)*.o $(build_dir)a
