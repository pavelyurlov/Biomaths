CXX = clang++
CXXFLAGS = -std=c++14 -Wall -Wextra -pedantic -Wconversion -lm -lgsl -lgslcblas -g3 -march=native -O2 -isystem /usr/include/eigen3 # -fsanitize=memory

SRC = src/
build_dir = build/

vpath %.cpp $(SRC)
vpath %.h $(SRC)

all: main
		./build/a
		gnuplot plot_data

main: $(addprefix $(build_dir),main.o)
		$(CXX) $^ $(CXXFLAGS) -o $(build_dir)a

$(build_dir)%.o: %.cpp %.h
		$(CXX) $(CXXFLAGS) -c $< -o $@

$(build_dir)%.o: %.cpp
		$(CXX) $(CXXFLAGS) -c $< -o $@

clean:
		$(RM) -rf $(build_dir)*.o $(build_dir)a
