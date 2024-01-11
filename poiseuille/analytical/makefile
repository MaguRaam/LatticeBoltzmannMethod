CC = g++
CXXFLAGS = -std=c++17 -O3 -Wall
LIB =
SRC = analytical.cpp
EXE = a.out

$(EXE) : $(SRC)
	$(CC) $(CXXFLAGS) $(SRC) $(LIB)

clean:
	@rm -f $(EXE)

distclean:
	@rm -f *.vtk *.dat *.png $(EXE)
