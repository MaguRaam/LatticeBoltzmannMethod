CC = g++
CXXFLAGS = -std=c++17 -O3 -Wall -Wextra
LIB = -ltbb
SRC = poiseuille.cpp
EXE = a.out

$(EXE) : $(SRC)
	$(CC) $(CXXFLAGS) $(SRC) $(LIB)

clean:
	@rm -f $(EXE)

distclean:
	@rm -f *.vtk $(EXE)