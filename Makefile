CC = g++
CXXFLAGS = -std=c++17
LIB = -ltbb
SRC = test.cpp
EXE = a.out

$(EXE) : $(SRC)
	$(CC) $(CXXFLAGS) $(SRC) $(LIB)

clean:
	@rm -f $(EXE)

distclean:
	@rm -f *.vtk $(EXE)
