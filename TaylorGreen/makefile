CC = gcc
CFLAGS = -O3 -Wall -pedantic
EXE = a.out
SRC = cpu.c

$(EXE) : $(SRC)
	$(CC) $(CFLAGS) $^ -o $@ -lm

clean:
	@rm -f $(EXE)

distclean:
	@rm -f $(EXE) *.vtk