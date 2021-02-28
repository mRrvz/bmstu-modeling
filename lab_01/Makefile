CC=ghc
MAIN=Main.hs
WORKDIR=src
EXENAME=app.exe
FLAGS=-package boxes -O3 -threaded

build:
	cd $(WORKDIR) && $(CC) $(FLAGS) -o ../$(EXENAME) --make $(MAIN)

run: build
	./$(EXENAME) +RTS -N4

clean:
	rm $(WORKDIR)/*.o $(WORKDIR)/*.hi *.exe *.out