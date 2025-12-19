# Makefile for building mesher, ttimes, and gravity binaries

# Compiler and linker
CC := gcc
MPICC := mpicc

# Compiler flags
FLAGS := \
	-Wall \
	-Wextra \
	-std=c99 \
	-O3

# Directories
SRC := src
SHARED := $(SRC)/shared
SETUP := setup
OBJ := obj
BIN := bin

# Includes and libraries
INC := -I$(SETUP) -I$(SRC)/headers
LIB := -lm

# Object files
OBJECTS := \
	io.o \
	exmath.o \
	coordinates.o \
	tesselation.o \
	progress.o

# Binary targets
DEFAULT := \
	mesher \
	ttimes \
	gravity

# Complete paths for binaries and objects
DFT := $(patsubst %, $(BIN)/%, $(DEFAULT))
OBJS := $(patsubst %.o, $(OBJ)/%.o, $(OBJECTS))

# Command for cleaning
RM := rm -rf

#
# Compilation and linking
#
all: objDirectory binDirectory $(DFT)
	@ echo 'Finished building binaries!'

$(BIN)/mesher: $(OBJS) $(OBJ)/mesher.o
	@ echo 'Building binary using $(CC) linker: $@'
	$(CC) $(FLAGS) $(INC) $^ -o $@ $(LIB)
	@ echo 'Finished building binary: $@'
	@ echo ' '

$(BIN)/ttimes: $(OBJS) $(OBJ)/ttimes.o
	@ echo 'Building binary using $(CC) linker: $@'
	$(CC) $(FLAGS) $(INC) $^ -o $@ $(LIB)
	@ echo 'Finished building binary: $@'
	@ echo ' '

$(BIN)/gravity: $(OBJS) $(OBJ)/gravity.o
	@ echo 'Building binary using $(MPICC) linker: $@'
	$(MPICC) $(FLAGS) $(INC) $^ -o $@ $(LIB)
	@ echo 'Finished building binary: $@'
	@ echo ' '

$(OBJ)/io.o: $(SHARED)/io.c $(SETUP)
	@ echo 'Building target using $(CC) compiler: $@'
	$(CC) $(FLAGS) $(INC) -c $< -o $@
	@ echo ' '

$(OBJ)/exmath.o: $(SHARED)/exmath.c $(SETUP)
	@ echo 'Building target using $(CC) compiler: $@'
	$(CC) $(FLAGS) $(INC) -c $< -o $@
	@ echo ' '

$(OBJ)/coordinates.o: $(SHARED)/coordinates.c $(SETUP)
	@ echo 'Building target using $(CC) compiler: $@'
	$(CC) $(FLAGS) $(INC) -c $< -o $@
	@ echo ' '

$(OBJ)/tesselation.o: $(SHARED)/tesselation.c $(SETUP)
	@ echo 'Building target using $(CC) compiler: $@'
	$(CC) $(FLAGS) $(INC) -c $< -o $@
	@ echo ' '

$(OBJ)/progress.o: $(SHARED)/progress.c $(SETUP)
	@ echo 'Building target using $(CC) compiler: $@'
	$(CC) $(FLAGS) $(INC) -c $< -o $@
	@ echo ' '

$(OBJ)/mesher.o: $(SRC)/mesher.c $(SETUP)
	@ echo 'Building target using $(CC) compiler: $@'
	$(CC) $(FLAGS) $(INC) -c $< -o $@
	@ echo ' '

$(OBJ)/ttimes.o: $(SRC)/ttimes.c $(SETUP)
	@ echo 'Building target using $(CC) compiler: $@'
	$(CC) $(FLAGS) $(INC) -c $< -o $@
	@ echo ' '

$(OBJ)/gravity.o: $(SRC)/gravity.c $(SETUP)
	@ echo 'Building target using $(MPICC) compiler: $@'
	$(MPICC) $(FLAGS) $(INC) -c $< -o $@
	@ echo ' '

objDirectory:
	@ mkdir -p $(OBJ)

binDirectory:
	@ mkdir -p $(BIN)

clean:
	$(RM) $(OBJ)/ $(BIN)/

.PHONY: all clean
