CC = g++

CFLAGS = -Wall \
         -Wno-reorder \
         -Wwrite-strings \
         -Wno-strict-aliasing \
         -Wno-unknown-pragmas

CDEBUG = -g

COPTIM =

CINC = -I$(HOME)/include \
       -I$(HOME)/software/include \
       -I.

CLIB = \
      -L$(HOME)/software/lib \
      -lACE \
      -lz \
      $(HOME)/software/lib/libniftiio.a \
      $(HOME)/software/lib/libznz.a \


OBJ_DIR = obj

.DELETE_ON_ERROR:

.SUFFIXES:
.SUFFIXES: .o .cpp

$(OBJ_DIR)/%.o: %.cpp %.h
	echo '[make: building $@]'
	$(CC) $(CFLAGS) $(CDEBUG) $(COPTIM) $(CINC) -o $@ -c $<

HDR_FILES = $(wildcard *.h)
ALL_SRC_FILES = $(wildcard ./*.cpp)
SRC_FILES = $(filter-out ./vsend_nii.cpp, $(ALL_SRC_FILES))
TMP_FILES = $(patsubst ./%,$(OBJ_DIR)/%,$(SRC_FILES))
OBJ_FILES = $(TMP_FILES:.cpp=.o)

default: prep $(OBJ_FILES) vsend_nii

prep:
	mkdir -p $(OBJ_DIR)

vsend_nii: $(OBJ_FILES) vsend_nii.cpp
	$(CC) $(CFLAGS) $(CDEBUG) $(COPTIM) $(CINC) -o $@ $@.cpp $(OBJ_FILES) $(CLIB)

clean:
	-rm $(OBJ_FILES)
	-rm vsend_nii
