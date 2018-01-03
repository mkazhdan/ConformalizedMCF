TARGET=ConformalizedMCF
TARGET_SRC=ConformalizedMCF/
TARGET_SOURCE=$(TARGET_SRC)ConformalizedMCF.cpp

VIEWER=ConformalizedMCFOrbifoldVisualization
VIEWER_SRC=ConformalizedMCFOrbifoldVisualization/
VIEWER_SOURCE=$(VIEWER_SRC)glew.c $(VIEWER_SRC)ConformalizedMCFOrbifoldVisualization.cpp

CPPFLAGS += -fpermissive -fopenmp -Wno-deprecated -std=c++11
CFLAGS += -fopenmp -Wno-deprecated
#LFLAGS += -lgomp /usr/lib64/libcholmod.so.1.7.3
LFLAGS += -lgomp

CPPFLAGS_DEBUG = -DDEBUG -g3
CFLAGS_DEBUG = -DDEBUG -g3
LFLAGS_DEBUG =

CPPFLAGS_RELEASE = -O3 -DRELEASE -funroll-loops -ffast-math
CLAGS_RELEASE = -O3 -DRELEASE -funroll-loops -ffast-math
LFLAGS_RELEASE = -O3

BIN = Bin/Linux/
INCLUDE =-I. -I/usr/include/ -ILibrary/

CC=gcc
CXX=g++
MD=mkdir

TARGET_OBJECTS=$(addprefix $(BIN), $(addsuffix .o, $(basename $(notdir $(TARGET_SOURCE)))))
VIEWER_OBJECTS=$(addprefix $(BIN), $(addsuffix .o, $(basename $(notdir $(VIEWER_SOURCE)))))

all: CPPFLAGS += $(CPPFLAGS_RELEASE)
all: CFLAGS += $(CFLAGS_RELEASE)
all: LFLAGS += $(LFLAGS_RELEASE)
all: $(BIN)
all: $(BIN)$(TARGET)
all: $(BIN)$(VIEWER)

debug: CPPFLAGS += $(CPPFLAGS_DEBUG)
debug: CFLAGS += $(CFLAGS_DEBUG)
debug: LFLAGS += $(LFLAGS_DEBUG)
debug: $(BIN)
debug: $(BIN)$(TARGET)
debug: $(BIN)$(VIEWER)

clean:
	rm -f $(BIN)$(TARGET)
	rm -f $(BIN)$(VIEWER)
	rm -f $(OBJECTS)

$(BIN):
	$(MD) -p $(BIN)
	
$(BIN)$(TARGET): $(TARGET_OBJECTS)
	$(CXX) -o $@ $(TARGET_OBJECTS) $(LFLAGS)

$(BIN)$(VIEWER): $(VIEWER_OBJECTS)
	$(CXX) -o $@ $(VIEWER_OBJECTS) $(LFLAGS)
	
$(BIN)%.o: $(TARGET_SRC)%.c
	$(CC) -c -o $@ $(CFLAGS) $(INCLUDE) $<

$(BIN)%.o: $(VIEWER_SRC)%.c
	$(CC) -c -o $@ $(CFLAGS) $(INCLUDE) $<

$(BIN)%.o: $(TARGET_SRC)%.cpp
	$(CXX) -c -o $@ $(CPPFLAGS) $(INCLUDE) $<

$(BIN)%.o: $(VIEWER_SRC)%.cpp
	$(CXX) -c -o $@ $(CPPFLAGS) $(INCLUDE) $<
