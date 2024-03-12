CC = g++
 
# compiler flags:
#  -g     - this flag adds debugging information to the executable file
#  -Wall  - this flag is used to turn on most compiler warnings
CFLAGS  = -g -Wall

# The build target 
TARGET = project3
SRC = $(TARGET).cpp
EXEC = $(TARGET).exe

all: $(EXEC)

$(EXEC): $(SRC)
		$(CC) $(CFLAGS) -o $(EXEC) $(SRC)

clean:
		$(RM) $(EXEC)