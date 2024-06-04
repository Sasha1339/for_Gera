# Define variables
CC = gcc  # Compiler
CFLAGS = $(shell pkg-config gtk+-2.0 --cflags) -I/usr/local/include  # Compilation flags using pkg-config and audio include directory
LDFLAGS = $(shell pkg-config gtk+-2.0 --libs) -lmpg123 -lao -lsndfile  # Linking flags using pkg-config and audio libraries

TARGET = main  # Output executable name

# Source files
SRCS = main.c buttons.c audio.c staff.c z_fft.c z_math.c # Source files including audio.c
OBJS = $(SRCS:.c=.o)  # Corresponding object files

all: clean $(TARGET)

# Default target: compile and link the program
all: $(TARGET)

# Rule to build the executable
$(TARGET): $(OBJS)
	$(CC) $(OBJS) -o $(TARGET) $(LDFLAGS)

# Rule to compile each source file into object files
%.o: %.c
	$(CC) $(CFLAGS) -c $< -o $@

# Clean up build artifacts
.PHONY: clean
clean:
	rm -f $(TARGET) $(OBJS) audio