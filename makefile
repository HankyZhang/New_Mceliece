CC = gcc
CFLAGS = -Wall -Wextra -O2 -std=c99 -g
TARGET = mceliece
SOURCES = main.c mceliece_utils.c mceliece_gf.c mceliece_shake.c \
          mceliece_random.c mceliece_matrix_ops.c mceliece_berlekamp.c mceliece_kem.c mceliece_test.c

OBJECTS = $(SOURCES:.c=.o)
HEADERS = mceliece_types.h

all: $(TARGET)

$(TARGET): $(OBJECTS)
	$(CC) $(OBJECTS) -o $(TARGET) -lm

%.o: %.c $(HEADERS)
	$(CC) $(CFLAGS) -c $< -o $@

clean:
	rm -f $(OBJECTS) $(TARGET)

test: $(TARGET)
	./$(TARGET) test

fulltest: $(TARGET)
	./$(TARGET) fulltest

demo: $(TARGET)
	./$(TARGET) demo

keygen: $(TARGET)
	./$(TARGET) keygen

bench: $(TARGET)
	./$(TARGET) bench

debug: CFLAGS += -DDEBUG -O0
debug: $(TARGET)

install: $(TARGET)
	cp $(TARGET) /usr/local/bin/

.PHONY: all clean test demo keygen bench debug install