
CC      := cc
CFLAGS  := -std=c11  -g -O0 -Wall -Wextra -I./include
BUILD   := build
SRC     := src
BIN     := $(BUILD)/app

OBJS := $(BUILD)/cyclotomic_cosets.o $(BUILD)/mymath.o $(BUILD)/temporal.o

$(BIN): $(OBJS)
	$(CC) $(OBJS) -o $@ -lm

$(BUILD)/%.o: $(SRC)/%.c | $(BUILD)
	$(CC) $(CFLAGS) -c $< -o $@

$(BUILD):
	mkdir -p $(BUILD)

.PHONY: clean
clean:
	rm -rf $(BUILD)
