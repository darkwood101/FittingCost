# Compiler flags, compile command
CXX := g++
CXXFLAGS := -W -Wextra -Wall -Wno-shadow -O2 -std=c++17
cxxcompile = @echo " " $(2) && $(CXX) $(CXXFLAGS) $(DEPCXXFLAGS) $(INC) $(1)

# Build target
all: main generate_data

# Generic run command
run = @$(if $(2), echo " " $(2) $(3) &&,) $(1) $(3)

# Dependencies
DEPSDIR := .deps
DEPFILES := $(wildcard $(DEPSDIR)/*.d)
include $(DEPFILES)
DEPCXXFLAGS = -MD -MF $(DEPSDIR)/$(@F).d -MP

# Header files
INC := -Iinclude

# Source files
SRCFILES := $(wildcard *.cc)

# Object files
OBJDIR := obj
OBJFILES := $(patsubst %.cc, $(OBJDIR)/%.o, $(SRCFILES))

# Create directories if they don't exist
$(OBJFILES): | $(OBJDIR) $(DEPSDIR)

$(OBJDIR):
	$(call run, mkdir -p $@, CREATE $@/)

$(DEPSDIR):
	$(call run, mkdir -p $@, CREATE $@/)

# How to make object files
$(OBJDIR)/%.o: %.cc
	$(call cxxcompile, -o $@ -c $<, COMPILE $<)

# How to make the executable
main: $(OBJDIR)/main.o
	$(call cxxcompile, -o $@ $^, LINK $^)
	$(call run, true, OUTPUT $@)

generate_data: $(OBJDIR)/generate_data.o
	$(call cxxcompile, -o $@ $^, LINK $^)
	$(call run, true, OUTPUT $@)

# Remove dependencies, object files, and the executable
clean:
	$(call run, rm -rf $(OBJDIR) $(DEPSDIR) \
					   main \
					   generate_data \
					   *.out, CLEAN)

.PHONY: all clean

