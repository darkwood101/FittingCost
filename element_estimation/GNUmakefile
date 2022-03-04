# Compiler flags, compile command
CXX := g++
CXXFLAGS := -W -Wextra -Wall -Wno-shadow -O2 -std=c++17
cxxcompile = @echo " " $(2) && $(CXX) $(CXXFLAGS) $(DEPCXXFLAGS) $(INC) $(1)

# Build target
all: model_uniform model_uniform_noise

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
model_uniform: $(OBJDIR)/model.o $(OBJDIR)/uniform.o
	$(call cxxcompile, -o $@ $^, LINK)
	$(call run, true, "Build successful!")

model_uniform_noise: $(OBJDIR)/model.o $(OBJDIR)/uniform_noise.o
	$(call cxxcompile, -o $@ $^, LINK)
	$(call run, true, "Build successful!")

# Remove dependencies, object files, and the executable
clean:
	$(call run, rm -rf $(OBJDIR) $(DEPSDIR) model_uniform, CLEAN)

.PHONY: all clean
