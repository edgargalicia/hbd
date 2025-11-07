CXX := g++
CXXFLAGS := -std=c++14 -Wall -Wextra -O2
LDFLAGS :=
TARGET := build/hbd

SRC_DIRS := . Math
SRCS := $(foreach dir,$(SRC_DIRS),$(wildcard $(dir)/*.cpp))
OBJS := $(patsubst %.cpp,build/%.o,$(SRCS))
DEPS := $(OBJS:.o=.d)

CTAGS := ctags
CTAGS_FLAGS := -R --exclude='build*' --exclude='.git'
CTAGS_FILE := tags

.PHONY: all clean dirs tags

all: $(TARGET) tags

$(TARGET): $(OBJS)
	@echo "Linking $@"
	$(CXX) $(OBJS) -o $@ $(LDFLAGS)

build/%.o: %.cpp | dirs
	@echo "Compiling $<"
	@mkdir -p $(dir $@)
	$(CXX) $(CXXFLAGS) -MMD -MP -c $< -o $@

dirs:
	@mkdir -p $(addprefix build/,$(SRC_DIRS))

tags:
	@echo "Updating ctags..."
	@$(CTAGS) $(CTAGS_FLAGS) -f $(CTAGS_FILE) $(SRC_DIRS)

-include $(DEPS)

clean:
	rm -rf build
