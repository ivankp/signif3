CXX := g++
STD := -std=c++14
CPPFLAGS := $(STD) -Isrc
CXXFLAGS := $(STD) -Wall -O3 -Isrc -fmax-errors=3
# CXXFLAGS := $(STD) -Wall -g -Iinclude -fmax-errors=3
LDFLAGS := $(STD) -O3
LDLIBS :=

SRC := src
BIN := bin
BLD := .build
EXT := .cc

.PHONY: all clean

ifeq (0, $(words $(findstring $(MAKECMDGOALS), clean)))

LDFLAGS += $(shell sed -r 's/([^:]+)(:|$$)/ -L\1/g' <<< "$$LIBRARY_PATH")

ROOT_CXXFLAGS := $(shell root-config --cflags)
ROOT_LDFLAGS  := $(shell root-config --ldflags)
ROOT_LDLIBS   := $(shell root-config --libs)

# RPATH
rpath_script := ldd $(shell root-config --libdir)/libTreePlayer.so \
  | sed -nr 's|.*=> (.+)/.+\.so[.0-9]* \(.*|\1|p' \
  | sort -u \
  | sed -nr '/^(\/usr)?\/lib/!s/^/-Wl,-rpath=/p'
ROOT_LDLIBS += $(shell $(rpath_script))

CXXFLAGS += $(ROOT_CXXFLAGS)
LDFLAGS += $(ROOT_LDFLAGS)
LDLIBS += $(ROOT_LDLIBS) -lTreePlayer

SRCS := $(shell find $(SRC) -type f -name '*$(EXT)')
DEPS := $(patsubst $(SRC)/%$(EXT),$(BLD)/%.d,$(SRCS))

GREP_EXES := grep -rl '^ *int \+main *(' $(SRC) --include='*$(EXT)'
EXES := $(patsubst $(SRC)%$(EXT),$(BIN)%,$(shell $(GREP_EXES)))

all: $(EXES)

$(BIN)/test $(BIN)/signif $(BIN)/mig \
$(BIN)/superfine $(BIN)/optimize \
$(BIN)/simple_signif: $(BLD)/re_axes.o

-include $(DEPS)

.SECONDEXPANSION:

$(DEPS): $(BLD)/%.d: $(SRC)/%$(EXT) | $(BLD)/$$(dir %)
	$(CXX) $(CPPFLAGS) -MM -MT '$(@:.d=.o)' $< -MF $@

$(BLD)/%.o: | $(BLD)
	$(CXX) $(CXXFLAGS) $(C_$*) -c $(filter %$(EXT),$^) -o $@

$(BIN)/%: $(BLD)/%.o | $(BIN)
	$(CXX) $(LDFLAGS) $(filter %.o,$^) -o $@ $(LDLIBS) $(L_$*)

$(BIN):
	mkdir -p $@

$(BLD)/%/:
	mkdir -p $@

endif

clean:
	@rm -rfv $(BLD) $(BIN)
