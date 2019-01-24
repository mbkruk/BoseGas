LIB=lib/libbgcommon.a
ALL=$(LIB) bin/bgmc bin/bgev
DEBUG=

GPP_VERSION=$(shell g++ -dumpversion)
STD=$(shell if [ "$(GPP_VERSION)" \> "4" ]; then echo "-std=c++11"; else echo "-std=c++0x"; fi)
CFLAGS=-mavx -O3
LDFLAGS=$(shell [ "$(GPP_VERSION)" \> "4" ] && echo "-flto")

.PHONY: all
all: $(ALL)

build/%.cpp.o: %.cpp
	@mkdir -p $(@D)
	g++ -c $(DEBUG) $(STD) $(CFLAGS) $< -o $@

LIB_SOURCES=$(shell find -wholename "./source/BGCommon/*.cpp")
LIB_OBJECTS=$(LIB_SOURCES:%=build/%.o)

MC_SOURCES=$(shell find -wholename "./source/BGMC/*.cpp")
MC_OBJECTS=$(MC_SOURCES:%=build/%.o)

EV_SOURCES=$(shell find -wholename "./source/BGEV/*.cpp")
EV_OBJECTS=$(EV_SOURCES:%=build/%.o)

$(LIB): $(LIB_OBJECTS)
	@mkdir -p $(@D)
	ar rcs $@ $^

bin/bgmc: $(MC_OBJECTS) $(LIB)
	@mkdir -p $(@D)
	g++ -s $(LDFLAGS) $(DEBUG) $^ -o $@

bin/bgev: $(EV_OBJECTS) $(LIB)
	@mkdir -p $(@D)
	g++ -s $(LDFLAGS) $(DEBUG) $^ -o $@

.PHONY: clean
clean:
	rm -rRf build
	rm -f $(ALL)
	rmdir lib bin
