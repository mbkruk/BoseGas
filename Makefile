LIB=lib/libbgcommon.a
ALL=$(LIB) bin/bgmc bin/bgev
DEBUG=

.PHONY: all
all: $(ALL)

build/%.cpp.o: %.cpp
	@mkdir -p $(@D)
	g++ -c $(DEBUG) -mavx -Ofast -std=c++11 $< -o $@

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
	g++ -s -flto $(DEBUG) $^ -o $@

bin/bgev: $(EV_OBJECTS) $(LIB)
	@mkdir -p $(@D)
	g++ -s -flto $(DEBUG) $^ -o $@

.PHONY: clean
clean:
	rm -rRf build
	rm -f $(ALL)
	rmdir lib bin
