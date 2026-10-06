# Build and run targets for the joint zero-difference experiment.
# Use from a POSIX shell (Linux, macOS, or Git Bash / MSYS2 / WSL on Windows).
#
#   make              build ./joint_exp
#   make test         check the cipher implementation (seconds)
#   make table3       Table 3 from the stored paper data (seconds)
#   make quick        4 runs x 2^24 quartets (< 1 minute)
#   make medium       16 runs x 2^28 quartets, in parallel (minutes)
#   make full         100 runs x 2^30 quartets, in parallel (about 1 h on 16 cores)
#
# With the default SEED, the outputs of quick/medium/full are also compared
# with the stored reference outputs in Results/ (all columns except elapsed_sec).
# Variables: JOBS (parallel processes), SEED, PYTHON, CC, CFLAGS, CPPFLAGS.
# Other round numbers: make -B CPPFLAGS=-DNUM_ROUNDS=5 (with 5 rounds the mu-check
# fails); run `make -B` afterwards to rebuild the default 4-round ./joint_exp.
#   make clean        remove the binaries (also joint_exp5) and out/quick.csv
#   make distclean    also remove out/

CC      = gcc
CFLAGS  = -O2 -std=c99 -Wall -Wextra
CPPFLAGS =
PYTHON ?= $(shell python3 -c "import sys" >/dev/null 2>&1 && echo python3 || echo python)
JOBS   ?= $(shell getconf _NPROCESSORS_ONLN 2>/dev/null || echo 4)
SEED   ?= 20261002
EXE     = $(if $(filter Windows_NT,$(OS)),.exe,)
BIN     = joint_exp$(EXE)
SRC     = Code/joint_property_experiment.c

# reference options for the default seed, nothing otherwise
ref  = $(if $(filter 20261002,$(SEED)),--reference Results/$(1),)
pref = $(if $(filter 20261002,$(SEED)),-r Results/$(1),)
# fail early if Python is missing
pycheck = @$(PYTHON) -c "import sys; sys.exit(sys.version_info < (3, 6))" \
	|| { echo "Python >= 3.6 not found; set PYTHON=..."; exit 1; }

.PHONY: all test table3 quick medium full clean distclean

all: $(BIN)

$(BIN): $(SRC)
	$(CC) $(CPPFLAGS) $(CFLAGS) -o $@ $(SRC) -lm

test: $(BIN)
	./$(BIN) --self-test

table3:
	$(PYTHON) scripts/summarize.py Results/results.csv

quick: $(BIN)
	$(pycheck)
	@mkdir -p out
	./$(BIN) --quick --seed $(SEED) --force -o out/quick.csv
	$(PYTHON) scripts/summarize.py out/quick.csv $(call ref,quick_seed20261002.csv)

medium: $(BIN)
	PYTHON=$(PYTHON) sh scripts/run_parallel.sh -j $(JOBS) -k 16 -s $(SEED) -o out/medium \
	  $(call pref,medium_seed20261002.csv) -- --log2-trials 28

full: $(BIN)
	PYTHON=$(PYTHON) sh scripts/run_parallel.sh -j $(JOBS) -k 100 -s $(SEED) -o out/full \
	  $(call pref,results_v2.csv) -- --mu-trials 0

clean:
	rm -f joint_exp joint_exp.exe joint_exp5 joint_exp5.exe out/quick.csv

distclean: clean
	rm -rf out
