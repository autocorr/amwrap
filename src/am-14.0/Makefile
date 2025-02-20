# Makefile for am-14.0
# S. Paine rev. 2024 August 6


# The default target is a help message:

help:
	@echo ""
	@echo "Building (with gcc) and installing am:"
	@echo "  make am         # build parallel am (OpenMP, recommended)"
	@echo "  make serial     # build single-threaded am"
	@echo "  make install    # install to $(INSTALL_DIR) (as root, or via sudo)"
	@echo "  make clean      # clean up build files (required before changing"
	@echo "                  # the build target or compiler)"
	@echo "  make devel      # list additional build targets for development"
	@echo ""
	@echo "Building with other compilers:"
	@echo "  make nvc-omp    # Nvidia nvc C compiler"
	@echo "  make nvc-serial"
	@echo "  make icx-omp    # Intel icx oneAPI C compiler"
	@echo "  make icx-serial"
	@echo ""
	@echo "Use EXTRA_CFLAGS, EXTRA_LDFLAGS, and EXTRA_LDLIBS to add compiler"
	@echo "and linker flags.  Some examples:"
	@echo ""
	@echo "  Generate optimization and vectorization reports:"
	@echo "    make gcc-omp 'EXTRA_CFLAGS = -fopt-info'"
	@echo "    make nvc-omp 'EXTRA_CFLAGS = -Minfo -Mneginfo'"
	@echo "    make icx-omp 'EXTRA_CFLAGS = -qopt-report'"
	@echo ""
	@echo "  Build with generic x86-64 support, useful on heterogeneous clusters:"
	@echo "    make nvc-omp 'EXTRA_CFLAGS = -tp px'"
	@echo ""
	@echo "  Tune for 16K 4-way L1 cache, overriding default (32K 8-way):"
	@echo "    make gcc-omp 'EXTRA_CFLAGS = -DL1_CACHE_BYTES=0x4000 -DL1_CACHE_WAYS=4'"
	@echo ""


# make devel lists targets for development:

devel:
	@echo ""
	@echo "Default gcc build targets:"
	@echo "  make gcc-omp           # gcc, OpenMP multi-threaded"
	@echo "  make gcc-serial        # gcc, single-threaded"
	@echo ""
	@echo "Build without optimization:"
	@echo "  make gcc-no-opt        # gcc, single-threaded with no optimization"
	@echo "  make gcc-no-opt-debug  # gcc, single-threaded for debugging"
	@echo ""
	@echo "Compile with instrumentation for the GNU gcov coverage testing tool:"
	@echo "  make gcc-no-opt 'EXTRA_CFLAGS = --coverage' 'EXTRA_LDFLAGS = --coverage'"
	@echo "  make clean-gcov        # clean up gcov files in source directory"
	@echo ""
	@echo "Make a set of Windows native (crlf) source files, in a new directory"
	@echo "one level up from this one (requires sed):"
	@echo "  make crlf-source"
	@echo ""


clean:
	rm -f $(OBJS) version.o am

clean-gcov:
	rm -f *.gcov *.gcda *.gcno

# INSTALL_DIR is the installation directory used by the 'install'
# target.  /usr/local/bin is the appropriate directory on modern
# Unix-like systems that follow the Filesystem Hierarchy Standard.
#
# On macOS, the directory /usr/local/bin is listed in the file
# /etc/paths which is read by the path_helper utility that constructs
# the shell PATH environment variable.  However, /usr/local/bin may
# not actually exist on a fresh macOS installation. The 'install'
# target will create it if necessary.

INSTALL_DIR = /usr/local/bin/

install:
	install -v -d $(INSTALL_DIR) && install -v am $(INSTALL_DIR)


# The crlf-source target creates a set of source files for use on Windows.

crlf-source:
	$(MAKE) clean;\
	DEST=$$(pwd)-crlf;\
	rm -rf $${DEST};\
	mkdir $${DEST};\
	for f in *; do sed "s/$$//" $$f > $${DEST}/$$f; done


# If invoked directly, the 'am' target builds am with the default
# options here.  Recursive invocations build this target with modified
# macros that override these defaults.  The 'omp' and 'serial' targets
# point to the respective default build targets for the OpenMP and
# serial builds.

CC      = $(COMPILER_GCC)
CFLAGS  = $(CFLAGS_GCC_OPENMP)  $(OPT_FLAGS_GCC) $(EXTRA_CFLAGS)
LDFLAGS = $(LDFLAGS_GCC_OPENMP) $(EXTRA_LDFLAGS)
LDLIBS  = $(LDLIBS_GCC_OPENMP)  $(EXTRA_LDLIBS)


# User-definable extra flags from the command line.  Some examples:
#
#   Detailed optimization and vectorization reporting:
#     make gcc-omp "EXTRA_CFLAGS = -fopt-info"
#     make nvc-omp "EXTRA_CFLAGS = -Minfo -Mneginfo"
#     make icc-omp "EXTRA_CFLAGS = -qopt-report"
#
#   Multi-architecture support, useful on heterogeneous clusters:
#     make pgi-omp "EXTRA_CFLAGS = -tp x64"

EXTRA_CFLAGS  =
EXTRA_LDFLAGS = 
EXTRA_LDLIBS  = 


am: $(OBJS) version.o
	$(CC) -o am $(OBJS) version.o $(LDLIBS) $(LDFLAGS)

omp: gcc-omp

serial: gcc-serial


# Build targets for the gcc compiler.  Tested with gcc 14.2.0
# (2024 August 12).

COMPILER_GCC       = gcc

CFLAGS_GCC         = -std=c99 -Wall -W -Wshadow -Wwrite-strings $(EXTRA_CFLAGS)
CFLAGS_GCC_OPENMP  = $(CFLAGS_GCC) -fopenmp

LDFLAGS_GCC        = $(EXTRA_LDFLAGS)
LDFLAGS_GCC_OPENMP = $(LDFLAGS_GCC) 

LDLIBS_GCC         = -lm $(EXTRA_LDLIBS)
LDLIBS_GCC_OPENMP  = $(LDLIBS_GCC) -lgomp

OPT_FLAGS_GCC      = -O3

DEBUG_FLAGS_GCC    = -g -O0

gcc-omp:
	$(MAKE) am \
		"CC      = $(COMPILER_GCC)" \
		"CFLAGS  = $(CFLAGS_GCC_OPENMP) $(OPT_FLAGS_GCC)" \
		"LDFLAGS = $(LDFLAGS_GCC_OPENMP)" \
		"LDLIBS  = $(LDLIBS_GCC_OPENMP)"

gcc-serial:
	$(MAKE) am \
		"CC      = $(COMPILER_GCC)" \
		"CFLAGS  = $(CFLAGS_GCC) $(OPT_FLAGS_GCC)" \
		"LDFLAGS = $(LDFLAGS_GCC)" \
		"LDLIBS  = $(LDLIBS_GCC)"

gcc-no-opt:
	$(MAKE) am \
		"CC      = $(COMPILER_GCC)" \
		"CFLAGS  = $(CFLAGS_GCC)" \
		"LDFLAGS = $(LDFLAGS_GCC)" \
		"LDLIBS  = $(LDLIBS_GCC)"

gcc-no-opt-debug:
	$(MAKE) am \
		"CC      = $(COMPILER_GCC)" \
		"CFLAGS  = $(CFLAGS_GCC) $(DEBUG_FLAGS_GCC)" \
		"LDFLAGS = $(LDFLAGS_GCC)" \
		"LDLIBS  = $(LDLIBS_GCC)"


# Build targets for the Nvidia nvc compiler.  Tested with nvc 24.3-0
# (2024 July 24).  Note that OpenMP nested parallelism support was
# removed in pgcc/nvc following the acquisiton of PGI by Nvidia, and
# has not been restored as of this version.

COMPILER_NVC       = nvc

CFLAGS_NVC         = $(EXTRA_CFLAGS)
CFLAGS_NVC_OPENMP  = $(CFLAGS_NVC) -mp

LDFLAGS_NVC        = $(EXTRA_LDFLAGS)
LDFLAGS_NVC_OPENMP = $(LDFLAGS_NVC) -mp=allcores

LDLIBS_NVC         = $(EXTRA_LDLIBS)
LDLIBS_NVC_OPENMP  = $(LDLIBS_NVC)

OPT_FLAGS_NVC      = -fast -O3 -Msafeptr

nvc-omp:
	$(MAKE) am \
		"CC      = $(COMPILER_NVC)" \
		"CFLAGS  = $(CFLAGS_NVC_OPENMP) $(OPT_FLAGS_NVC)" \
		"LDFLAGS = $(LDFLAGS_NVC_OPENMP)" \
		"LDLIBS  = $(LDLIBS_NVC_OPENMP)"

nvc-serial:
	$(MAKE) am \
		"CC      = $(COMPILER_NVC)" \
		"CFLAGS  = $(CFLAGS_NVC) $(OPT_FLAGS_NVC)" \
		"LDFLAGS = $(LDFLAGS_NVC)" \
		"LDLIBS  = $(LDLIBS_NVC)"


# Build targets for the Intel icx compiler, the C driver for the
# Intel OneAPI compiler.  This replaces the deprecated icc compiler.
# Tested with OneAPI 2024.0.2 (2024 July 24).

COMPILER_ICX       = icx

CFLAGS_ICX         = $(EXTRA_CFLAGS)
CFLAGS_ICX_OPENMP  = $(CFLAGS_ICX) -qopenmp

LDFLAGS_ICX        = $(EXTRA_LDFLAGS)
LDFLAGS_ICX_OPENMP = $(LDFLAGS_ICX) -qopenmp

LDLIBS_ICX         = $(EXTRA_LDLIBS)
LDLIBS_ICX_OPENMP  = $(LDLIBS_ICX)

OPT_FLAGS_ICX      = -O3 -fno-alias -xHost

icx-omp:
	$(MAKE) am \
		"CC      = $(COMPILER_ICX)" \
		"CFLAGS  = $(CFLAGS_ICX_OPENMP) $(OPT_FLAGS_ICX)" \
		"LDFLAGS = $(LDFLAGS_ICX_OPENMP)" \
		"LDLIBS  = $(LDLIBS_ICX_OPENMP)"

icx-serial:
	$(MAKE) am \
		"CC      = $(COMPILER_ICX)" \
		"CFLAGS  = $(CFLAGS_ICX) $(OPT_FLAGS_ICX)" \
		"LDFLAGS = $(LDFLAGS_ICX)" \
		"LDLIBS  = $(LDLIBS_ICX)"


# hno3.o and o3.o are listed first among all the object files because
# they take longest to compile.  This saves time if make is invoked
# with the -j option to run parallel jobs.  Note that version.o is not
# included in this list, to avoid creating a circular dependency.

OBJS = hno3.o o3.o abscoeff.o am_alloc.o am_sysdep.o am_types.o\
 ch4.o ch3cn.o ch3oh.o clo.o cia.o co.o co2.o column.o config.o\
 dcache.o doc.o errlog.o fileops.o fit.o hbr.o hcl.o hcn.o h2co.o\
 hf.o h2o.o h2o_continuum.o h2o_ice.o h2o_liquid.o h2o_psat.o h2o2.o\
 h2s.o ho2.o hocl.o ils.o interp.o jacobian.o kcache.o layer.o\
 linesum.o main.o mapping.o model.o mt_ckd.o nh3.o nscale.o n2o.o\
 no.o no2.o o.o o2.o ocs.o oh.o oneline.o output.o planck.o\
 rayleigh.o rt.o simplex.o so2.o specfunc.o spectra.o tags.o\
 transform.o units.o

am: $(OBJS) version.o

$(OBJS): am_types.h version.h

version.o: $(OBJS) version.h

abscoeff.o: abscoeff.h am_alloc.h cia.h dcache.h errlog.h\
 h2o_continuum.h h2o_ice.h h2o_liquid.h kcache.h linesum.h output.h\
 phys_const.h units.h ch4.h ch3cn.h ch3oh.h co.h co2.h clo.h hbr.h\
 hcn.h hcl.h hf.h hno3.h h2co.h h2o.h h2o2.h ho2.h hocl.h h2s.h nh3.h\
 n2o.h no.h no2.h o.h o2.h o3.h ocs.h oh.h so2.h oneline.h

am_alloc.o: abscoeff.h am_alloc.h column.h errlog.h kcache.h layer.h\
 model.h molecules.h nscale.h output.h simplex.h units.h

am_sysdep.o: am_sysdep.h

am_types.o: abscoeff.h column.h fit.h ils.h layer.h model.h\
 molecules.h output.h phys_const.h units.h

cia.o: am_sysdep.h errlog.h math_const.h phys_const.h specfunc.h

column.o: abscoeff.h cia.h column.h errlog.h layer.h molecules.h\
 phys_const.h spectra.h units.h

config.o: abscoeff.h am_alloc.h column.h config.h errlog.h fileops.h\
 fit.h ils.h layer.h linesum.h mapping.h math_const.h model.h\
 nscale.h output.h phys_const.h simplex.h tags.h units.h

dcache.o: abscoeff.h dcache.h errlog.h fileops.h phys_const.h

doc.o: doc.h

errlog.o: errlog.h

fileops.o: am_sysdep.h fileops.h

fit.o: abscoeff.h am_alloc.h am_sysdep.h config.h dcache.h errlog.h\
 fileops.h fit.h jacobian.h kcache.h model.h nscale.h output.h\
 simplex.h spectra.h units.h

h2o_continuum.o: errlog.h h2o_continuum.h mt_ckd.h phys_const.h

h2o_ice.o: errlog.h h2o_ice.h math_const.h molecules.h phys_const.h\
 rayleigh.h

h2o_liquid.o: errlog.h h2o_liquid.h math_const.h molecules.h\
 phys_const.h rayleigh.h

h2o_psat.o: column.h errlog.h h2o_psat.h phys_const.h

ils.o: errlog.h ils.h math_const.h output.h phys_const.h specfunc.h\
 transform.h

jacobian.o: errlog.h jacobian.h mapping.h model.h output.h simplex.h

kcache.o: errlog.h kcache.h

layer.o: abscoeff.h am_alloc.h column.h errlog.h h2o_psat.h interp.h\
 layer.h mapping.h model.h nscale.h phys_const.h units.h

linesum.o: am_sysdep.h errlog.h linesum.h math_const.h phys_const.h\
 specfunc.h

main.o: am_alloc.h am_sysdep.h config.h dcache.h doc.h errlog.h fit.h\
 jacobian.h kcache.h model.h nscale.h output.h planck.h simplex.h\
 tags.h transform.h units.h

mapping.o: errlog.h mapping.h

model.o: am_sysdep.h column.h errlog.h ils.h interp.h layer.h\
 math_const.h model.h output.h phys_const.h rt.h spectra.h\
 units.h

mt_ckd.o: mt_ckd.h

nscale.o: nscale.h

output.o: abscoeff.h column.h dcache.h errlog.h fit.h ils.h kcache.h\
 layer.h linesum.h mapping.h model.h nscale.h output.h phys_const.h\
 simplex.h tags.h units.h

planck.o: am_sysdep.h phys_const.h planck.h specfunc.h

rt.o: am_sysdep.h column.h errlog.h layer.h model.h phys_const.h\
 planck.h specfunc.h spectra.h rt.h

rayleigh.o: errlog.h math_const.h phys_const.h rayleigh.h

simplex.o: errlog.h simplex.h

spectra.o: abscoeff.h column.h errlog.h math_const.h model.h output.h\
	phys_const.h planck.h specfunc.h transform.h units.h

transform.o: am_sysdep.h math_const.h transform.h

units.o: errlog.h units.h

ch4.o: ch4.h molecules.h

ch3cn.o: ch3cn.h molecules.h

ch3oh.o: ch3oh.h molecules.h

clo.o: clo.h molecules.h

co.o: co.h molecules.h

co2.o: co2.h molecules.h

hbr.o: hbr.h molecules.h

hcl.o: hcl.h molecules.h

hcn.o: hcn.h molecules.h

hf.o: hf.h molecules.h

h2co.o: h2co.h molecules.h

hno3.o: hno3.h molecules.h

h2o.o: h2o.h molecules.h

h2o2.o: h2o2.h molecules.h

ho2.o: ho2.h molecules.h

hocl.o: hocl.h molecules.h

h2s.o: h2s.h molecules.h

nh3.o: nh3.h molecules.h

n2o.o: n2o.h molecules.h

no.o: no.h molecules.h

no2.o: no2.h molecules.h

o.o: o.h molecules.h

o2.o: o2.h molecules.h

o3.o: o3.h molecules.h

ocs.o: ocs.h molecules.h

oh.o: oh.h molecules.h

so2.o: so2.h molecules.h

oneline.o: oneline.h molecules.h
