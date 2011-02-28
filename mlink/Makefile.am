
# Let automake guess compiler options for dummy target
EXTRA_PROGRAMS = PJFryDummy
PJFryDummy_SOURCES = PJFry.cpp
PJFryDummy_LDADD = ../src/libpjfry.la
EXTRA_PJFryDummy_SOURCES = PJFry.tm
PJFryDummy_CXXFLAGS = -Wno-write-strings

EXTRA_DIST = PJFry.cpp PJFry.tm
DISTCLEANFILES = PJFry.cpp
CLEANFILES = PJFry$(EXEEXT) PJFry_auto.cpp

if HAVE_MCC
     # Pretend that PJFry is a script
     nodist_bin_SCRIPTS = PJFry$(EXEEXT)
endif

# Trick automake into using libtool to install scripts
INSTALL_SCRIPT=$(INSTALL_PROGRAM_ENV) $(LIBTOOL) $(AM_LIBTOOLFLAGS) $(LIBTOOLFLAGS) --mode=install $(INSTALL_PROGRAM)

if HAVE_MPREP
PJFry.cpp: PJFry.tm
	$(MPREP) $? -o $@
endif

PJFRYSRC=$(shell if test -f 'PJFry.cpp'; then $(CYGPATH_W) 'PJFry.cpp'; else $(CYGPATH_W) '$(srcdir)/PJFry.cpp'; fi)
PJFry_auto.cpp: $(PJFRYSRC)
	ln -s -f $< $@

PJFry$(EXEEXT): CXXLD=$(MCC)
PJFry$(EXEEXT): PJFry_auto.cpp $(PJFryDummy_DEPENDENCIES)
	@rm -f PJFry$(EXEEXT)
	$(AM_V_CXXLD)$(PJFryDummy_LINK) $(DEFS) $(DEFAULT_INCLUDES) $(INCLUDES) $(AM_CPPFLAGS) $(CPPFLAGS) \
	-MT PJFryDummy-PJFry.o -MD -MP -MF $(DEPDIR)/PJFryDummy-PJFry.Tpo  $(PJFryDummy_LDADD) $(LIBS) $<
# 		`test -f '$<' || echo '$(srcdir)/'`$<
	$(AM_V_at)$(am__mv) $(DEPDIR)/PJFryDummy-PJFry.Tpo $(DEPDIR)/PJFryDummy-PJFry.Po
