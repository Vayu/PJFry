

check_PROGRAMS = demo

demo_LDADD = ../src/libpjfry.la
demo_SOURCES = demo.f
EXTRA_demo_SOURCES = pjfry.f
EXTRA_DIST = demo.txt

if EMU_GOLEM
  noinst_LTLIBRARIES = libgolem_demo_module.la libpjfry95pg_module.la
  libgolem_demo_module_la_SOURCES = golem_demo.f90

  nodist_libpjfry95pg_module_la_SOURCES = pjfry95pg.f90
  EXTRA_libpjfry95pg_module_la_SOURCES = pjfry95.py pjfry90.f90

  check_PROGRAMS += golem_4point golem_5point golem_6point

  golem_4point_LDADD = ../src/libpjfry.la libgolem_demo_module.la libpjfry95pg_module.la
  golem_4point_SOURCES = golem_4point.f90

  golem_5point_LDADD = ../src/libpjfry.la libgolem_demo_module.la libpjfry95pg_module.la
  golem_5point_SOURCES = golem_5point.f90

  golem_6point_LDADD = ../src/libpjfry.la libgolem_demo_module.la
  golem_6point_SOURCES = golem_6point.f90

pjfry95pg.f90: pjfry95.py
	python $<

DISTCLEANFILES=pjfry95pg.f90 pjfry95pj.f90

endif

clean-local:
	$(RM) *.mod
