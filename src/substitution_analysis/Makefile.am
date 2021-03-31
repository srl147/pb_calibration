bin_PROGRAMS = substitution_analysis

substitution_analysis_SOURCES = substitution_analysis.c

substitution_analysis_LDADD = @IO_LIB_LDFLAGS@ @SAMTOOLS_LDFLAGS@ -L../shared -lhelper_funcs

AM_CPPFLAGS = @IO_LIB_CFLAGS@ @SAMTOOLS_CFLAGS@ -I@top_srcdir@/shared
