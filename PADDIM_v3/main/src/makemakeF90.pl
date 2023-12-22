#! /usr/bin/perl
#
# Usage: makemake {<program name>}
#
# Generate a Makefile from the sources in the current directory.  The source
# files may be in either C, FORTRAN 77, Fortran 90 or some combination of
# these languages.  If the F90 compiler specified is cray or parasoft, then
# the Makefile generated will conform to the conventions of these compilers.
# To run makemake, it will be necessary to modify the first line of this script
# to point to the actual location of Perl on your system.
#
# Adapted from a perl code written by Michael Wester 
#                                       <wester@math.unm.edu> February 16, 1995
# 
open(MAKEFILE, "> Makefile");
#
print MAKEFILE "# +------------------------------------------------------------------------------+\n";
print MAKEFILE "# | Makefile for PADDI ([PA]rallel [D]ouble [DI]ffusion) Code                    |\n"; 
print MAKEFILE "# +------------------------------------------------------------------------------+\n";
print MAKEFILE "#\n";
print MAKEFILE "# To compile, the following flags can be set:\n";
print MAKEFILE "# -------------------------------------------\n";
print MAKEFILE "#\n";
print MAKEFILE "#   -DMPI_MODULE         Set if the mpi implementation provides a f90 MPI Module.\n";
print MAKEFILE "#                        (most implementations do nowadays...)                   \n";
print MAKEFILE "#   -AB_BDF3             Use AB_BDF3 time stepping scheme.\n";
print MAKEFILE "#                        (different time stepping scheme used to exist in pervious versions.)\n";
print MAKEFILE "#                        (Currently, only AB_BDF3 is supported, so always use this flag!)\n";
print MAKEFILE "#   -DTEMPERATURE_FIELD  Include Temperature field\n";
print MAKEFILE "#   -DCHEMICAL_FIELD     Include Chemical field\n";
print MAKEFILE "#   -DTWO_DIMENSIONAL    Build 2d code version\n";
print MAKEFILE "#                        (if not specified, the 3D version is generated)\n";
print MAKEFILE "#   -DSINGLE_PRECISION   Build code in single precision\n";
print MAKEFILE "#   -DDOUBLE_PRECISION   Build code in double precision\n";
print MAKEFILE "DEFS          = -DSINGLE_PRECISION -DMPI_MODULE -DAB_BDF3 -DTEMPERATURE_FIELD -DCHEMICAL_FIELD\n\n";
print MAKEFILE "FC            = mpif90 \n\n";
print MAKEFILE "F90           = \$(FC) \n\n";
print MAKEFILE "FFLAGS        = -O3 -xSSSE3 -align all -warn all \$(DEFS) \n\n";
print MAKEFILE "F90FLAGS      = -I../../stuff_needed/include -cpp -align all -O3 -xSSSE3 -warn all -stand f95 \$(DEFS) \n\n";
print MAKEFILE "LD	      = \$(FC)\n\n";
print MAKEFILE "LDFLAGS	      = -L../../stuff_needed/lib -xT\n\n";
print MAKEFILE "LIBS	      = \n\n";
print MAKEFILE "ADDLIBS       = -lfftw3 -lfftw3f -ljc -ljpeg -lpnetcdf \n\n";
#
# Source listing
#
print MAKEFILE "SRCS	      = ";
@srcs = <*.f90 *.f *.F *.c>;
&PrintWords(8, 1, @srcs);
print MAKEFILE "\n\n";
#
# Program name
#
print MAKEFILE "PROGRAM	      = $ARGV[0]\n\n";
#
# Object listing
#
print MAKEFILE "OBJS          = ";
@objs = @srcs;
foreach (@objs) { s/\.[^.]+$/.o/ };
&PrintWords(8, 1, @objs);
print MAKEFILE "\n\n";
#
# make
#
print MAKEFILE "all:\t \$(PROGRAM)\n\n";
print MAKEFILE "\$(PROGRAM):     \$(OBJS) \$(LIBS)\n";
print MAKEFILE "\t\t\@echo \"Linking \$(PROGRAM) ...\"\n";
print MAKEFILE "\t\t\@\$(LD) \$(LDFLAGS) \$(OBJS) \$(LIBS) \$(ADDLIBS) -o \$(PROGRAM)\n";
print MAKEFILE "\t\t\@echo \"done\"\n\n";
#
#print MAKEFILE "all: \$(PROG)\n\n";
#print MAKEFILE "\$(PROG): \$(OBJS)\n";
#print MAKEFILE "\t\$(", &LanguageCompiler(mpif90, @srcs);
#print MAKEFILE ") \$(LDFLAGS) -o \$@ \$(OBJS) \$(LIBS) \$(LDADD)\n\n";
#
# make clean
#
print MAKEFILE "clean:\n";
print MAKEFILE "\trm -f \$(PROGRAM) \$(OBJS) *.mod core\n\n";
#
# make clean
#
print MAKEFILE "totalclean:\n";
print MAKEFILE "\trm -f \$(PROG) \$(OBJS) *.mod *.~ *.txt *.gp *.out *.eps\n\n";
#
# Make .f90 a valid suffix
#
print MAKEFILE ".SUFFIXES: \$(SUFFIXES) .f90\n\n";
#
# .f90 -> .o
#
print MAKEFILE ".f90.o:\n";
print MAKEFILE "\t\$(F90) \$(F90FLAGS) -c \$<\n\n";
#
# stop interpreting .mod files as modular 2 stuff
#
print MAKEFILE ".mod.o:\n";
print MAKEFILE "\t\$(F90) \$(F90FLAGS) -c \$*.f90\n\n";
#
# Dependency listings
#
&MakeDependsf90(mpif90);
#&MakeDepends("*.f *.F", '^\s*include\s+["\']([^"\']+)["\']');
&MakeDepends("*.f90",   '^\s*include\s+["\']([^"\']+)["\']');
&MakeDepends("*.f90",   '^\s*#\s*include\s+["\']([^"\']+)["\']');
&MakeDepends("*.c",     '^\s*#\s*include\s+["\']([^"\']+)["\']');

#
# &PrintWords(current output column, extra tab?, word list); --- print words
#    nicely
#
sub PrintWords {
   local($columns) = 78 - shift(@_);
   local($extratab) = shift(@_);
   local($wordlength);
   #
   print MAKEFILE @_[0];
   $columns -= length(shift(@_));
   foreach $word (@_) {
      $wordlength = length($word);
      if ($wordlength + 1 < $columns) {
         print MAKEFILE " $word";
         $columns -= $wordlength + 1;
         }
      else {
         #
         # Continue onto a new line
         #
         if ($extratab) {
            print MAKEFILE " \\\n\t\t$word";
            $columns = 62 - $wordlength;
            }
         else {
            print MAKEFILE " \\\n\t$word";
            $columns = 70 - $wordlength;
            }
         }
      }
   }

#
# &LanguageCompiler(compiler, sources); --- determine the correct language
#    compiler
#
sub LanguageCompiler {
   local($compiler) = &toLower(shift(@_));
   local(@srcs) = @_;
   #
   if (length($compiler) > 0) {
      CASE: {
         grep(/^$compiler$/, ("fc", "f77")) &&
            do { $compiler = "FC"; last CASE; };
         grep(/^$compiler$/, ("cc", "c"))   &&
            do { $compiler = "CC"; last CASE; };
         $compiler = "F90";
         }
      }
   else {
      CASE: {
         grep(/\.f90$/, @srcs)   && do { $compiler = "f90"; last CASE; };
         grep(/\.(f|F)$/, @srcs) && do { $compiler = "FC";  last CASE; };
         grep(/\.c$/, @srcs)     && do { $compiler = "CC";  last CASE; };
         $compiler = "???";
         }
      }
   $compiler;
   }

#
# &toLower(string); --- convert string into lower case
#
sub toLower {
   local($string) = @_[0];
   $string =~ tr/A-Z/a-z/;
   $string;
   }

#
# &uniq(sorted word list); --- remove adjacent duplicate words
#
sub uniq {
   local(@words);
   foreach $word (@_) {
      if ($word ne $words[$#words]) {
         push(@words, $word);
         }
      }
   @words;
   }

#
# &MakeDepends(language pattern, include file sed pattern); --- dependency
#    maker
#
sub MakeDepends {
   local(@incs);
   local($lang) = @_[0];
   local($pattern) = @_[1];
   #
   foreach $file (<${lang}>) {
      open(FILE, $file) || warn "Cannot open $file: $!\n";
      while (<FILE>) {
         /$pattern/i && push(@incs, $1);
         }
      if (defined @incs) {
         $file =~ s/\.[^.]+$/.o/;
         print MAKEFILE "$file: ";
         &PrintWords(length($file) + 2, 0, @incs);
         print MAKEFILE "\n";
         undef @incs;
         }
      }
   }

#
# &MakeDependsf90(f90 compiler); --- FORTRAN 90 dependency maker
#
sub MakeDependsf90 {
   local($compiler) = &toLower(@_[0]);
   local(@dependencies);
   local(%filename);
   local(@incs);
   local(@modules);
   local($objfile);
   #
   # Associate each module with the name of the file that contains it
   #
   foreach $file (<*.f90>) {
      open(FILE, $file) || warn "Cannot open $file: $!\n";
      while (<FILE>) {
         /^\s*module\s+([^\s!]+)/i &&
            ($filename{&toLower($1)} = $file) =~ s/\.f90$/.o/;
         }
      }
   #
   # Print the dependencies of each file that has one or more include's or
   # references one or more modules
   #
   foreach $file (<*.f90>) {
      open(FILE, $file);
      while (<FILE>) {
         /^\s*include\s+["\']([^"\']+)["\']/i && push(@incs, $1);
         /^\s*use\s+([^\s,!]+)/i && push(@modules, &toLower($1));
         }
      if (defined @incs || defined @modules) {
         ($objfile = $file) =~ s/\.f90$/.o/;
         print MAKEFILE "$objfile: ";
         undef @dependencies;
         foreach $module (@modules) {
            push(@dependencies, $filename{$module});
            }
         @dependencies = &uniq(sort(@dependencies));
         &PrintWords(length($objfile) + 2, 0,
                     @dependencies, &uniq(sort(@incs)));
         print MAKEFILE "\n";
         undef @incs;
         undef @modules;
         #
         # Cray F90 compiler
         #
         if ($compiler eq "cray") {
            print MAKEFILE "\t\$(F90) \$(F90FLAGS) -c ";
            foreach $depend (@dependencies) {
               push(@modules, "-p", $depend);
               }
            push(@modules, $file);
            &PrintWords(30, 1, @modules);
            print MAKEFILE "\n";
            undef @modules;
            }
         #
         # ParaSoft F90 compiler
         #
         if ($compiler eq "parasoft") {
            print MAKEFILE "\t\$(F90) \$(F90FLAGS) -c ";
            foreach $depend (@dependencies) {
               $depend =~ s/\.o$/.f90/;
               push(@modules, "-module", $depend);
               }
            push(@modules, $file);
            &PrintWords(30, 1, @modules);
            print MAKEFILE "\n";
            undef @modules;
            }
         }
      }
   }
