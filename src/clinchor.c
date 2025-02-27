/*************************************************************************\
 * Copyright (c) 2002 The University of Chicago, as Operator of Argonne
 * National Laboratory.
 * Copyright (c) 2002 The Regents of the University of California, as
 * Operator of Los Alamos National Laboratory.
 * This file is distributed subject to a Software License Agreement found
 * in the file LICENSE that is included with this distribution. 
\*************************************************************************/

/* 
   $Log: clinchor.c,v $
   Revision 1.25  2011/06/23 08:38:16  lemery
   Removed a debugging print statement.

   Revision 1.24  2011/04/27 04:34:28  lemery
   Changed a long declaration to int32_t to fix a problem with setting a value
   by reference in a SDDS library.

   Revision 1.23  2011/04/27 02:17:00  lemery
   Shang added macro capability for the command line.
   Emery removed the unnecessary cm_show command for longitudinal plane.

   Revision 1.22  2010/10/08 14:59:39  soliday
   Fixed an issue on 64bit computers.

   Revision 1.21  2010/05/19 15:23:26  soliday
   Updated to fix a few compiler warnings.

   Revision 1.20  2010/04/13 20:27:52  lemery
   Added comment that function cg could be replaced by a (presumed optimized) lapack function.
  
   Revision 1.19  2010/04/13 16:03:40  soliday
   All of the complex number code in the SDDS libraries has been ported to C++.
   Until this code can also be ported I have added the needed code to get this
   to compile again.
  
   Revision 1.18  2007/09/13 17:36:38  emery
   Forgot to add a SDDS_FindColumn call for the new column "Fixed."
  
   Revision 1.17  2007/09/13 17:31:05  emery
   Added a column Fixed for the HOM input files
   to allow frequencies to be fixed, i.e. not staggered or randomized.
  
   Revision 1.16  2007/09/13 16:31:11  emery
   Removed the default case from the return test of Frequency column check units.
  
   Revision 1.15  2007/09/13 16:20:32  emery
   Added better warning message forHOM  Frequency column units.
  
   Revision 1.14  2007/06/21 06:11:38  emery
   Corrected the form factor for longitudinal growth as well.
  
   Revision 1.13  2007/05/21 05:58:41  emery
   Regrouped terms in complex exponentials to prevent overflow
   in numerators when Q is low.
  
   Revision 1.12  2007/05/17 08:05:20  emery
   Removed the factor two in the bunch form factor.
   Wrote comment explaining the physical reason for this.
  
   Revision 1.11  2006/11/03 20:39:08  emery
   Added namelist command randomizeBunchPattern.
   Changed the interpretation of "spread" for randomized quantities.
   Spread should be the expected difference between minimum and
   maximum possible values. I had previously made the spread equal
   to the difference between the original value and the extreme values.
  
   Revision 1.10  2006/10/26 03:19:33  emery
   Change variable name bunchLengthByBunch with name bunchDurationByBunch
   to make the point that the bunch length quantity has time dimension.
  
   Revision 1.9  2006/09/15 22:55:24  emery
   Added bunch length data file to ring parameter for reading so that
   individual bunches get a possible different length for form-factor
   calculation. Bunch length used are written to file.
   Fixed checks on units in HOM file read-in.
   Fixed output to file of eigenvectors. Form factor had an error
   of 2 in the argument.
  
   clinchor program. Louis Emery ANL 1992
   * The general organization and some code is taken from program elegant
   * written by M. Borland
   * file clinchor.c (March 1995) is the version with SDDS input and output.
   */
#include "mdb.h"
#include "scan.h"
#include "namelist.h"
#include "match_string.h"
#include "clinchor.h" /* generated from clinchor.nl */
#include "constants.h"
#include "SDDS.h"
#if defined(MKL)
#include "mkl.h"
#endif
#if defined(LAPACK)
#include "zgeev.c"
#endif

typedef struct {
  double r, i;
} COMPLEX, *COMPTR;
typedef struct {
  double **ar, **ai;
  int n, m;
} CMATRIX;

#define SWAP_INT32(x, y) {int32_t tmp_swap_int32; tmp_swap_int32=(x); (x)=(y); (y)=tmp_swap_int32; }

#define USAGE "Usage: clinchor <inputFile> [-verbose] [-macro=<tag>=<value>] [-threads=<num>]"
void loadDataFromTwissFile() ;
void loadDataFromBunchLengthFile() ;


#define CLO_VERBOSE 0
#define CLO_MACRO 1
#define CLO_THREADS 2
#define CLO_OPTIONS 3

char *option[CLO_OPTIONS] = {
  "verbose", "macro", "threads"
};

#define MONOPOLE 0
#define DIPOLE 1
#define MULTIPOLE_TYPES 2
char *multipoleTypeList[MULTIPOLE_TYPES]={
  "monopole",
  "dipole"};

#define RING_PARAMETERS 0
#define SYMMETRIC_BUNCH_PATTERN 1
#define BUNCH_TRAIN 2
#define GENERAL_BUNCH_PATTERN 3
#define MONOPOLE_HOMS 4
#define DIPOLE_HOMS 5
#define SWEEP_HOM_FREQUENCY 6
#define RANDOMIZE_HOM_FREQUENCIES 7
#define DO_LONGITUDINAL_MOTION 8
#define DO_TRANSVERSE_MOTION 9
#define STOP 10
#define BUNCH_PATTERN_FROM_FILE 11
#define RANDOMIZE_BUNCH_CURRENT 12
#define SEMAPHORES 13
#define COMMANDS 14
char *commandNames[COMMANDS]={
  "ringParameters",
  "symmetricBunchPattern",
  "bunchTrain",
  "generalBunchPattern",
  "monopoleHOMs",
  "dipoleHOMs",
  "sweepHOMFrequency",
  "randomizeHOMFrequencies",
  "doLongitudinalMotion",
  "doTransverseMotion",
  "stop",
  "bunchPatternFromFile",
  "randomizeBunchCurrent",
  "semaphores",
};

/* The variable array description does not have to follow the same order
the above arrays.  They just have to have the same number of entries */
char *description[COMMANDS] = {
  "ringParameters        define some ring parameters",
  "symmetricBunchPattern define a pattern of symmetrically placed\n\
                       and equally charged bunches",
  "bunchTrain            define a pattern of equally spaced and equally\n\
                       charged bunches",
  "generalBunchPattern   define a general bunch pattern and current",
  "bunchPatternFromFile  define multiply bunch patterns from a file which\n\
                       contains bucket and current columns",
  "randomizeBunchPattern define multiply bunch patterns from a file which\n\
                       contains bucket and current columns",
  "monopoleHOMs          read in a file describing monopole HOMs",
  "dipoleHOMs            read in a file describing dipole HOMs",
  "sweepHOMFrequency     allows one to select a single HOM type whose\n\
                       resonant frequency will be swept",
  "randomizeHOMFrequencies adds a random component to the resonant\n\
                       frequencies of any number of HOMS.",
  "doLongitudinalMotion  action command to calculated normal modes and\n\
                       growth rates in the longitudinal plane",
  "doTransverseMotion    action command to calculated normal modes and\n\
                       growth rates in the transverse plane",
  "semaphores          use 0-length output files to indicate status of job",
  "stop                  command to stop reading the input file and\n\
                       to stop excution"
};

typedef struct {
  double unperturbedFrequency, shuntImpedance, Q, RoQ;
  double *frequency;
  double *power,*powerPlus,*powerMinus; /* power estimated from the Vb of bunch 0 */
  double deQFactor;
  int32_t cavities;
  double staggeringStep;
  int32_t shiftToResonance;
  int32_t fixed;
} HOM_RESONATOR;


HOM_RESONATOR *monopoleHOMResonator;
long monopoleHOMTypes=0;
HOM_RESONATOR *dipoleHOMResonator;
long dipoleHOMTypes=0;
long totalFilledBuckets=0;
int32_t *bucketSelected;
/* All currents in variables and in SDDS file shall be mA unit (rather than A), 
   because most rings use mA for current units. */
double *bunchCurrent;
long sweepFlag=0;
long randomizeHOMFlag=0;
long randomizeBunchCurrentFlag=0;
long printEigenvectorsFlag=1;
long verbose;

/* some physics quantities */
double synchrotronFrequency, synchrotronAngFrequency;
double revolutionFrequency, revolutionAngFrequency;
double tune, betatronFrequency, betatronAngFrequency;
double betaAtRFCavs;
double synchronousPhase;
double rfAcceptance; /* in fraction */
/* bunch length is used to compute the form factor for the induced
   voltage of HOM. bunchLength has dimension of
   length. bunchDurationByBunch has dimension of time */
double bunchLength, *bunchDurationByBunch;
double rfFrequency;
double dampingFactor;
lapack_complex_double *lapack_eigenvalues;
double beamCurrentMA = 0;

long bunchLengthValues;
double *bunchCurrentData, *bunchChargeData, *bunchLengthPicoSecondData;
char *bunchLengthFile;

void DefaultPhysicalParameters();
void calculateSomePhysicalParameters();
void setupSymmetricBunchPattern(NAMELIST_TEXT *namelistText);
void setupBunchTrain(NAMELIST_TEXT *namelistText);
void setupGeneralBunchPattern(NAMELIST_TEXT *namelistText);
void setupBunchPatternFromFile(NAMELIST_TEXT *namelistText);
void setupRandomizeBunchCurrent(NAMELIST_TEXT *namelistText);
void doRandomizeBunchCurrent();
void checkForBucketDuplication();
void determineBunchLengthByBunch();
void printBunchLengthsUsed ();
void inputHOMs(NAMELIST_TEXT *namelistText, char *multipoleType);
void setupSweepHOMFrequency(NAMELIST_TEXT *namelistText);
void setupRandomizeHOMFrequencies(NAMELIST_TEXT *namelistText);
void calculateLongitudinalGrowthRate(NAMELIST_TEXT *namelistText);
void calculateTransverseGrowthRate(NAMELIST_TEXT *namelistText);
double powerOfResonators();
COMPLEX frequencyOfMaxLongGrowthRate();
COMPLEX frequencyOfMaxTransGrowthRate();

void do_semaphore_setup(char **semaphoreFile, NAMELIST_TEXT *nltext, char *rootname);
void createSemaphoreFile(char *filename);
static char *semaphoreFile[3];
void bombClinchor(char *errorMessage);
char *compose_filename(char *template, char *root_name);


COMPLEX csub(COMPLEX a, COMPLEX b)
{
  static COMPLEX c;
  c.r = a.r - b.r;
  c.i = a.i - b.i;
  return(c);
}
COMPLEX cmul(COMPLEX a, COMPLEX b)
{
  static COMPLEX c;
  static double cr, ci;
  cr = a.r*b.r - a.i*b.i;
  ci = a.r*b.i + a.i*b.r;
  c.r = cr;
  c.i = ci;
  return(c);
}
COMPLEX cconj(COMPLEX a)
{
  static COMPLEX b;
  b.r = a.r;
  b.i = -a.i;
  return(b);
}
COMPLEX cmulr(COMPLEX a, double r)
{
  static COMPLEX b;
  b.r = a.r*r;
  b.i = a.i*r;
  return(b);
}
COMPLEX cadd(COMPLEX a, COMPLEX b)
{
  static COMPLEX c;
  c.r = a.r + b.r;
  c.i = a.i + b.i;
  return(c);
}
COMPLEX cexp_oag(COMPLEX a)
{
  static COMPLEX b;
  static double exp_x;
  exp_x = exp(a.r);
  b.r  = exp_x*cos(a.i);
  b.i  = exp_x*sin(a.i);
  return(b);
}
/*
  COMPLEX cdiv(COMPLEX a, COMPLEX b)
  {
  static COMPLEX c;
  static double b_mod, cr, ci;
  if ((b_mod = b.r*b.r + b.i*b.i) == 0.0) {
  bombClinchor("division by zero in cdiv()");
  }
  cr = (a.r*b.r + a.i*b.i)/b_mod;
  ci = (-a.r*b.i + a.i*b.r)/b_mod;
  c.r = cr;
  c.i = ci;
  return(c);
  }
*/
double cmod(COMPLEX a)
{
  return(sqrt(sqr(a.r)+sqr(a.i)));
}

#if !defined(MKL)
void vdExp(const int n, const double *a, double *y) {
  int i;
  for (i=0;i<n;i++)
    y[i]=exp(a[i]);
}
void vdSinCos(const int n, const double *a, double *y, double *z) {
  int i;
  for (i=0;i<n;i++) {
    y[i]=sin(a[i]);
    z[i]=cos(a[i]);
  }
}
#endif

long threads=1;

int main( int argc, char **argv)
{
  SCANNED_ARG *scanned;
  char s[1024];
  FILE *fpIn = NULL;
  char *inputfile;
  long i,j;
  char *parameterFile;
  SDDS_DATASET parameterDataSet;
  long macros;
  char **macroTag, **macroValue;
  char *rootname;
  char *ptr;
  NAMELIST_TEXT namelist_text;

  rootname = NULL;
  semaphoreFile[0] = semaphoreFile[1] = semaphoreFile[2] = NULL;

  macros = 0;
  macroTag = macroValue = NULL;
#if defined(VAX_VMS) || defined(UNIX) || defined(SUN4)
  init_stats();
#endif
  argc = scanargs(&scanned, argc, argv);
  if (argc<2 || argc>(2+CLO_OPTIONS))
    bomb(NULL, USAGE);

  inputfile = NULL;
  for (i=1; i<argc; i++) {
    if (scanned[i].arg_type==OPTION) {
      switch( match_string(scanned[i].list[0], option , CLO_OPTIONS, 0)) {
      case CLO_VERBOSE:
        verbose=1;
        break;
      case CLO_MACRO:
	if ((scanned[i].n_items-=1)<0) {
          bombClinchor("Invalid -macro syntax");
        }
        if (!(macroTag=SDDS_Realloc(macroTag, sizeof(*macroTag)*(macros+scanned[i].n_items))) ||
            !(macroValue=SDDS_Realloc(macroValue, sizeof(*macroValue)*(macros+scanned[i].n_items))))
          bomb("memory allocation failure (-macro)", NULL);
        else {
          for (j=0; j<scanned[i].n_items; j++) {
            macroTag[macros] = scanned[i].list[j+1];
            if (!(macroValue[macros] = strchr(macroTag[macros], '='))) {
              bombClinchor("Invalid -macro syntax");
            }
            macroValue[macros][0] = 0;
            macroValue[macros] += 1;
            macros++;
          }
        }
        scanned[i].n_items +=1;
	break;
      case CLO_THREADS:
	if (scanned[i].n_items!=2) {
          bombClinchor("Invalid -threads syntax");
        }
        if (sscanf(scanned[i].list[1], "%ld", &threads)!=1 || threads<=0) {
          bombClinchor("Invalid -threads syntax");
        }
	break;
      default:
        bomb("unknown option given.", USAGE);
        break;
      }
    }
    else {
      if (!inputfile)
        fpIn = fopen_e(inputfile = scanned[i].list[0], "r", 0);
      else
        bomb("Too many files listed", USAGE);
    }
  }
  
  if (!inputfile)
    bomb("No input file given", USAGE);

  /* extract the root filename from the input filename */
  strcpy_ss(s, inputfile);
  if (rootname==NULL) {
    clean_filename(s);
    if ((ptr=strrchr(s, '.')))
      *ptr = 0;
    cp_str(&rootname, s);
  }

  DefaultPhysicalParameters();
  while( get_namelist(s, 1024, fpIn) ) {
    substituteTagValue(s, 1024, macroTag, macroValue, macros);
    scan_namelist(&namelist_text, s);
    switch( match_string(namelist_text.group_name, commandNames, COMMANDS, EXACT_MATCH) ) {
    case RING_PARAMETERS:
      set_namelist_processing_flags(0);
      twissFile = NULL;
      process_namelist(&ringParameters, &namelist_text);
      print_namelist(stdout, &ringParameters);
      /* check for validity of input */
      if (!twissFile) {
        if (!energyGeV)
          bomb("Variable energyGeV not specified or was set to zero", NULL);
        if (!circumference)
          bomb("Variable circumference not specified or was set to zero", NULL);
        if (!rfVoltageMV)
          bomb("Variable rfVoltageMV not specified or was set to zero", NULL);
        if (!momentumCompaction)
          bomb("Variable momentumCompaction not specified or was set to zero", NULL);
        if (!bunchLengtheningFactor)
          bomb("Variable bunchLengtheningFactor not specified or was set to zero", NULL);
      } else {
        loadDataFromTwissFile();

      }
      
      if (bunchLengthTableFile) {
        bunchLengthTableFile = compose_filename(bunchLengthTableFile, inputfile);
        loadDataFromBunchLengthFile();
      }
      
      calculateSomePhysicalParameters();
      /* write to parameter file */
      if (parameters) {
        /* note that parameters is a keyword in a namlist command. One must
           make a file name with a different memory if one wants to expand
           the %s string into something else. Same for other files in a
           namelist command.*/
        parameterFile  = compose_filename(parameters, inputfile);
        if (!SDDS_InitializeOutput(&parameterDataSet, SDDS_BINARY, 1,
                                   "Physical parameters of clinchor run",
                                   "Physical parameters", parameterFile))
          SDDS_PrintErrors(stderr, SDDS_EXIT_PrintErrors|SDDS_VERBOSE_PrintErrors);
        if (0>SDDS_DefineParameter(&parameterDataSet, "Energy", "E", "GeV",
                                   NULL, NULL, SDDS_DOUBLE, 0) ||
            0>SDDS_DefineParameter(&parameterDataSet, "Circumference", "C", "m",
                                   NULL, NULL, SDDS_DOUBLE, 0) ||
            0>SDDS_DefineParameter(&parameterDataSet, "RFVoltage", "V$bRF$a", "MV",
                                   NULL, NULL, SDDS_DOUBLE, 0) ||
            0>SDDS_DefineParameter(&parameterDataSet, "MomentumCompaction", "$ga$r$bc$n", NULL,
                                   NULL, NULL, SDDS_DOUBLE, 0)||
            0>SDDS_DefineParameter(&parameterDataSet, "BunchLengtheningFactor", "$gs/s$r$b0$n", NULL,
                                   NULL, NULL, SDDS_DOUBLE, 0)||
            0>SDDS_DefineParameter(&parameterDataSet, "EnergySpread", "$gs$be$r$n", NULL,
                                   NULL, NULL, SDDS_DOUBLE, 0)||
            0>SDDS_DefineParameter(&parameterDataSet, "LongDampingRate", "$gt$be$r$n", NULL,
                                   NULL, NULL, SDDS_DOUBLE, 0)||
            0>SDDS_DefineParameter(&parameterDataSet, "HorTransDampingRate", "$gt$r$bx$n", NULL,
                                   NULL, NULL, SDDS_DOUBLE, 0)||
            0>SDDS_DefineParameter(&parameterDataSet, "VertTransDampingRate", "$gt$r$by$n", NULL,
                                   NULL, NULL, SDDS_DOUBLE, 0)||
            0>SDDS_DefineParameter(&parameterDataSet, "HarmonicNumber", "h", NULL,
                                   NULL, NULL, SDDS_LONG, 0)||
            0>SDDS_DefineParameter(&parameterDataSet, "HorizontalTune", "$gn$r$bx$n", NULL,
                                   NULL, NULL, SDDS_DOUBLE, 0)||
            0>SDDS_DefineParameter(&parameterDataSet, "VerticalTune", "$gn$r$by$n", NULL,
                                   NULL, NULL, SDDS_DOUBLE, 0)||
            0>SDDS_DefineParameter(&parameterDataSet, "BetaXatRfCavities", "$gb$r$bx,RF$n", NULL,
                                   NULL, NULL, SDDS_DOUBLE, 0)||
            0>SDDS_DefineParameter(&parameterDataSet, "BetaYatRfCavities", "$gb$r$by,RF$n", NULL,
                                   NULL, NULL, SDDS_DOUBLE, 0)||
            0>SDDS_DefineParameter(&parameterDataSet, "RevolutionFrequency", "f$brev$n", NULL,
                                   NULL, NULL, SDDS_DOUBLE, 0)||
            0>SDDS_DefineParameter(&parameterDataSet, "RFFrequency", "f$bRF$n", NULL,
                                   NULL, NULL, SDDS_DOUBLE, 0)||
            0>SDDS_DefineParameter(&parameterDataSet, "SynchronousPhase", "$gf$r$bs$n", NULL,
                                   NULL, NULL, SDDS_DOUBLE, 0)||
            0>SDDS_DefineParameter(&parameterDataSet, "SynchrotronFrequency", "f$bs$n", NULL,
                                   NULL, NULL, SDDS_DOUBLE, 0)||
            0>SDDS_DefineParameter(&parameterDataSet, "SynchrotronTune", "$gn$r$bs$n", NULL,
                                   NULL, NULL, SDDS_DOUBLE, 0)||
            0>SDDS_DefineParameter(&parameterDataSet, "RFacceptance", NULL, NULL,
                                   NULL, NULL, SDDS_DOUBLE, 0)||
            0>SDDS_DefineParameter(&parameterDataSet, "BunchLength", "$gs$r$bz$n", NULL,
                                   NULL, NULL, SDDS_DOUBLE, 0)||
            !SDDS_WriteLayout(&parameterDataSet))
          SDDS_PrintErrors(stderr, SDDS_EXIT_PrintErrors|SDDS_VERBOSE_PrintErrors);
        if (!SDDS_StartPage(&parameterDataSet, 0))
          SDDS_PrintErrors(stderr, SDDS_EXIT_PrintErrors|SDDS_VERBOSE_PrintErrors);
        if (!SDDS_SetParameters(&parameterDataSet, SDDS_BY_NAME|SDDS_PASS_BY_VALUE,
                                "Energy", energyGeV,
                                "Circumference", circumference,
                                "RFVoltage", rfVoltageMV,
                                "MomentumCompaction", momentumCompaction,
                                "BunchLengtheningFactor", bunchLengtheningFactor,
                                "EnergySpread", relativeEnergySpread,
                                "LongDampingRate", longDampingRate,
                                "HorTransDampingRate", horTransDampingRate,
                                "VertTransDampingRate", vertTransDampingRate,
                                "HarmonicNumber", harmonicNumber,
                                "HorizontalTune", horizontalTune,
                                "VerticalTune", verticalTune,
                                "BetaXatRfCavities", betaxAtRFCavities,
                                "BetaYatRfCavities", betayAtRFCavities,
                                "RevolutionFrequency", revolutionFrequency,
                                "RFFrequency", rfFrequency,
                                "SynchronousPhase", synchronousPhase,
                                "SynchrotronFrequency", synchrotronFrequency,
                                "SynchrotronTune", synchrotronFrequency/revolutionFrequency,
                                "RFacceptance", rfAcceptance,
                                "BunchLength", bunchLength, NULL)||
            !SDDS_WritePage(&parameterDataSet) || !SDDS_Terminate(&parameterDataSet))
          SDDS_PrintErrors(stderr, SDDS_EXIT_PrintErrors|SDDS_VERBOSE_PrintErrors);
      }
      /* create a bunch length file name for later use when writing bunch length
         data */
      if (bunchLengthUsedFile) {
        bunchLengthFile = compose_filename(bunchLengthUsedFile, inputfile);
      }
      
      break;
    case SYMMETRIC_BUNCH_PATTERN:
      setupSymmetricBunchPattern(&namelist_text);
      break;
    case BUNCH_TRAIN:
      setupBunchTrain(&namelist_text);
      break;
    case GENERAL_BUNCH_PATTERN:
      setupGeneralBunchPattern(&namelist_text);
      break;
    case BUNCH_PATTERN_FROM_FILE:
      setupBunchPatternFromFile(&namelist_text);
      break;
    case RANDOMIZE_BUNCH_CURRENT:
      setupRandomizeBunchCurrent(&namelist_text);
      break;
    case MONOPOLE_HOMS:
      inputHOMs(&namelist_text, "monopole");
      break;
    case DIPOLE_HOMS:
      inputHOMs(&namelist_text, "dipole");
      break;
    case SWEEP_HOM_FREQUENCY:
      setupSweepHOMFrequency(&namelist_text);
      break;
    case RANDOMIZE_HOM_FREQUENCIES:
      setupRandomizeHOMFrequencies(&namelist_text);
      break;
    case DO_LONGITUDINAL_MOTION:
      calculateLongitudinalGrowthRate(&namelist_text);
#if defined(VAX_VMS) || defined(UNIX) || defined(SUN4)
      report_stats(stdout, "statistics:");
      fflush(stdout);
#endif
      break;
    case DO_TRANSVERSE_MOTION:
      calculateTransverseGrowthRate(&namelist_text);
#if defined(VAX_VMS) || defined(UNIX) || defined(SUN4)
      report_stats(stdout, "statistics:");
      fflush(stdout);
#endif
      break;
    case SEMAPHORES:
      /* if one is going to use semaphore files, then the namelist
	 has to be the first one in the command file.
	 There's no automatic enforcement of this. 
	 This functions writes to string array semaphoreFile */
      do_semaphore_setup(semaphoreFile, &namelist_text, rootname);
      /* start */
      if (semaphoreFile[0]) 
	createSemaphoreFile(semaphoreFile[0]);
      break;
    case STOP:
      /* done */
      if (semaphoreFile[1]) 
	createSemaphoreFile(semaphoreFile[1]);
      exit(0);
      break;
    default:
      printf("Unknown namelist %s given. Known namelists are:\n", namelist_text.group_name);
      for (i=0; i<COMMANDS; i++)
        printf("%s\n", description[i]);
      bombClinchor(NULL);
      break;
    }              
  }
#if defined(VAX_VMS) || defined(UNIX) || defined(SUN4)
  report_stats(stdout, "statistics:");
  fflush(stdout);
#endif
  /* done */
  if (semaphoreFile[1]) 
    createSemaphoreFile(semaphoreFile[1]);
  exit(0);
}

void inputHOMs(NAMELIST_TEXT *namelistText, char *multipoleType) {
  long offset;
  HOM_RESONATOR *resonator = NULL;
  SDDS_DATASET inputDataSet;
  long iDataSet, DataSetRows, row;
  char *FrequencyColumnName;
  char *ImpedanceColumnName;
  char *QColumnName;
  char *RQColumnName;
  char *deQFactorColumnName;
  char *StaggeringStepColumnName;
  char *ResonanceColumnName;
  char *FixedColumnName;
  char *CavitiesParameterName;
  char *filename;
  double frequencyConversionFactor;
  int32_t cavities;
  
  FrequencyColumnName=NULL;
  ImpedanceColumnName=NULL;
  QColumnName=NULL;
  RQColumnName=NULL;
  deQFactorColumnName =NULL;
  CavitiesParameterName=NULL;
  StaggeringStepColumnName=NULL;
  ResonanceColumnName=NULL;
  FixedColumnName=NULL;
  filename=NULL;
  frequencyConversionFactor = 1;

  
  switch(match_string(multipoleType, multipoleTypeList, MULTIPOLE_TYPES, 0)){
  case MONOPOLE:
    monopoleHOMs_struct.clearPreviousMonopoleHOMs = 0;
    set_namelist_processing_flags(0);
    process_namelist(&monopoleHOMs, namelistText);
    print_namelist(stdout, &monopoleHOMs);
    filename=monopoleHOMs_struct.filename;
    break;
  case DIPOLE:
    dipoleHOMs_struct.clearPreviousDipoleHOMs = 0;
    set_namelist_processing_flags(0);
    process_namelist(&dipoleHOMs, namelistText);
    print_namelist(stdout, &dipoleHOMs);
    filename=dipoleHOMs_struct.filename;
    break;
  default:
    bomb("unknown multipole type in function inputHOMs.", NULL);
  }
  
  if (!SDDS_InitializeInput(&inputDataSet, filename))
    SDDS_PrintErrors(stderr, SDDS_EXIT_PrintErrors|SDDS_VERBOSE_PrintErrors);
  while(0<(iDataSet=SDDS_ReadPage(&inputDataSet))){
    if (iDataSet==1){
      if (!(FrequencyColumnName=SDDS_FindColumn(&inputDataSet, FIND_NUMERIC_TYPE, "Frequency", "f", NULL))){
        bombClinchor("Problem with frequency column. Allowed names are Frequency or f");
      }
      switch(SDDS_CheckColumn(&inputDataSet, FrequencyColumnName, "Hz", SDDS_ANY_FLOATING_TYPE, NULL)){
      case SDDS_CHECK_WRONGUNITS:
        switch(SDDS_CheckColumn(&inputDataSet, FrequencyColumnName, "MHz", SDDS_DOUBLE, stderr)){
        case SDDS_CHECK_OKAY:
          frequencyConversionFactor=1e6;
          break;
        case SDDS_CHECK_WRONGUNITS:
          bombClinchor("Invalid value for units. Permitted values are MHz and Hz.");
        }
      case SDDS_CHECK_OKAY:
        break;
      default:
        fprintf(stderr,"Something wrong with column %s in file %s. Wrong type?\n", FrequencyColumnName, filename);
        bombClinchor(NULL);
      }
      ImpedanceColumnName=SDDS_FindColumn(&inputDataSet, FIND_NUMERIC_TYPE,
                                          "R", "Rs", "Rt", "ShuntImpedance", NULL);
      if (ImpedanceColumnName){
        switch(match_string(multipoleType, multipoleTypeList, MULTIPOLE_TYPES, 0)){ 
        case MONOPOLE:
          switch(SDDS_CheckColumn(&inputDataSet, ImpedanceColumnName, "Ohm", SDDS_DOUBLE, stderr)){
          case SDDS_CHECK_OKAY:
            break;
          case SDDS_CHECK_WRONGUNITS:
            switch(SDDS_CheckColumn(&inputDataSet, ImpedanceColumnName, "", SDDS_DOUBLE, stderr)){
            case SDDS_CHECK_OKAY:
              break;
            case SDDS_CHECK_WRONGUNITS:
              bombClinchor("Invalid value for units. Permitted value is Ohm.");
              break;
            }
            break;
          default:
            printf("Something wrong with column %s in file %s. Wrong type?\n", ImpedanceColumnName, filename);
            bombClinchor(NULL);
          }
          break;
        case DIPOLE:
          switch(SDDS_CheckColumn(&inputDataSet, ImpedanceColumnName, "Ohm/m", SDDS_DOUBLE, stderr)){
          case SDDS_CHECK_OKAY:
            break;
          case SDDS_CHECK_WRONGUNITS:
            switch(SDDS_CheckColumn(&inputDataSet, ImpedanceColumnName, "", SDDS_DOUBLE, stderr)){
            case SDDS_CHECK_OKAY:
              break;
            case SDDS_CHECK_WRONGUNITS:
              bombClinchor("Invalid value for units. Permitted values are Ohm and Ohm/m.");
              break;
            }
            break;
          default:
            printf("Something wrong with column %s in file %s. Wrong type?\n", ImpedanceColumnName, filename);
            bombClinchor(NULL);
          }
          break;
        }
      }
      
      RQColumnName=SDDS_FindColumn(&inputDataSet, FIND_NUMERIC_TYPE, "R/Q", "ROQ", "RoQ", "ROverQ","RoverQ", NULL);
      if (RQColumnName){
        switch(match_string(multipoleType, multipoleTypeList, MULTIPOLE_TYPES, 0)){ 
	case MONOPOLE:
	  switch(SDDS_CheckColumn(&inputDataSet, RQColumnName, "Ohm", SDDS_DOUBLE, stderr)){
	  case SDDS_CHECK_OKAY:
	    break;
	  case SDDS_CHECK_WRONGUNITS:
	    switch(SDDS_CheckColumn(&inputDataSet, RQColumnName, "", SDDS_DOUBLE, stderr)){
	    case SDDS_CHECK_OKAY:
	      break;
	    case SDDS_CHECK_WRONGUNITS:
	      bombClinchor("Invalid value for units. Permitted value is Ohm.");
	      break;
	    }
	    break;
	  default:
	    printf("Something wrong with column %s in file %s. Wrong type?\n", RQColumnName, filename);
	    bombClinchor(NULL);
	  }
	  break;
	case DIPOLE:
	  switch(SDDS_CheckColumn(&inputDataSet, RQColumnName, "Ohm/m", SDDS_DOUBLE, stderr)){
	  case SDDS_CHECK_OKAY:
	    break;
	  case SDDS_CHECK_WRONGUNITS:
	    switch(SDDS_CheckColumn(&inputDataSet, RQColumnName, "", SDDS_DOUBLE, stderr)){
	    case SDDS_CHECK_OKAY:
	      break;
	    case SDDS_CHECK_WRONGUNITS:
	      bombClinchor("Invalid value for units. Permitted value is Ohm/m.");
	      break;
	    }
	    break;
	  default:
	    printf("Something wrong with column %s in file %s. Wrong type?\n", RQColumnName, filename);
	    bombClinchor(NULL);
	  }
	  break;
	}
      }
      QColumnName=SDDS_FindColumn(&inputDataSet, FIND_NUMERIC_TYPE, "Q", "QualityFactor", NULL);
      if (QColumnName){
        switch(SDDS_CheckColumn(&inputDataSet, QColumnName, NULL, SDDS_DOUBLE, stderr)){
        case SDDS_CHECK_OKAY:
          break;
        default:
          printf("Something wrong with column %s in file %s. Wrong type?\n", QColumnName, filename);
          bombClinchor(NULL);
        }
      }
      deQFactorColumnName=SDDS_FindColumn(&inputDataSet, FIND_NUMERIC_TYPE,
                                          "deQFactor", "DeQFactor", NULL);
      if (deQFactorColumnName){
        switch(SDDS_CheckColumn(&inputDataSet, deQFactorColumnName, NULL, SDDS_DOUBLE, stderr)){
        case SDDS_CHECK_OKAY:
          break;
        default:
          printf("Something wrong with column %s in file %s. Wrong type?\n", deQFactorColumnName, filename);
          bombClinchor(NULL);
        }
      }
      CavitiesParameterName=SDDS_FindParameter(&inputDataSet, FIND_NUMERIC_TYPE,
                                               "NumberOfCavities", "NCavities",
                                               "Cavities", "cavities", NULL);
      if (CavitiesParameterName){
        switch(SDDS_CheckParameter(&inputDataSet, CavitiesParameterName, NULL, SDDS_LONG, stderr)){
        case SDDS_CHECK_OKAY:
          break;
        default:
          printf("Something wrong with parameter %s in file %s. Wrong type?\n", CavitiesParameterName, filename);
          bombClinchor(NULL);
        }
      }
      StaggeringStepColumnName=SDDS_FindColumn(&inputDataSet, FIND_NUMERIC_TYPE,
                                               "StaggeringStep", "staggeringStep",
                                               "DeltaFrequency", "deltaFrequency", "df", NULL);
      if (StaggeringStepColumnName){
        switch(SDDS_CheckColumn(&inputDataSet, StaggeringStepColumnName, "Hz", SDDS_DOUBLE, stderr)){
        case SDDS_CHECK_OKAY:
          break;
        default:
          printf("Something wrong with column %s in file %s. Wrong type?\n", StaggeringStepColumnName, filename);
          bombClinchor(NULL);
        }
      }
      ResonanceColumnName=SDDS_FindColumn(&inputDataSet, FIND_NUMERIC_TYPE,
                                          "ShiftToResonance", "shiftToResonance",
                                          "Resonance", "resonance", NULL);
      if (ResonanceColumnName){
        switch(SDDS_CheckColumn(&inputDataSet, ResonanceColumnName, NULL, SDDS_LONG, stderr)){
        case SDDS_CHECK_OKAY:
          break;
        default:
          printf("Something wrong with column %s in file %s. Wrong type?\n", ResonanceColumnName, filename);
          bombClinchor(NULL);
        }
      }
      FixedColumnName=SDDS_FindColumn(&inputDataSet, FIND_NUMERIC_TYPE,
                                      "Fixed", NULL);
      if (FixedColumnName){
        switch(SDDS_CheckColumn(&inputDataSet, FixedColumnName, NULL, SDDS_LONG, stderr)){
        case SDDS_CHECK_OKAY:
          break;
        default:
          printf("Something wrong with column %s in file %s. Wrong type?\n", FixedColumnName, filename);
          bombClinchor(NULL);
        }
      }

      /* if only one of ImpedanceColumnName RQColumnName QColumnName exists, then there is insufficient data */
      if (((!ImpedanceColumnName)&&(!RQColumnName))||
          ((!ImpedanceColumnName)&&(!QColumnName))||
          ((!RQColumnName)&&(!QColumnName))){
        printf("Error: Missing columns in file %s. Expect two of R, Q, or R/Q.\n", filename);
        bombClinchor(NULL);
      }
    }
    DataSetRows=SDDS_CountRowsOfInterest(&inputDataSet);
    /* allocate memory for resonator structure */
    switch(match_string(multipoleType, multipoleTypeList, MULTIPOLE_TYPES, 0)){
    case MONOPOLE:
      if (monopoleHOMs_struct.clearPreviousMonopoleHOMs || !monopoleHOMTypes) {
        offset = 0;
        monopoleHOMTypes = DataSetRows;
        monopoleHOMResonator = (HOM_RESONATOR *) tmalloc(sizeof(HOM_RESONATOR)*monopoleHOMTypes);
      }
      else {
        /* reallocate */
        offset = monopoleHOMTypes;
        monopoleHOMTypes += DataSetRows;
        monopoleHOMResonator = (HOM_RESONATOR *) trealloc(monopoleHOMResonator,
                                                          sizeof(HOM_RESONATOR)*monopoleHOMTypes);
	/* need to initialize the new memory :( */
	SDDS_ZeroMemory(&monopoleHOMResonator[offset], DataSetRows * sizeof(HOM_RESONATOR));	
      }
      if (verbose)
	printf("Filling monopole structure array with file %s, page %ld ...\n", monopoleHOMs_struct.filename, iDataSet);
      /* resonator will be the pointer to the next resonator to process in the next section */
      resonator=&monopoleHOMResonator[offset];
      break;   
    case DIPOLE:
      if (dipoleHOMs_struct.clearPreviousDipoleHOMs || !dipoleHOMTypes) {
        offset = 0;
        dipoleHOMTypes = DataSetRows;
        dipoleHOMResonator = (HOM_RESONATOR *) tmalloc(sizeof(HOM_RESONATOR)*dipoleHOMTypes);
      }
      else {
        /* reallocate */
        offset = dipoleHOMTypes;
        dipoleHOMTypes += DataSetRows;
        dipoleHOMResonator = (HOM_RESONATOR *) trealloc(dipoleHOMResonator,
                                                        sizeof(HOM_RESONATOR)*dipoleHOMTypes);
	/* need to initialize the new memory :( */
	SDDS_ZeroMemory(&dipoleHOMResonator[offset], DataSetRows *sizeof(HOM_RESONATOR)); 
      }
      if (verbose)
        printf("Filling dipole structure array with file %s ...\n", dipoleHOMs_struct.filename);
      /* resonator will be the pointer to the next resonator to process in the next section */
      resonator=&dipoleHOMResonator[offset];
      break;
    }

    if (CavitiesParameterName){
      if (!SDDS_GetParameter(&inputDataSet, CavitiesParameterName, &cavities))
        SDDS_PrintErrors(stderr, SDDS_EXIT_PrintErrors|SDDS_VERBOSE_PrintErrors);
    }
    else {
      cavities = 1;
    }

    for (row=0; row<DataSetRows; row++){
      resonator[row].cavities = cavities;
      if (!SDDS_GetValue(&inputDataSet, FrequencyColumnName, row, &resonator[row].unperturbedFrequency))
        SDDS_PrintErrors(stderr, SDDS_EXIT_PrintErrors|SDDS_VERBOSE_PrintErrors);
      resonator[row].unperturbedFrequency *= frequencyConversionFactor;
      if (ImpedanceColumnName){
        if (!SDDS_GetValue(&inputDataSet, ImpedanceColumnName, row, &resonator[row].shuntImpedance))
          SDDS_PrintErrors(stderr, SDDS_EXIT_PrintErrors|SDDS_VERBOSE_PrintErrors);
      }
      if (QColumnName){
        if (!SDDS_GetValue(&inputDataSet, QColumnName, row, &resonator[row].Q))
          SDDS_PrintErrors(stderr, SDDS_EXIT_PrintErrors|SDDS_VERBOSE_PrintErrors);
      }
      if (RQColumnName){
        if (!SDDS_GetValue(&inputDataSet, RQColumnName, row, &resonator[row].RoQ))
          SDDS_PrintErrors(stderr, SDDS_EXIT_PrintErrors|SDDS_VERBOSE_PrintErrors);
      }
      if (!ImpedanceColumnName)
        resonator[row].shuntImpedance=resonator[row].RoQ*resonator[row].Q;
      if (!QColumnName)
        resonator[row].Q=resonator[row].shuntImpedance/resonator[row].RoQ;
      if (!RQColumnName)
        resonator[row].RoQ=resonator[row].shuntImpedance/resonator[row].Q;
      if (deQFactorColumnName){
        if (!SDDS_GetValue(&inputDataSet, deQFactorColumnName, row, &resonator[row].deQFactor))
          SDDS_PrintErrors(stderr, SDDS_EXIT_PrintErrors|SDDS_VERBOSE_PrintErrors);
      }    
      else {
        resonator[row].deQFactor=1.0;
      }
      if (StaggeringStepColumnName){
        if (!SDDS_GetValue(&inputDataSet, StaggeringStepColumnName, row, &resonator[row].staggeringStep))
          SDDS_PrintErrors(stderr, SDDS_EXIT_PrintErrors|SDDS_VERBOSE_PrintErrors);
      }
      else {
        resonator[row].staggeringStep=0;
      }
      if (ResonanceColumnName){
        if (!SDDS_GetValue(&inputDataSet, ResonanceColumnName, row, &resonator[row].shiftToResonance))
          SDDS_PrintErrors(stderr, SDDS_EXIT_PrintErrors|SDDS_VERBOSE_PrintErrors);
      }
      else {
        resonator[row].shiftToResonance=0;
      }  
      if (FixedColumnName){
        if (!SDDS_GetValue(&inputDataSet, FixedColumnName, row, &resonator[row].fixed))
          SDDS_PrintErrors(stderr, SDDS_EXIT_PrintErrors|SDDS_VERBOSE_PrintErrors);
      }
      else {
        resonator[row].fixed=0;
      }  
    }
  }
}

void setupSymmetricBunchPattern(NAMELIST_TEXT *namelistText) {
#include "symmetricBunchPattern.h"
  long offset;
  long i;
  long bucketInterval;

  set_namelist_processing_flags(STICKY_NAMELIST_DEFAULTS);
  process_namelist(&symmetricBunchPattern, namelistText);
  print_namelist(stdout, &symmetricBunchPattern);
  if ( clearPreviousPatterns ) {
    if (bucketSelected)
      free(bucketSelected);
    if (bunchCurrent)
      free(bunchCurrent);
    totalFilledBuckets = 0;
    beamCurrentMA = 0;
  }
  offset = totalFilledBuckets;
  if ( clearPreviousPatterns ) {
    totalFilledBuckets = bunches;
    bucketSelected = (int32_t *)tmalloc(sizeof(int32_t)*bunches);
    bunchCurrent = (double *)tmalloc(sizeof(double)*bunches);
  }
  else {
    totalFilledBuckets += bunches;
    bucketSelected = (int32_t *)trealloc(bucketSelected, sizeof(int32_t)*totalFilledBuckets);
    bunchCurrent = (double *)trealloc(bunchCurrent, sizeof(double)*totalFilledBuckets);
  }
  if ( !currentPerBucketMA ) {
    if ( !totalCurrentMA )
      printf("Warning: buckets are to be assigned zero currents\n");
    else 
      currentPerBucketMA = totalCurrentMA/ bunches;
  }

  bucketInterval =  harmonicNumber/ bunches;
  
  if ( startBucket + (bunches-1)*bucketInterval >= harmonicNumber )
    printf("Warning: bucket number to be assigned exceeds harmonic number of %ld.\n The bucket number will be adjusted appropriately", harmonicNumber);
  for (i=0; i<bunches; i++) {
    bucketSelected[i+offset] = (startBucket + i * bucketInterval) % harmonicNumber;
    bunchCurrent[i+offset] = currentPerBucketMA;
    beamCurrentMA += currentPerBucketMA;
  }
  checkForBucketDuplication();
  if (bunchLengthTableFile) {
    determineBunchLengthByBunch();
    printBunchLengthsUsed();
  }
}

void setupBunchTrain(NAMELIST_TEXT *namelistText) {
#include "bunchTrain.h"
  long offset;
  long i;
  
  set_namelist_processing_flags(STICKY_NAMELIST_DEFAULTS);
  process_namelist(&bunchTrain, namelistText);
  print_namelist(stdout, &bunchTrain);
  if ( clearPreviousPatterns ) {
    if (bucketSelected)
      free(bucketSelected);
    if (bunchCurrent)
      free(bunchCurrent);
    totalFilledBuckets = 0;
    beamCurrentMA = 0;
  }
  offset = totalFilledBuckets;
  if ( clearPreviousPatterns ) {
    totalFilledBuckets = bunches;
    bucketSelected = (int32_t *)tmalloc(sizeof(int32_t)*bunches);
    bunchCurrent = (double *)tmalloc(sizeof(double)*bunches);
  }
  else {
    totalFilledBuckets += bunches;
    bucketSelected = (int32_t *)trealloc(bucketSelected, sizeof(int32_t)*totalFilledBuckets);
    bunchCurrent = (double *)trealloc(bunchCurrent, sizeof(double)*totalFilledBuckets);
  }
  if ( !currentPerBucketMA ) {
    if ( !totalCurrentMA )
      printf("Warning: buckets are to be assigned zero currents\n");
    else 
      currentPerBucketMA = totalCurrentMA/ bunches;
  }

  if (!bucketInterval) 
    printf("Warning: variable bucketInterval was set to 0. Now setting it to 1.\n");
  
  if ( startBucket + (bunches-1)*bucketInterval >= harmonicNumber )
    printf("Warning: bucket number to be assigned exceeds harmonic number of %ld.\n The bucket number will be adjusted appropriately", harmonicNumber);
  for (i=0; i<bunches; i++) {
    bucketSelected[i+offset] = (startBucket + i * bucketInterval) % harmonicNumber;
    bunchCurrent[i+offset] = currentPerBucketMA;
    beamCurrentMA += currentPerBucketMA;
  }
  checkForBucketDuplication();
  if (bunchLengthTableFile) {
    determineBunchLengthByBunch();
    printBunchLengthsUsed();
  }
}

void setupGeneralBunchPattern(NAMELIST_TEXT *namelistText) {
#include "generalBunchPattern.h"
  long offset;
  long bucket;
  double current;
  long  unequalNumEntries, warningPrinted;
  long bunches;
  long i;

  unequalNumEntries=warningPrinted=0;
  currentPerBucketMA = NULL;
  bucketSelection = NULL;
  totalCurrentMA = 0.0;
  
  set_namelist_processing_flags(STICKY_NAMELIST_DEFAULTS);
  process_namelist(&generalBunchPattern, namelistText);
  print_namelist(stdout, &generalBunchPattern);
  /* check presence of one current variable in the namelist */
  if ( bucketSelection == NULL )
    bomb("String variable bucketSelection not specified.\n", NULL);
  if ( currentPerBucketMA == NULL && !totalCurrentMA)
    printf("Warning : String variable currentPerBucketMA not specified and variable totalCurrentMA has not been set or has been set to zero.\n");
  /* parse strings containing bucket selection */
  if ( clearPreviousPatterns ) {
    if (bucketSelected)
      free(bucketSelected);
    if (bunchCurrent)
      free(bunchCurrent);
    totalFilledBuckets = 0;
    beamCurrentMA = 0;
  }        
  offset = totalFilledBuckets;
  while (get_long(&bucket, bucketSelection)) {
    if (!totalFilledBuckets) {
      /* first allocation */
      totalFilledBuckets++;
      bucketSelected = (int32_t *)tmalloc(sizeof(int32_t));
      bunchCurrent = (double *)tmalloc(sizeof(double));
    }
    else {
      totalFilledBuckets++;
      bucketSelected = (int32_t *)trealloc(bucketSelected, sizeof(int32_t)*totalFilledBuckets);
      bunchCurrent = (double *)trealloc(bunchCurrent, sizeof(double)*totalFilledBuckets);
    }
    if ( bucket>=harmonicNumber && !warningPrinted) {
      warningPrinted=1;
      printf("Warning: bucket number to be assigned exceeds harmonic number minus one.\n");
      printf("The bucket number will be adjusted appropriately.\n");
    }
    bucketSelected[totalFilledBuckets-1] = bucket%harmonicNumber;
  }
  bunches = totalFilledBuckets - offset;
  /* scan for current values */
  if (currentPerBucketMA == NULL) 
    for (i=offset; i<totalFilledBuckets; i++)
      bunchCurrent[i] = totalCurrentMA/ bunches;
  else {
    for (i=offset;  i<totalFilledBuckets; i++) {
      if ( get_double(&current, currentPerBucketMA) )
        bunchCurrent[i] = current;
      else {
        unequalNumEntries = 1;
        bunchCurrent[i] = 0.0;
      }
      beamCurrentMA += current;
    }
    if ( get_double(&current, currentPerBucketMA) || !unequalNumEntries)
      printf("Warning: Unequal number of entries in current_per_bucket_mA and bucket_selection.\nThe longest list is truncated.\n");
  }
  checkForBucketDuplication();
  if (bunchLengthTableFile) {
    determineBunchLengthByBunch();
    printBunchLengthsUsed();
  }
}
  
void setupBunchPatternFromFile(NAMELIST_TEXT *namelistText) {
#include "bunchPatternFromFile.h"
  long offset;
  int32_t *bucket;
  double *current;
  long  warningPrinted;
  long rows;
  long i;
  SDDS_TABLE bunchTable;

  filename=NULL;
  bucket=NULL;
  current=NULL;
  warningPrinted=0;

  set_namelist_processing_flags(STICKY_NAMELIST_DEFAULTS);
  process_namelist(&bunchPatternFromFile, namelistText);
  print_namelist(stdout, &bunchPatternFromFile);
  /* check presence of filename variable in the namelist */
  if ( filename == NULL )
    bomb("String variable filename not specified.\n", NULL);
  
  if ( clearPreviousPatterns ) {
    if (bucketSelected)
      free(bucketSelected);
    if (bunchCurrent)
      free(bunchCurrent);
    totalFilledBuckets = 0;
    beamCurrentMA = 0;
    bucketSelected=NULL;
    bunchCurrent=NULL;
  }        
  offset = totalFilledBuckets;
  if (!SDDS_InitializeInput(&bunchTable, filename)) {
    SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors);
    bombClinchor(NULL);
  }
  while ((rows=SDDS_ReadPage(&bunchTable))>0) {
    if (!(bucket=(int32_t*)SDDS_GetColumnInLong(&bunchTable,"bucket")))
      SDDS_PrintErrors(stderr, SDDS_EXIT_PrintErrors|SDDS_VERBOSE_PrintErrors);
    if (!(current=(double*)SDDS_GetColumnInDoubles(&bunchTable,"current")))
      SDDS_PrintErrors(stderr, SDDS_EXIT_PrintErrors|SDDS_VERBOSE_PrintErrors);
    if (!totalFilledBuckets) {
      /* first allocation */
      totalFilledBuckets=rows;
      bucketSelected = (int32_t *)tmalloc(sizeof(int32_t)*totalFilledBuckets);
      bunchCurrent = (double *)tmalloc(sizeof(double)*totalFilledBuckets);
    }
    else {
      totalFilledBuckets +=rows;
      bucketSelected = (int32_t *)trealloc(bucketSelected, sizeof(int32_t)*totalFilledBuckets);
      bunchCurrent = (double *)trealloc(bunchCurrent, sizeof(double)*totalFilledBuckets);
    }
    for (i=0;i<rows;i++) {
      if ( bucket[i]>=harmonicNumber && !warningPrinted) {
        warningPrinted=1;
        printf("Warning: bucket number to be assigned exceeds harmonic number minus one.\n");
        printf("The bucket number will be adjusted appropriately.\n");
      }
      bucketSelected[i+offset] = bucket[i]%harmonicNumber;
    }
    /* scan for current values */
    for (i=0;  i<rows; i++) {
      bunchCurrent[i+offset] = current[i];
      beamCurrentMA += current[i];
    }
    free(bucket);
    free(current);
    bucket=NULL;
    current=NULL;
    offset=totalFilledBuckets;
  }
  if (!SDDS_Terminate(&bunchTable))
    SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors);
  checkForBucketDuplication();
  if (bunchLengthTableFile) {
    determineBunchLengthByBunch();
    printBunchLengthsUsed();
  }
}

void setupRandomizeBunchCurrent(NAMELIST_TEXT *namelistText) {
  set_namelist_processing_flags(0);
  process_namelist(&randomizeBunchCurrent, namelistText);
  print_namelist(stdout, &randomizeBunchCurrent);

  if (!randomizeBunchCurrent_struct.absoluteSpreadMA && !randomizeBunchCurrent_struct.relativeSpread) {
    fprintf(stderr, "Variable absoluteSpreadMA and relativeSpread in namelist randomizeBunchCurrent are not specified or are set to zero in the input file.\nTo make the command meaningful, the value of spread is set to 10%%.\n");
    randomizeBunchCurrent_struct.relativeSpread = 0.10 ;

  }
  random_2(randomizeBunchCurrent_struct.seed);
  randomizeBunchCurrentFlag = 1;
}

void doRandomizeBunchCurrent() {
  /* Array bunchCurrent will be overwritten.
     the calling function should make a save of these values 
     for recovery */
  double sum1,sum2;
  long iBunch;
  
  sum1 = 0; /* total current */
  sum2 = 0; /* total current after randomization */
  for (iBunch=0; iBunch<totalFilledBuckets; iBunch++) {
    sum1 += bunchCurrent[iBunch];
    if (randomizeBunchCurrent_struct.absoluteSpreadMA) {
      bunchCurrent[iBunch] += randomizeBunchCurrent_struct.absoluteSpreadMA * (random_2(0)-0.5);
    }
    else {
      bunchCurrent[iBunch] *= 1 + randomizeBunchCurrent_struct.relativeSpread * (random_2(0)-0.5);
    }
    sum2 += bunchCurrent[iBunch];
  }
  for (iBunch=0; iBunch<totalFilledBuckets; iBunch++)
    bunchCurrent[iBunch] *= (sum1/sum2); /* normalize to original total current */
  
}


void checkForBucketDuplication(){
  long iBunch, jBunch, kBunch;
  for (iBunch=0; iBunch<totalFilledBuckets; iBunch++) {
    for (jBunch=iBunch+1; jBunch<totalFilledBuckets; jBunch++) {
      if (bucketSelected[iBunch] == bucketSelected[jBunch]) {
        /* remove jBunch entry and combine currents, and reduce number of bunches */
        printf("Warning bucket no. %d is found to be populated with more than one bunch. "\
               "The bunches will be combined into one.\n", bucketSelected[iBunch]);
        bunchCurrent[iBunch] += bunchCurrent[jBunch];
        totalFilledBuckets--;
        for ( kBunch=jBunch;  kBunch<totalFilledBuckets; kBunch++) {
          bunchCurrent[kBunch] = bunchCurrent[kBunch+1];
          bucketSelected[kBunch] = bucketSelected[kBunch+1];
        }
        /* for clarity, the last element is set to zero */
        bunchCurrent[totalFilledBuckets] = 0.0;
        bucketSelected[totalFilledBuckets] = 0;
      }
    }
  }
  /* sort the bunches possibly created by different commands */
  /* the bunches must be sorted so that the eigenvectors can be viewed
     with the order of elements corresponding to the same order of
     the bunches in the ring */
  for (iBunch=0; iBunch<totalFilledBuckets; iBunch++) {
    for (jBunch=iBunch+1; jBunch<totalFilledBuckets; jBunch++) {
      if (bucketSelected[iBunch] > bucketSelected[jBunch]) {
        SWAP_INT32(bucketSelected[iBunch], bucketSelected[jBunch]);
        SWAP_DOUBLE(bunchCurrent[iBunch], bunchCurrent[jBunch]);
      }
    }
  }
}

void determineBunchLengthByBunch(){
  long i, iBunch;
  unsigned long code;
  OUTRANGE_CONTROL aboveRange, belowRange;
  
  aboveRange.flags = belowRange.flags = OUTRANGE_SATURATE;
  bunchDurationByBunch = (double *) malloc( totalFilledBuckets * sizeof(double));
  if (!bunchCurrentData) {
    bunchCurrentData = (double *) malloc( bunchLengthValues * sizeof(double));
    for (i=0; i<bunchLengthValues; i++)
      /*     mA                   nC                     Hz        */ 
      bunchCurrentData[i] = bunchChargeData[i] * revolutionFrequency * 1e-6;
  }
    
  for (iBunch=0; iBunch<totalFilledBuckets; iBunch++) {
    bunchDurationByBunch[iBunch] = 1e-12 * interpolate( bunchLengthPicoSecondData, bunchCurrentData,
                                                        bunchLengthValues,
                                                        bunchCurrent[iBunch], &belowRange, &aboveRange,
                                                        2, &code, 0);
  }
}

void printBunchLengthsUsed () {
  SDDS_TABLE bunchLengths;

  /* bunchLengthFile was declared global variable */
  if (!SDDS_InitializeOutput(&bunchLengths, SDDS_BINARY, 1,
                             "Bunch lengths used for calculation of clinchor run",
                             "Bunch Lengths", bunchLengthFile))
    SDDS_PrintErrors(stderr, SDDS_EXIT_PrintErrors|SDDS_VERBOSE_PrintErrors);
  if (0>SDDS_DefineColumn(&bunchLengths, "Bucket", "Bucket", NULL,
                          NULL, NULL, SDDS_LONG, 0) ||
      0>SDDS_DefineColumn(&bunchLengths, "BunchCurrent", "I$bb$n", "mA",
                          NULL, NULL, SDDS_DOUBLE, 0) ||
      0>SDDS_DefineColumn(&bunchLengths, "Length", "$gs$n$bt$n$r", "ps",
                          NULL, NULL, SDDS_DOUBLE, 0) ||
      !SDDS_WriteLayout(&bunchLengths))
    SDDS_PrintErrors(stderr, SDDS_EXIT_PrintErrors|SDDS_VERBOSE_PrintErrors);
  if (!SDDS_StartPage(&bunchLengths, totalFilledBuckets))
    SDDS_PrintErrors(stderr, SDDS_EXIT_PrintErrors|SDDS_VERBOSE_PrintErrors);
  if (!SDDS_SetColumn(&bunchLengths, SDDS_SET_BY_NAME|SDDS_PASS_BY_REFERENCE,
                      bucketSelected, totalFilledBuckets, "Bucket") ||
      !SDDS_SetColumn(&bunchLengths, SDDS_SET_BY_NAME|SDDS_PASS_BY_REFERENCE,
                      bunchCurrent, totalFilledBuckets, "BunchCurrent") ||
      !SDDS_SetColumn(&bunchLengths, SDDS_SET_BY_NAME|SDDS_PASS_BY_REFERENCE,
                      bunchDurationByBunch, totalFilledBuckets, "Length") ||
      !SDDS_WritePage(&bunchLengths) || !SDDS_Terminate(&bunchLengths))
    SDDS_PrintErrors(stderr, SDDS_EXIT_PrintErrors|SDDS_VERBOSE_PrintErrors);
}


/* setup flags for sweeping HOM frequencies */
void setupSweepHOMFrequency(NAMELIST_TEXT *namelistText) {
  set_namelist_processing_flags(0);
  process_namelist(&sweepHOMFrequency, namelistText);
  print_namelist(stdout, &sweepHOMFrequency);
  if (!sweepHOMFrequency_struct.frequencyRange) {
    printf("Warning: Variable frequencyRange in namelist sweepHOMFrequency is unspecified or set to zero "
           "in the input file.\nVariable frequencyRange is set to %lg, the revolution frequency.\n", revolutionFrequency);
    sweepHOMFrequency_struct.frequencyRange = revolutionFrequency;
  }
  if (!sweepHOMFrequency_struct.points)
    bomb("Variable points in namelist sweep_HOMFrequency is set to zero in the input file.", NULL);
  if (sweepHOMFrequency_struct.filename == NULL)
    bomb("String filename in namelist sweep_HOMFrequency is not specified.", NULL);
  sweepFlag = 1;
  randomizeHOMFlag = 0;
}

/* setup flags for randomizing HOM frequencies */
void setupRandomizeHOMFrequencies(NAMELIST_TEXT *namelistText) {
  set_namelist_processing_flags(0);
  process_namelist(&randomizeHOMFrequencies, namelistText);
  print_namelist(stdout, &randomizeHOMFrequencies);
  if (randomizeHOMFrequencies_struct.CBMFrequencyFilename == NULL)
    bomb("String filename in namelist randomizeHOMFrequencies is not specified.", NULL);
  if (!randomizeHOMFrequencies_struct.spread) {
    printf("Variable spread in namelist randomize_HOM_frequencies is not specified or is set to zero in the input file.\n"
           "To make the command meaningful, the value of spread is set to the revolution frequency %lf Hz.\n", revolutionFrequency);
    randomizeHOMFrequencies_struct.spread = revolutionFrequency;
  }
  if (!randomizeHOMFrequencies_struct.uniform) {
    printf("Variable uniform in namelist randomizeHOMFrequencies is set to 1. (Only uniform distribution is available.)\n");
    randomizeHOMFrequencies_struct.uniform = 1;
  }
  random_1(randomizeHOMFrequencies_struct.seed);  
  randomizeHOMFlag = 1;
  sweepFlag = 0;
}

void calculateLongitudinalGrowthRate(NAMELIST_TEXT *namelistText) {
  COMPLEX f, f0, fShift, fDiff;
  long i, iBunch, iSample, iRes, iCav, iFundamental=0;
  double shift;
  long sweptResonator;
  HOM_RESONATOR *resonator;
  SDDS_DATASET CBMFrequencyDataSet, maxGrowthRateDataSet, randHOMFrequencyDataSet, CBMEigenvalues;
  long iHOM, HOMs;
  double totalCurrent;
  double *bunchCurrentSave=NULL;
  double power;

  dampingFactor = longDampingRate;
  f0.r = synchrotronAngFrequency;
  f0.i = -dampingFactor;

  set_namelist_processing_flags(0);
  process_namelist(&doLongitudinalMotion, namelistText);
  print_namelist(stdout, &doLongitudinalMotion);
  if (!doLongitudinalMotion_struct.normalModes)
    printf("Warning: Program can only do normal modes. Value of variable normalModes in namelist doLongitudinalMotion is set to 1.\n");
  if (doLongitudinalMotion_struct.doLaplace)
    printf("Warning: Program doesn't do laplace transforms. Value of variable doLaplace in namelist doLongitudinalMotion is set to 0.\n");
  if (doLongitudinalMotion_struct.eigenvectorFilename && (sweepFlag || randomizeHOMFlag))
    printf("Warning: Request for printing matrices and eigenvectors ignored when multiple calculations are performed.\n");
  
  /* allocate modified frequencies array, and assign staggered frequencies */
  if (!monopoleHOMTypes) {
    bombClinchor("No monopole HOMs defined.");
  }

  if (monopoleHOMs_struct.detuneFundamental) {
    double deviation, minDeviation = DBL_MAX;
    iFundamental = -1;
    /* find the fundamental */
    for (iRes=0; iRes<monopoleHOMTypes; iRes++) {
      if ((deviation=fabs(revolutionFrequency*harmonicNumber-monopoleHOMResonator[iRes].unperturbedFrequency))<minDeviation) {
        iFundamental = iRes;
        minDeviation = deviation;
      }
    }
    if (iFundamental>-1) {
      double beta, deltaF;
      
      resonator = &monopoleHOMResonator[iFundamental];
      printf("Fundamental mode is apparently at %e Hz\n", resonator->unperturbedFrequency);
      printf("Shifting to exact harmonic and detuning for optimum generator matching.\n");
      resonator->unperturbedFrequency = revolutionFrequency*harmonicNumber;
      beta = 1 + (beamCurrentMA/1e3)*(2*resonator->shuntImpedance)*sin(synchronousPhase)/(1e6*rfVoltageMV);
      printf("Optimal beta is %e\n", beta);
      deltaF = -(beta-1)*resonator->unperturbedFrequency/(2*resonator->Q*tan(synchronousPhase));
      printf("Detuning cavity by %e Hz\n", deltaF);
      resonator->unperturbedFrequency += deltaF;
    }
  }
  for (iRes=0; iRes<monopoleHOMTypes; iRes++) {
    resonator = &monopoleHOMResonator[iRes];
    if (resonator->cavities) {
      if (!resonator->frequency) { 
        resonator->frequency = (double *) tmalloc(sizeof(double)*resonator->cavities);
        resonator->power = (double *) tmalloc(sizeof(double)*resonator->cavities);
        resonator->powerPlus = (double *) tmalloc(sizeof(double)*resonator->cavities);
        resonator->powerMinus = (double *) tmalloc(sizeof(double)*resonator->cavities);
      }
      else {
	resonator->frequency = (double *) trealloc(resonator->frequency, sizeof(double)*resonator->cavities);
	resonator->power = (double *) trealloc(resonator->power, sizeof(double)*resonator->cavities);
	resonator->powerPlus = (double *) trealloc(resonator->powerPlus, sizeof(double)*resonator->cavities);
	resonator->powerMinus = (double *) trealloc(resonator->powerMinus, sizeof(double)*resonator->cavities);
      }
      for (iCav=0; iCav<resonator->cavities; iCav++) {
	/* request for fixed frequency overrules shift to resonance */
        if (resonator->fixed) {
          resonator->frequency[iCav] = resonator->unperturbedFrequency;
        } else {
          resonator->frequency[iCav] = resonator->unperturbedFrequency + iCav * resonator->staggeringStep * monopoleHOMs_struct.staggeringStepMultiplier;
          /* apply shift to nearest revolution harmonic upper sideband, if requested. */
          if ( resonator->shiftToResonance )
            resonator->frequency[iCav] =
              ((long)(((resonator->frequency[iCav]-synchrotronFrequency)/revolutionFrequency) + 0.5)) *
              revolutionFrequency + synchrotronFrequency;
        }
      }
    }
  }
  
  /* setup allocation and pointers for printing eigenvalues and eigenvectors */
  if (!lapack_eigenvalues) {
    lapack_eigenvalues = (lapack_complex_double*) tmalloc(sizeof(lapack_complex_double)*totalFilledBuckets);
  }
  else {
    lapack_eigenvalues = (lapack_complex_double*) trealloc(lapack_eigenvalues, sizeof(lapack_complex_double)*totalFilledBuckets);
  }

  totalCurrent = 0;
  for (i=0;i<totalFilledBuckets;i++) {
    totalCurrent += bunchCurrent[i];
  }
  
  /* if sweep */
  if ( sweepFlag ) {
    sweptResonator = sweepHOMFrequency_struct.resonatorIndex;
    resonator = &monopoleHOMResonator[sweptResonator];
    if (!SDDS_InitializeOutput(&CBMFrequencyDataSet, SDDS_BINARY, 1,
                               "Coherent longitudinal oscillation frequency as a function of HOM Frequency",
                               "Mode frequencies", sweepHOMFrequency_struct.filename)||
        (0>SDDS_DefineParameter(&CBMFrequencyDataSet, "TotalCurrent", "I$btot$n", "mA",
                                NULL, NULL, SDDS_DOUBLE, 0)) ||
        (0>SDDS_DefineParameter(&CBMFrequencyDataSet, "Bunches", "N$bb$n", NULL,
                                NULL, NULL, SDDS_LONG, 0)) ||
        (0>SDDS_DefineColumn(&CBMFrequencyDataSet, "HOMFrequency", "HOM Frequency", "Hz", NULL, NULL, SDDS_DOUBLE, 0))||
        (0>SDDS_DefineColumn(&CBMFrequencyDataSet, "RealCBMFrequency", "Re{$gw$r}", "rad/sec",
                             "Real part of mode frequency", NULL, SDDS_DOUBLE, 0))||
        (0>SDDS_DefineColumn(&CBMFrequencyDataSet, "ImagCBMFrequency", "Im{$gw$r}", "rad/sec",
                             "Growth rate", NULL, SDDS_DOUBLE, 0))||
        (0>SDDS_DefineColumn(&CBMFrequencyDataSet, "RealDeltaCBMFrequency", "Re{$gDw$r}", "rad/sec",
                             "Real part of mode frequency shift", NULL, SDDS_DOUBLE, 0))||
        (0>SDDS_DefineColumn(&CBMFrequencyDataSet, "ImagDeltaCBMFrequency", "Im{$gDw$r}", "rad/sec",
                             "Imaginary part of mode frequency shift", NULL, SDDS_DOUBLE, 0))||
        (0>SDDS_DefineColumn(&CBMFrequencyDataSet, "MagDeltaCBMFrequency", "$sb$e$gDw$r$sb$e", "rad/sec",
                             "Magnitude of frequency shift", NULL, SDDS_DOUBLE, 0))||        
	(0>SDDS_DefineColumn(&CBMFrequencyDataSet, "Power", "P", "W",
                             "Instantaneous Power at bunch 0", NULL, SDDS_DOUBLE, 0))||
        !SDDS_WriteLayout(&CBMFrequencyDataSet)||
        !SDDS_StartPage(&CBMFrequencyDataSet, sweepHOMFrequency_struct.points))
      SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
    if( !SDDS_SetParameters(&CBMFrequencyDataSet, SDDS_BY_NAME|SDDS_PASS_BY_VALUE,
                            "TotalCurrent", totalCurrent,
                            "Bunches", totalFilledBuckets, NULL))
      SDDS_PrintErrors(stderr, SDDS_EXIT_PrintErrors|SDDS_VERBOSE_PrintErrors);
    
    for (i=0; i<sweepHOMFrequency_struct.points; i++) {
      shift = 2*i*sweepHOMFrequency_struct.frequencyRange/ (sweepHOMFrequency_struct.points - 1)
        - sweepHOMFrequency_struct.frequencyRange;
      /* all cavities of the same hom type are shifted the same way */
      for (iCav=0; iCav<resonator->cavities; iCav++) {
        resonator->frequency[iCav] = resonator->unperturbedFrequency + shift;
      }

      /* randomize bunch current if requested */
      if (randomizeBunchCurrentFlag) {
        bunchCurrentSave = (double *) malloc( totalFilledBuckets * sizeof(double));
        for (iBunch=0; iBunch<totalFilledBuckets; iBunch++)
          bunchCurrentSave[iBunch] = bunchCurrent[iBunch];
        doRandomizeBunchCurrent();
      }
      power = powerOfResonators();
      f = frequencyOfMaxLongGrowthRate();
      
      if (randomizeBunchCurrentFlag) {
        for (iBunch=0; iBunch<totalFilledBuckets; iBunch++)
          bunchCurrent[iBunch] = bunchCurrentSave[iBunch];
        free(bunchCurrentSave);      
      }
      
      fShift = csub(f, f0);
      if (!SDDS_SetRowValues(&CBMFrequencyDataSet, SDDS_BY_NAME|SDDS_PASS_BY_VALUE, i,
                             "HOMFrequency", resonator->frequency[0],
                             "RealCBMFrequency", f.r,
                             "ImagCBMFrequency", f.i,
                             "RealDeltaCBMFrequency", fShift.r,
                             "ImagDeltaCBMFrequency", fShift.i,
                             "MagDeltaCBMFrequency", cmod(fShift), 
			     "Power", power, NULL))
        SDDS_PrintErrors(stderr, SDDS_EXIT_PrintErrors|SDDS_VERBOSE_PrintErrors);
    }
    if (!SDDS_WritePage(&CBMFrequencyDataSet)||!SDDS_Terminate(&CBMFrequencyDataSet))
      SDDS_PrintErrors(stderr, SDDS_EXIT_PrintErrors|SDDS_VERBOSE_PrintErrors);
  }
  else if (randomizeHOMFlag) {
    /* if random */
    if (!SDDS_InitializeOutput(&maxGrowthRateDataSet, SDDS_BINARY, 1,
                               "Maximum longitudinal CBM Growth Rates for randomized HOM frequencies", "CBM Growth Rates",
                               randomizeHOMFrequencies_struct.CBMFrequencyFilename)||
        (0>SDDS_DefineParameter(&maxGrowthRateDataSet, "TotalCurrent", "I$btot$n", "mA",
                                NULL, NULL, SDDS_DOUBLE, 0)) ||
        (0>SDDS_DefineParameter(&maxGrowthRateDataSet, "Bunches", "N$bb$n", NULL,
                                NULL, NULL, SDDS_LONG, 0)) ||
        (0>SDDS_DefineColumn(&maxGrowthRateDataSet, "MaxGrowthRate", "(1/$gt$r)$bmax$n", "1/sec",
                             "Maximum growth rate among all possible CBMs", NULL, SDDS_DOUBLE, 0))||
        (0>SDDS_DefineColumn(&maxGrowthRateDataSet, "Power", "P", "W",
                             "Instantaneous Power at bunch 0", NULL, SDDS_DOUBLE, 0))||
        !SDDS_WriteLayout(&maxGrowthRateDataSet)||
        !SDDS_StartPage(&maxGrowthRateDataSet, randomizeHOMFrequencies_struct.samples))
      SDDS_PrintErrors(stderr, SDDS_EXIT_PrintErrors|SDDS_VERBOSE_PrintErrors);
    if( !SDDS_SetParameters(&maxGrowthRateDataSet, SDDS_BY_NAME|SDDS_PASS_BY_VALUE,
                            "TotalCurrent", totalCurrent,
                            "Bunches", totalFilledBuckets, NULL))
      SDDS_PrintErrors(stderr, SDDS_EXIT_PrintErrors|SDDS_VERBOSE_PrintErrors);
    
    if (randomizeHOMFrequencies_struct.HOMFrequencyFilename) {
      iHOM=0;
      for (iRes=0; iRes<monopoleHOMTypes; iRes++) {
        resonator = &monopoleHOMResonator[iRes];
        if (resonator->cavities) 
          for (iCav=0; iCav<resonator->cavities; iCav++) {
            iHOM++;
          }
      }
      HOMs=iHOM;
      if (!SDDS_InitializeOutput(&randHOMFrequencyDataSet, SDDS_BINARY, 1,
                                 "Randomized HOM Frequencies", "Randomized HOM Frequencies",
                                 randomizeHOMFrequencies_struct.HOMFrequencyFilename)||
          (0>SDDS_DefineParameter(&randHOMFrequencyDataSet, "TotalCurrent", "I$btot$n", "mA",
                                  NULL, NULL, SDDS_DOUBLE, 0)) ||
          (0>SDDS_DefineParameter(&randHOMFrequencyDataSet, "Bunches", "N$bb$n", NULL,
                                  NULL, NULL, SDDS_LONG, 0)) ||
          (0>SDDS_DefineParameter(&randHOMFrequencyDataSet, "Sample", "Sample", NULL,
                                  NULL, NULL, SDDS_LONG, 0))||
          (0>SDDS_DefineColumn(&randHOMFrequencyDataSet, "HOMTypeIndex", "HOM type index", NULL,
                               NULL, NULL, SDDS_LONG, 0)) ||
          (0>SDDS_DefineColumn(&randHOMFrequencyDataSet, "CavityIndex", "Cavity index", NULL,
                               NULL, NULL, SDDS_LONG, 0)) ||
          (0>SDDS_DefineColumn(&randHOMFrequencyDataSet, "HOMFrequency", "HOM Frequency", "Hz",
                               NULL, NULL, SDDS_DOUBLE, 0)) ||
          !SDDS_WriteLayout(&randHOMFrequencyDataSet)||
          !SDDS_StartPage(&randHOMFrequencyDataSet, HOMs))
        SDDS_PrintErrors(stderr, SDDS_EXIT_PrintErrors|SDDS_VERBOSE_PrintErrors);
      if( !SDDS_SetParameters(&randHOMFrequencyDataSet, SDDS_BY_NAME|SDDS_PASS_BY_VALUE,
                              "TotalCurrent", totalCurrent,
                              "Bunches", totalFilledBuckets, NULL))
        SDDS_PrintErrors(stderr, SDDS_EXIT_PrintErrors|SDDS_VERBOSE_PrintErrors);
    }
    for (iSample=0; iSample<randomizeHOMFrequencies_struct.samples; iSample++) {
      iHOM=0;
      /* calculating all resonator frequencies */
      for (iRes=0; iRes<monopoleHOMTypes; iRes++) {
        resonator = &monopoleHOMResonator[iRes];
        if (resonator->cavities) {
          for (iCav=0; iCav<resonator->cavities; iCav++) {
            if (iRes==iFundamental && monopoleHOMs_struct.keepFundamentalFixed) {
              resonator->frequency[iCav] = resonator->unperturbedFrequency;
            } else if (resonator->fixed) {
              resonator->frequency[iCav] = resonator->unperturbedFrequency;
            } else {
              resonator->frequency[iCav] = resonator->unperturbedFrequency + iCav * resonator->staggeringStep * monopoleHOMs_struct.staggeringStepMultiplier;
              /*  the shift to nearest revolution harmonic upper sideband is applied before randomization */
              if ( resonator->shiftToResonance )
                resonator->frequency[iCav] =
                  ((long)(((resonator->frequency[iCav]-synchrotronFrequency)/revolutionFrequency) + 0.5)) *
                  revolutionFrequency + synchrotronFrequency;
              resonator->frequency[iCav] += 
                randomizeHOMFrequencies_struct.spread * (random_1(0)-0.5);
            }
            if (randomizeHOMFrequencies_struct.HOMFrequencyFilename) {
              if (!SDDS_SetRowValues(&randHOMFrequencyDataSet, SDDS_BY_NAME|SDDS_PASS_BY_VALUE, iHOM,
                                     "HOMTypeIndex", iRes,
                                     "CavityIndex", iCav,
                                     "HOMFrequency", resonator->frequency[iCav], NULL))
                SDDS_PrintErrors(stderr, SDDS_EXIT_PrintErrors|SDDS_VERBOSE_PrintErrors);
              iHOM++;
            }
          }
        }
      }
      if (randomizeHOMFrequencies_struct.HOMFrequencyFilename)
        if (!SDDS_SetParameters(&randHOMFrequencyDataSet, SDDS_BY_NAME|SDDS_PASS_BY_VALUE,
                                "Sample", iSample, NULL) ||
            !SDDS_WritePage(&randHOMFrequencyDataSet))
          SDDS_PrintErrors(stderr, SDDS_EXIT_PrintErrors|SDDS_VERBOSE_PrintErrors);

      /* randomize bunch current if requested */
      if (randomizeBunchCurrentFlag) {
        bunchCurrentSave = (double *) malloc( totalFilledBuckets * sizeof(double));
        for (iBunch=0; iBunch<totalFilledBuckets; iBunch++)
          bunchCurrentSave[iBunch] = bunchCurrent[iBunch];
        doRandomizeBunchCurrent();
      }

      power = powerOfResonators();
      f = frequencyOfMaxLongGrowthRate();
      
      if (randomizeBunchCurrentFlag) {
        for (iBunch=0; iBunch<totalFilledBuckets; iBunch++)
          bunchCurrent[iBunch] = bunchCurrentSave[iBunch];
        free(bunchCurrentSave);      
      }
      
      if (!SDDS_SetRowValues(&maxGrowthRateDataSet, SDDS_BY_NAME|SDDS_PASS_BY_VALUE, iSample,
                             "MaxGrowthRate", f.i, 
			     "Power", power,  NULL))
        SDDS_PrintErrors(stderr, SDDS_EXIT_PrintErrors|SDDS_VERBOSE_PrintErrors);
    }
    if (!SDDS_WritePage(&maxGrowthRateDataSet) || !SDDS_Terminate(&maxGrowthRateDataSet))
      SDDS_PrintErrors(stderr, SDDS_EXIT_PrintErrors|SDDS_VERBOSE_PrintErrors);
    if (randomizeHOMFrequencies_struct.HOMFrequencyFilename)
      if (!SDDS_Terminate(&randHOMFrequencyDataSet))
        SDDS_PrintErrors(stderr, SDDS_EXIT_PrintErrors|SDDS_VERBOSE_PrintErrors);
  }
  else {
    /*  one time straight calculation.  */
    /* randomize bunch current if requested */
    if (randomizeBunchCurrentFlag) {
      bunchCurrentSave = (double *) malloc( totalFilledBuckets * sizeof(double));
      for (iBunch=0; iBunch<totalFilledBuckets; iBunch++)
        bunchCurrentSave[iBunch] = bunchCurrent[iBunch];
      doRandomizeBunchCurrent();
    }
    power = powerOfResonators();
    f = frequencyOfMaxLongGrowthRate();
    
    if (randomizeBunchCurrentFlag) {
      for (iBunch=0; iBunch<totalFilledBuckets; iBunch++)
        bunchCurrent[iBunch] = bunchCurrentSave[iBunch];
      free(bunchCurrentSave);      
    }
    
    printf("Maximum growth rate is %lf 1/sec.\n", f.i);
    printf("Instantaneous HOM power at bunch 0 is %lf W.\n", power);
    /* list the frequencies of all modes */
    if (!doLongitudinalMotion_struct.CBMFrequencyFilename) {
      bomb("Filename CBMFrequencyFilename not specified in namelist doLongitudinalMotion.", NULL);
    }   

    if (!SDDS_InitializeOutput(&CBMEigenvalues, SDDS_BINARY, 1,
                               "CBM eigenvalues", "CBM eigenvalues",
                               doLongitudinalMotion_struct.CBMFrequencyFilename)||
        (0>SDDS_DefineParameter(&CBMEigenvalues, "TotalCurrent", "I$btot$n", "mA",
                                NULL, NULL, SDDS_DOUBLE, 0)) ||
        (0>SDDS_DefineParameter(&CBMEigenvalues, "Bunches", "N$bb$n", NULL,
                                NULL, NULL, SDDS_LONG, 0)) ||
        (0>SDDS_DefineColumn(&CBMEigenvalues, "CBMIndex", "CBM Index", NULL, "CBM index", NULL, SDDS_LONG, 0))||
        (0>SDDS_DefineColumn(&CBMEigenvalues, "RealCBMFrequency", "Re{$gw$r}", NULL,
                             "Real part of CBM frequency", NULL, SDDS_DOUBLE, 0))||
        (0>SDDS_DefineColumn(&CBMEigenvalues, "ImagCBMFrequency", "Im{$gw$r}", NULL,
                             "Imag part of CBM frequency", NULL, SDDS_DOUBLE, 0))||
        (0>SDDS_DefineColumn(&CBMEigenvalues, "RealDeltaCBMFrequency", "Re{$gDw$r}", "rad/sec",
                             "Real part of CBM frequency shift", NULL, SDDS_DOUBLE, 0))||
        (0>SDDS_DefineColumn(&CBMEigenvalues, "ImagDeltaCBMFrequency", "Im{$gDw$r}", "rad/sec",
                             "Imaginary part of CBM frequency shift", NULL, SDDS_DOUBLE, 0))||
        (0>SDDS_DefineColumn(&CBMEigenvalues, "MagDeltaCBMFrequency", "$sb$e$gDw$r$sb$e", "rad/sec",
                             "Magnitude of frequency shift", NULL, SDDS_DOUBLE, 0))||
        !SDDS_WriteLayout(&CBMEigenvalues)||
        !SDDS_StartPage(&CBMEigenvalues, totalFilledBuckets))
      SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
    
    if( !SDDS_SetParameters(&CBMEigenvalues, SDDS_BY_NAME|SDDS_PASS_BY_VALUE,
                            "TotalCurrent", totalCurrent,
                            "Bunches", totalFilledBuckets, NULL))
      SDDS_PrintErrors(stderr, SDDS_EXIT_PrintErrors|SDDS_VERBOSE_PrintErrors);

    for (i=0; i<totalFilledBuckets; i++) {
      //fDiff=csub(eigenvalues[i], f0);
      fDiff.r = lapack_eigenvalues[i].real - f0.r;
      fDiff.i = lapack_eigenvalues[i].imag - f0.i;
      if (!SDDS_SetRowValues(&CBMEigenvalues, SDDS_BY_NAME|SDDS_PASS_BY_VALUE, i,
                             "CBMIndex", i,
                             "RealCBMFrequency", lapack_eigenvalues[i].real,
                             "ImagCBMFrequency", lapack_eigenvalues[i].imag,
                             "RealDeltaCBMFrequency", fDiff.r,
                             "ImagDeltaCBMFrequency", fDiff.i,
                             "MagDeltaCBMFrequency", cmod(fDiff), NULL))
        SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
    }
    if (!SDDS_WritePage(&CBMEigenvalues)||!SDDS_Terminate(&CBMEigenvalues))
      SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
  }
  /* reset randomize flag */
  randomizeHOMFlag = 0;
  /* reset sweep flag */
  sweepFlag = 0;
}


void calculateTransverseGrowthRate(NAMELIST_TEXT *namelistText) {
  COMPLEX f, f0, fShift, fDiff;
  long i, iBunch, iSample, iRes, iCav;
  double shift;
  long sweptResonator;
  HOM_RESONATOR *resonator;
  SDDS_DATASET CBMFrequencyDataSet, maxGrowthRateDataSet, randHOMFrequencyDataSet, CBMEigenvalues;
  long iHOM, HOMs;
  double tune=0, betatronFrequency;
  double totalCurrent;
  double *bunchCurrentSave=NULL;

  set_namelist_processing_flags(0);
  process_namelist(&doTransverseMotion, namelistText);
  print_namelist(stdout, &doTransverseMotion);

  if (!doTransverseMotion_struct.normalModes)
    printf("Warning: Program can only do normal modes. "
           "Value of variable normalModes in namelist doTransverseMotion is set to 1.\n");
  if (doTransverseMotion_struct.doLaplace)
    printf("Warning: Program doesn't do laplace transforms. "
           "Value of variable doLaplace in namelist doTransverseMotion is set to 0.\n");
  if (doTransverseMotion_struct.eigenvectorFilename && (sweepFlag || randomizeHOMFlag))
    printf("Warning: Request for printing matrices and eigenvectors ignored when multiple calculations are performed.\n");
  if (!dipoleHOMTypes) {
    bombClinchor("No dipole HOMs defined.");
  }
  /* check tune values */
  if ( !doTransverseMotion_struct.horizontalDirection && !doTransverseMotion_struct.verticalDirection)
    bomb("Neither horizontal_direction nor vertical_direction was specified in namelist doTransverseMotion"
         "\n or they were both set to zero in namelist.", NULL);
  if ( doTransverseMotion_struct.horizontalDirection ) {
    if (!horizontalTune)
      bomb("Value of horizontalTune in namelist ringParameters was not specified or was set to 0.0", NULL);
    if (!betaxAtRFCavities)
      bomb("Value of betaxAtRFCavities in namelist ringParameters was not specified or was set to 0.0", NULL);
    betaAtRFCavs = betaxAtRFCavities;
    tune = horizontalTune;
    dampingFactor = horTransDampingRate;
  }
  if (doTransverseMotion_struct.verticalDirection) {
    if (!verticalTune)
      bomb("Value of verticalTune in namelist ringParameters was not specified or was set to 0.0", NULL);
    if (!betayAtRFCavities)
      bomb("Value of betayAtRFCavities in namelist ringParameters was not specified or was set to 0.0", NULL);
    betaAtRFCavs = betayAtRFCavities;         
    tune = verticalTune;
    dampingFactor = vertTransDampingRate;
  }
  betatronFrequency = tune * revolutionFrequency;
  betatronAngFrequency = 2 * PI * betatronFrequency;
  f0.r = betatronAngFrequency;
  f0.i = -dampingFactor;

  /* allocate modified frequencies array, and assign staggered frequencies */
  for (iRes=0; iRes<dipoleHOMTypes; iRes++) {
    resonator = &dipoleHOMResonator[iRes];
    if (resonator->cavities) {
      if (!resonator->frequency) 
        resonator->frequency = (double *) tmalloc(sizeof(double)*resonator->cavities);
      else 
        resonator->frequency = 
          (double *) trealloc(resonator->frequency, sizeof(double)*resonator->cavities);
      for (iCav=0; iCav<resonator->cavities; iCav++) {
        if (resonator->fixed) {
          resonator->frequency[iCav] = resonator->unperturbedFrequency;
        } else {
          resonator->frequency[iCav] = resonator->unperturbedFrequency + iCav * resonator->staggeringStep * dipoleHOMs_struct.staggeringStepMultiplier;
          /* apply shift to nearest revolution harmonic lower betatron sideband, if necesary */
          if ( resonator->shiftToResonance )
            resonator->frequency[iCav] =
              ((long)(((resonator->frequency[iCav]+betatronFrequency)/revolutionFrequency) + 0.5)) *
              revolutionFrequency - betatronFrequency;
        }
      }
    }
  }

  /* setup allocation and pointers for printing eigenvalues and eigenvectors */
  if (!lapack_eigenvalues) {
    lapack_eigenvalues = (lapack_complex_double*) tmalloc(sizeof(lapack_complex_double)*totalFilledBuckets);
  } else {
    lapack_eigenvalues = (lapack_complex_double*) trealloc(lapack_eigenvalues, sizeof(lapack_complex_double)*totalFilledBuckets);
  }

  totalCurrent = 0;
  for (i=0;i<totalFilledBuckets;i++) {
    totalCurrent += bunchCurrent[i];
  }
  /* Introduce randomization of bunch current to work
     with both sweep and randomize HOM commands */

  /* if sweep */
  if ( sweepFlag ) {
    sweptResonator = sweepHOMFrequency_struct.resonatorIndex;
    resonator = &dipoleHOMResonator[sweptResonator];
    if (!SDDS_InitializeOutput(&CBMFrequencyDataSet, SDDS_BINARY, 1,
                               "Coherent transverse oscillation frequency as a function of HOM Frequency",
                               "CBM frequencies", sweepHOMFrequency_struct.filename)||
        (0>SDDS_DefineParameter(&CBMFrequencyDataSet, "TotalCurrent", "I$btot$n", "mA",
                                NULL, NULL, SDDS_DOUBLE, 0)) ||
        (0>SDDS_DefineParameter(&CBMFrequencyDataSet, "Bunches", "N$bb$n", NULL,
                                NULL, NULL, SDDS_LONG, 0)) ||
        (0>SDDS_DefineColumn(&CBMFrequencyDataSet, "HOMFrequency", "HOM Frequency", "Hz", NULL, NULL, SDDS_DOUBLE, 0))||
        (0>SDDS_DefineColumn(&CBMFrequencyDataSet, "RealCBMFrequency", "Re{$gw$r}", "rad/sec",
                             "Real part of CBM frequency", NULL, SDDS_DOUBLE, 0))||
        (0>SDDS_DefineColumn(&CBMFrequencyDataSet, "ImagCBMFrequency", "Im{$gw$r}", "rad/sec",
                             "Growth rate", NULL, SDDS_DOUBLE, 0))||
        (0>SDDS_DefineColumn(&CBMFrequencyDataSet, "RealDeltaCBMFrequency", "Re{$gDw$r}", "rad/sec",
                             "Real part of CBM frequency shift", NULL, SDDS_DOUBLE, 0))||
        (0>SDDS_DefineColumn(&CBMFrequencyDataSet, "ImagDeltaCBMFrequency", "Im{$gDw$r}", "rad/sec",
                             "Imaginary part of CBM frequency shift", NULL, SDDS_DOUBLE, 0))||
        (0>SDDS_DefineColumn(&CBMFrequencyDataSet, "MagDeltaCBMFrequency", "$sb$e$gDw$r$sb$e", "rad/sec",
                             "Magnitude of frequency shift", NULL, SDDS_DOUBLE, 0))||
        !SDDS_WriteLayout(&CBMFrequencyDataSet)||
        !SDDS_StartPage(&CBMFrequencyDataSet, sweepHOMFrequency_struct.points))
      SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
    if( !SDDS_SetParameters(&CBMFrequencyDataSet, SDDS_BY_NAME|SDDS_PASS_BY_VALUE,
                            "TotalCurrent", totalCurrent,
                            "Bunches", totalFilledBuckets, NULL))
      SDDS_PrintErrors(stderr, SDDS_EXIT_PrintErrors|SDDS_VERBOSE_PrintErrors);
    
    for (i=0; i<sweepHOMFrequency_struct.points; i++) {
      shift = 2 * i * sweepHOMFrequency_struct.frequencyRange/ (sweepHOMFrequency_struct.points - 1)
        - sweepHOMFrequency_struct.frequencyRange;
      /* all cavities of the same hom type are shifted the same way */
      for (iCav=0; iCav<resonator->cavities; iCav++) {
        resonator->frequency[iCav] = resonator->unperturbedFrequency + shift;
      }
      /* randomize bunch current if requested */
      if (randomizeBunchCurrentFlag) {
        bunchCurrentSave = (double *) malloc( totalFilledBuckets * sizeof(double));
        for (iBunch=0; iBunch<totalFilledBuckets; iBunch++)
          bunchCurrentSave[iBunch] = bunchCurrent[iBunch];
        doRandomizeBunchCurrent();
      }
      
      f = frequencyOfMaxTransGrowthRate();
      
      if (randomizeBunchCurrentFlag) {
        for (iBunch=0; iBunch<totalFilledBuckets; iBunch++)
          bunchCurrent[iBunch] = bunchCurrentSave[iBunch];
        free(bunchCurrentSave);      
      }
      
      fShift=csub(f, f0);
      if (!SDDS_SetRowValues(&CBMFrequencyDataSet, SDDS_BY_NAME|SDDS_PASS_BY_VALUE, i,
                             "HOMFrequency", resonator->frequency[0],
                             "RealCBMFrequency", f.r,
                             "ImagCBMFrequency", f.i,
                             "RealDeltaCBMFrequency", fShift.r,
                             "ImagDeltaCBMFrequency", fShift.i,
                             "MagDeltaCBMFrequency", cmod(fShift), NULL))
        SDDS_PrintErrors(stderr, SDDS_EXIT_PrintErrors|SDDS_VERBOSE_PrintErrors);
    }
    if (!SDDS_WritePage(&CBMFrequencyDataSet)||!SDDS_Terminate(&CBMFrequencyDataSet))
      SDDS_PrintErrors(stderr, SDDS_EXIT_PrintErrors|SDDS_VERBOSE_PrintErrors);
  }
  else if (randomizeHOMFlag) {
    /* if random */
    if (!SDDS_InitializeOutput(&maxGrowthRateDataSet, SDDS_BINARY, 1,
                               "Maximum transverse CBM Growth Rates for randomized HOM frequencies", "CBM Growth Rates",
                               randomizeHOMFrequencies_struct.CBMFrequencyFilename)||
        (0>SDDS_DefineParameter(&maxGrowthRateDataSet, "TotalCurrent", "I$btot$n", "mA",
                                NULL, NULL, SDDS_DOUBLE, 0)) ||
        (0>SDDS_DefineParameter(&maxGrowthRateDataSet, "Bunches", "N$bb$n", NULL,
                                NULL, NULL, SDDS_LONG, 0)) ||
        (0>SDDS_DefineColumn(&maxGrowthRateDataSet, "MaxGrowthRate", "(1/$gt$r)$bmax$n", "1/sec",
                             "Maximum growth rate among all possible CBMs", NULL, SDDS_DOUBLE, 0))||
        !SDDS_WriteLayout(&maxGrowthRateDataSet)||
        !SDDS_StartPage(&maxGrowthRateDataSet, randomizeHOMFrequencies_struct.samples))
      SDDS_PrintErrors(stderr, SDDS_EXIT_PrintErrors|SDDS_VERBOSE_PrintErrors);
    if( !SDDS_SetParameters(&maxGrowthRateDataSet, SDDS_BY_NAME|SDDS_PASS_BY_VALUE,
                            "TotalCurrent", totalCurrent,
                            "Bunches", totalFilledBuckets, NULL))
      SDDS_PrintErrors(stderr, SDDS_EXIT_PrintErrors|SDDS_VERBOSE_PrintErrors);
    if (randomizeHOMFrequencies_struct.HOMFrequencyFilename) {
      iHOM=0;
      for (iRes=0; iRes<dipoleHOMTypes; iRes++) {
        resonator = &dipoleHOMResonator[iRes];
        if (resonator->cavities) 
          for (iCav=0; iCav<resonator->cavities; iCav++) {
            iHOM++;
          }
      }
      HOMs=iHOM;
      if (!SDDS_InitializeOutput(&randHOMFrequencyDataSet, SDDS_BINARY, 1,
                                 "Randomized HOM Frequencies", "Randomized HOM Frequencies",
                                 randomizeHOMFrequencies_struct.HOMFrequencyFilename)||
          (0>SDDS_DefineParameter(&randHOMFrequencyDataSet, "TotalCurrent", "I$btot$n", "mA",
                                  NULL, NULL, SDDS_DOUBLE, 0)) ||
          (0>SDDS_DefineParameter(&randHOMFrequencyDataSet, "Bunches", "N$bb$n", NULL,
                                  NULL, NULL, SDDS_LONG, 0)) ||
          (0>SDDS_DefineParameter(&randHOMFrequencyDataSet, "Sample", "Sample", NULL,
                                  NULL, NULL, SDDS_LONG, 0))||
          (0>SDDS_DefineColumn(&randHOMFrequencyDataSet, "HOMTypeIndex", "HOM type index", NULL,
                               NULL, NULL, SDDS_LONG, 0))||
          (0>SDDS_DefineColumn(&randHOMFrequencyDataSet, "CavityIndex", "Cavity index", NULL,
                               NULL, NULL, SDDS_LONG, 0))||
          (0>SDDS_DefineColumn(&randHOMFrequencyDataSet, "HOMFrequency", "HOM Frequency", "Hz",
                               NULL, NULL, SDDS_DOUBLE, 0))||
          !SDDS_WriteLayout(&randHOMFrequencyDataSet)||
          !SDDS_StartPage(&randHOMFrequencyDataSet, HOMs))
        SDDS_PrintErrors(stderr, SDDS_EXIT_PrintErrors|SDDS_VERBOSE_PrintErrors);
      if( !SDDS_SetParameters(&randHOMFrequencyDataSet, SDDS_BY_NAME|SDDS_PASS_BY_VALUE,
                              "TotalCurrent", totalCurrent,
                              "Bunches", totalFilledBuckets, NULL))
        SDDS_PrintErrors(stderr, SDDS_EXIT_PrintErrors|SDDS_VERBOSE_PrintErrors);
    }
    for (iSample=0; iSample<randomizeHOMFrequencies_struct.samples; iSample++) {
      iHOM=0;
      /* calculating all resonator frequencies */
      for (iRes=0; iRes<dipoleHOMTypes; iRes++) {
        resonator = &dipoleHOMResonator[iRes];
        if (resonator->cavities) {
          for (iCav=0; iCav<resonator->cavities; iCav++) {
            if (resonator->fixed) {
              resonator->frequency[iCav] = resonator->unperturbedFrequency;
            } else {
              resonator->frequency[iCav] = resonator->unperturbedFrequency + iCav * resonator->staggeringStep * dipoleHOMs_struct.staggeringStepMultiplier;
              /*  the shift to nearest revolution harmonic upper sideband is applied before randomization */
              if ( resonator->shiftToResonance )
                resonator->frequency[iCav] =
                  ((long)(((resonator->frequency[iCav]+betatronFrequency)/revolutionFrequency) + 0.5)) *
                  revolutionFrequency - betatronFrequency;
              resonator->frequency[iCav] += 
                randomizeHOMFrequencies_struct.spread * (random_1(0)-0.5);
            }
            
            if (randomizeHOMFrequencies_struct.HOMFrequencyFilename) {
              if (!SDDS_SetRowValues(&randHOMFrequencyDataSet, SDDS_BY_NAME|SDDS_PASS_BY_VALUE, iHOM,
                                     "HOMTypeIndex", iRes,
                                     "CavityIndex", iCav,
                                     "HOMFrequency", resonator->frequency[iCav], NULL))
                SDDS_PrintErrors(stderr, SDDS_EXIT_PrintErrors|SDDS_VERBOSE_PrintErrors);
              iHOM++;
            }
          }
        }
      }
      if (randomizeHOMFrequencies_struct.HOMFrequencyFilename)
        if (!SDDS_SetParameters(&randHOMFrequencyDataSet, SDDS_BY_NAME|SDDS_PASS_BY_VALUE,
                                "Sample", iSample, NULL) ||
            !SDDS_WritePage(&randHOMFrequencyDataSet))
          SDDS_PrintErrors(stderr, SDDS_EXIT_PrintErrors|SDDS_VERBOSE_PrintErrors);
      
      /* randomize bunch current if requested */
      if (randomizeBunchCurrentFlag) {
        bunchCurrentSave = (double *) malloc( totalFilledBuckets * sizeof(double));
        for (iBunch=0; iBunch<totalFilledBuckets; iBunch++)
          bunchCurrentSave[iBunch] = bunchCurrent[iBunch];
        doRandomizeBunchCurrent();
      }
      
      f = frequencyOfMaxTransGrowthRate();
      if (randomizeBunchCurrentFlag) {
        for (iBunch=0; iBunch<totalFilledBuckets; iBunch++)
          bunchCurrent[iBunch] = bunchCurrentSave[iBunch];
        free(bunchCurrentSave);      
      }
      
      if (!SDDS_SetRowValues(&maxGrowthRateDataSet, SDDS_BY_NAME|SDDS_PASS_BY_VALUE, iSample,
                             "MaxGrowthRate", f.i, NULL))
        SDDS_PrintErrors(stderr, SDDS_EXIT_PrintErrors|SDDS_VERBOSE_PrintErrors);
    }
    if (!SDDS_WritePage(&maxGrowthRateDataSet)||!SDDS_Terminate(&maxGrowthRateDataSet))
      SDDS_PrintErrors(stderr, SDDS_EXIT_PrintErrors|SDDS_VERBOSE_PrintErrors);
    if (randomizeHOMFrequencies_struct.HOMFrequencyFilename)
      if (!SDDS_Terminate(&randHOMFrequencyDataSet))
        SDDS_PrintErrors(stderr, SDDS_EXIT_PrintErrors|SDDS_VERBOSE_PrintErrors);
  }
  else {
    /*  One time straight calculation.   */
    /* randomize bunch current if requested */
    if (randomizeBunchCurrentFlag) {
      bunchCurrentSave = (double *) malloc( totalFilledBuckets * sizeof(double));
      for (iBunch=0; iBunch<totalFilledBuckets; iBunch++)
        bunchCurrentSave[iBunch] = bunchCurrent[iBunch];
      doRandomizeBunchCurrent();
    }
    
    f = frequencyOfMaxTransGrowthRate();
      
    if (randomizeBunchCurrentFlag) {
      for (iBunch=0; iBunch<totalFilledBuckets; iBunch++)
        bunchCurrent[iBunch] = bunchCurrentSave[iBunch];
      free(bunchCurrentSave);      
    }
    
    printf("Maximum growth rate is %lf 1/sec.\n", f.i);
    /* list the frequencies of all CBMs */
    if (!doTransverseMotion_struct.CBMFrequencyFilename) {
      bomb("Filename CBMFrequencyFilename not specified in namelist doTransverseMotion.", NULL);
    }            
    if (!SDDS_InitializeOutput(&CBMEigenvalues, SDDS_BINARY, 1,
                               "CBM eigenvalues", "CBM eigenvalues",
                               doTransverseMotion_struct.CBMFrequencyFilename)||
        (0>SDDS_DefineParameter(&CBMEigenvalues, "TotalCurrent", "I$btot$n", "mA",
                                NULL, NULL, SDDS_DOUBLE, 0)) ||
        (0>SDDS_DefineParameter(&CBMEigenvalues, "Bunches", "N$bb$n", NULL,
                                NULL, NULL, SDDS_LONG, 0)) ||
        (0>SDDS_DefineColumn(&CBMEigenvalues, "CBMIndex", "CBM Index", NULL, "CBM index", NULL, SDDS_LONG, 0))||
        (0>SDDS_DefineColumn(&CBMEigenvalues, "RealCBMFrequency", "Re{$gw$r}", NULL,
                             "Real part of CBM frequency", NULL, SDDS_DOUBLE, 0))||
        (0>SDDS_DefineColumn(&CBMEigenvalues, "ImagCBMFrequency", "Im{$gw$r}", NULL,
                             "Imag part of CBM frequency", NULL, SDDS_DOUBLE, 0))||
        (0>SDDS_DefineColumn(&CBMEigenvalues, "RealDeltaCBMFrequency", "Re{$gDw$r}", "rad/sec",
                             "Real part of CBM frequency shift", NULL, SDDS_DOUBLE, 0))||
        (0>SDDS_DefineColumn(&CBMEigenvalues, "ImagDeltaCBMFrequency", "Im{$gDw$r}", "rad/sec",
                             "Imaginary part of CBM frequency shift", NULL, SDDS_DOUBLE, 0))||
        (0>SDDS_DefineColumn(&CBMEigenvalues, "MagDeltaCBMFrequency", "$sb$e$gDw$r$sb$e", "rad/sec",
                             "Magnitude of frequency shift", NULL, SDDS_DOUBLE, 0))||
        !SDDS_WriteLayout(&CBMEigenvalues)||
        !SDDS_StartPage(&CBMEigenvalues, totalFilledBuckets))
      SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
    if( !SDDS_SetParameters(&CBMEigenvalues, SDDS_BY_NAME|SDDS_PASS_BY_VALUE,
                            "TotalCurrent", totalCurrent,
                            "Bunches", totalFilledBuckets, NULL))
      SDDS_PrintErrors(stderr, SDDS_EXIT_PrintErrors|SDDS_VERBOSE_PrintErrors);
    
    for (i=0; i<totalFilledBuckets; i++) {
      //fDiff=csub(eigenvalues[i], f0);
      fDiff.r = lapack_eigenvalues[i].real - f0.r;
      fDiff.i = lapack_eigenvalues[i].imag - f0.i;
      if (!SDDS_SetRowValues(&CBMEigenvalues, SDDS_BY_NAME|SDDS_PASS_BY_VALUE, i,
                             "CBMIndex", i,
                             "RealCBMFrequency", lapack_eigenvalues[i].real,
                             "ImagCBMFrequency", lapack_eigenvalues[i].imag,
                             "RealDeltaCBMFrequency", fDiff.r,
                             "ImagDeltaCBMFrequency", fDiff.i,
                             "MagDeltaCBMFrequency", cmod(fDiff), NULL))
        SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
    }
    if (!SDDS_WritePage(&CBMEigenvalues)||!SDDS_Terminate(&CBMEigenvalues))
      SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
  }
  /* reset randomize flag */
  randomizeHOMFlag = 0;
  /* reset sweep flag */
  sweepFlag = 0;
}

void DefaultPhysicalParameters() {
  bunchLengtheningFactor=1;
  return;
}

void calculateSomePhysicalParameters() {
  double q; /* over voltage */
  revolutionFrequency = c_mks/ circumference * (1 - sqr(me_mev/(energyGeV*1e3)));
  rfFrequency = revolutionFrequency * harmonicNumber;
  q = rfVoltageMV/ energyLossPerTurnMeV;
  synchronousPhase = asin( 1/q );
  rfAcceptance = sqrt(2*energyLossPerTurnMeV/PI/momentumCompaction/harmonicNumber/(energyGeV*1e3) * (sqrt( sqr(q) - 1) - acos(1/q)  ));
  synchrotronAngFrequency = sqrt( 2*PI * momentumCompaction * harmonicNumber * 
                                  cos(synchronousPhase) *
                                  rfVoltageMV/ (energyGeV*1e3) ) * revolutionFrequency;
  synchrotronFrequency = synchrotronAngFrequency /2 /PI;
  bunchLength = relativeEnergySpread * momentumCompaction * c_mks/ synchrotronAngFrequency;
  printf("Some calculated physical quantities:\n" "\trevolution frequency (MHz):\t%lf\n" "\tRF frequency (MHz):\t%lf"
         "\n\tsynchronous phase (deg):\t%lf\n" "\tsynchrotron angular frequency (rad/sec):\t%lf\n"
         "\tsynchrotron frequency (Hz):\t%lf\n" "\tzero-current bunch length (m):\t%lf\n",
         revolutionFrequency/1e6, rfFrequency/1e6, synchronousPhase/PI*180, synchrotronAngFrequency, synchrotronFrequency,
         bunchLength);
  bunchLength *= bunchLengtheningFactor;
  printf("\tbunch length used (m):\t%lf\n", bunchLength);
}

double powerOfResonators() {
  long i, j, bunches;
  long iRes, iCav;
  HOM_RESONATOR *resonator;
  double *formFactor, wakePot;
  static COMPLEX cZero={0.0, 0.0};
  double Q, Tf, timeElapsed, charge, power;
  COMPLEX argument, Vb, delVb, Vbplus, Vbminus;

  /* shorter name for the same thing */
  bunches = totalFilledBuckets;
  formFactor = (double *) malloc(bunches*sizeof(double));
  power = 0.0 ;
  for (iRes=0; iRes<monopoleHOMTypes; iRes++) {
    resonator = &monopoleHOMResonator[iRes];
    if (!monopoleHOMs_struct.Q&&(!resonator->Q||!resonator->RoQ))
      continue;
    if (resonator->cavities) 
      for (iCav=0; iCav<resonator->cavities; iCav++) {
        /* form factor will be included in the bunch j-index loop */
        if (!bunchLengthPicoSecondData) {
          for (j=0; j<bunches; j++) {
            /* double bunchLength has dimension of length (units m) */
            formFactor[j] = exp( - sqr(bunchLength/c_mks) * sqr(2 * PI * resonator->frequency[iCav]));
          }
        }
        else {
          for (j=0; j<bunches; j++)
            /* double array bunchDurationByBunch has dimension of time (units sec) */
            /* why is there a factor of 0.5 here but not above? */
            formFactor[j] = exp( - 0.5 * sqr(bunchDurationByBunch[j]) * sqr(2 * PI * resonator->frequency[iCav]));
        }
        /*  K.Thompson's wake is defined as a negative number while I like to use positive numbers
            so I add a minus sign for the wakePot variable here.    */
        wakePot = - 2*PI*resonator->frequency[iCav] * resonator->RoQ ;
	/*	fprintf(stdout, "frequency %g  R/q %g\n", resonator->frequency[iCav], resonator->RoQ); */
	Vb = cZero;
	Vbplus = cZero;
	Vbminus = cZero;
        for (i=0; i<bunches; i++) {
          /* Here, inside the single loop I can sum up the voltage vectors of
             one (or more) turn of bunches. See Fig. 6.3 of SLAC-PUB-2884,
             Perry Wilson. I only sum one turn, each bunch contributing once.
             The formula 6.25a doesn't apply generally since we may have more
             complicated bunch pattern. For symmetric patterns the closed form 
             may be convenient, but here we can afford the little CPU time of 
             making explict sums. So we use the first part of 6.25a */
	  /* Need to calculate V_b as bunch-by-bunch passage as
	     represented in Figure 6.3 
	     to estimate power losses for each mode */
	  /* Tf: time to fill this particular mode */
	  if (!monopoleHOMs_struct.Q)
	    Q = resonator->Q/(resonator->deQFactor * monopoleHOMs_struct.deQFactorMultiplier);          
	  else
	    Q = monopoleHOMs_struct.Q;
	  Tf = Q / PI / resonator->frequency[iCav];
	  timeElapsed = (1.0 * bucketSelected[i]) / harmonicNumber / revolutionFrequency;
	  /*	  fprintf (stdout, "bucketSelected[%ld] %ld\n", i,  bucketSelected[i]); */
	  argument.r = - timeElapsed / Tf;
	  argument.i = 2 * PI * timeElapsed * (resonator->frequency[iCav] - rfFrequency);
	  /*	  fprintf (stdout, "argument (%lf, %lf)\n", argument.r, argument.i ); */
	  charge = (bunchCurrent[i] * 1e-3) / revolutionFrequency;
	  /*	  fprintf (stdout, "bunchCurrent[%ld] %lf mA \n", i,  bunchCurrent[i]); */
	  delVb = cexp_oag( argument);
	  delVb.r = - wakePot * formFactor[i] * charge * delVb.r;
	  delVb.i = - wakePot * formFactor[i] * charge * delVb.i;
	  /*	  fprintf (stdout, "bunch %ld wakePot %lf charge %g delVb.r %g delVb.i %g\n", i, wakePot, charge, delVb.r, delVb.i); 
	   */
	  if (i==0) {
	    Vb.r += 0.5 * delVb.r;   /* voltage at passage, i.e. beam-loading voltage */
	    Vb.i += 0.5 * delVb.i;
	    Vbplus.r += delVb.r;  
	    Vbplus.i += delVb.i;
	  } 
	  else {
	    Vb.r +=  delVb.r;
	    Vb.i +=  delVb.i;
	    Vbplus.r +=  delVb.r;   /* voltage after passage */
	    Vbplus.i +=  delVb.i;
	    Vbminus.r +=  delVb.r;  /* voltage before passage */
	    Vbminus.i +=  delVb.i;
	  }
	}
	/* each resonator (many in a cavity) will have its own power. The voltages of modes
	   are uncorrelated because they have wildly different frequencies, not to mention
	   wildly different spatial patterns. */
	resonator->power[iCav] = (sqr(Vb.r) + sqr(Vb.i))/( 2 * resonator->shuntImpedance);
	resonator->powerPlus[iCav] = (sqr(Vbplus.r) + sqr(Vbplus.i))/( 2 * resonator->shuntImpedance);
	resonator->powerMinus[iCav] = (sqr(Vbminus.r) + sqr(Vbminus.i))/( 2 * resonator->shuntImpedance);
	/*	fprintf (stdout, "iRes %ld iCav %ld Power ave %lf  PowerPlus %lf  PowerMinus %lf.\n", iRes, iCav, resonator->power[iCav],resonator->powerPlus[iCav],resonator->powerMinus[iCav]); */
	/* sum the powers */
	power += resonator->power[iCav];
	/* fprintf (stdout, "Power %lf W.\n", power); */
      }
  }
  /* fprintf (stdout, "Power %lf W.\n", power); */
  return power;
}

COMPLEX frequencyOfMaxLongGrowthRate() {
  COMPLEX t;
  lapack_complex_double *chi;
  HOM_RESONATOR *resonator;
  double *formFactor, wakePot;
  COMPLEX waveNumber, cArg, cArg1, denom1Term, denom2Term;
  COMPLEX numer1Term, numer2Term, numer3Term, numer4Term;
  COMPLEX ik, mik;
  static COMPLEX I={0.0, 1.0};
  static COMPLEX mI={0.0, -1.0};
  static COMPLEX one={1.0, 0.0};
  long iRes, iCav, iCBM, iMax;
  long bunches;
  long i, j;
  double maxGrowthRate;
  char modeName[32];
  SDDS_TABLE eigenvectorPage;
  double sum;
  long error_code;
  static int eigenvectorFileWritten=0;
  double tmp;
  lapack_complex_double *lapack_m;
  COMPLEX result;
  lapack_complex_double *lapack_eigenvectors=NULL;
  double *eigenvectors_real, *eigenvectors_imag;

  /* shorter name for the same thing */
  bunches = totalFilledBuckets;
  chi = malloc(sizeof(lapack_complex_double) * bunches * bunches);
#pragma omp parallel num_threads(threads) default(none) shared(chi) firstprivate(bunches) private(i)
  {
#pragma omp for
    for (i=0; i<bunches*bunches; i++) {
      chi[i].real=0;
      chi[i].imag=0;
    }
  }
  formFactor = (double *) malloc(bunches*sizeof(double));
  /* I assume that chi is initialized to zero values */    
  for (iRes=0; iRes<monopoleHOMTypes; iRes++) {
    resonator = &monopoleHOMResonator[iRes];
    if (!monopoleHOMs_struct.Q&&(!resonator->Q||!resonator->RoQ))
      continue;
    if (resonator->cavities) 
      for (iCav=0; iCav<resonator->cavities; iCav++) {
        /* form factor will be included in the bunch j-index loop */
        if (!bunchLengthPicoSecondData) {
          tmp =  exp( - sqr(bunchLength/c_mks) * sqr(2 * PI * resonator->frequency[iCav]));
#pragma omp parallel num_threads(threads) default(none) shared(formFactor) firstprivate(bunches,tmp) private(j)
          {
#pragma omp for
            for (j=0; j<bunches; j++) {
              /* double bunchLength has dimension of length (units m) */
              formFactor[j] = tmp;
            }
          }
        }
        else {
          tmp = sqr(2 * PI * resonator->frequency[iCav]);
#pragma omp parallel num_threads(threads) default(none) shared(formFactor,bunchDurationByBunch) firstprivate(bunches,tmp) private(j)
          {
#pragma omp for
            for (j=0; j<bunches; j++) {
              /* double array bunchDurationByBunch has dimension of time (units sec) */
              formFactor[j] = exp( - sqr(bunchDurationByBunch[j]) * tmp);
            }
          }
        }
        /*  K.Thompson's wake is defined as a negative number while I like to use positive numbers
            so I add a minus sign for the wakePot variable here.    */
        wakePot = - 2*PI*resonator->frequency[iCav] * resonator->RoQ ;
        waveNumber.r = 2*PI*resonator->frequency[iCav]/ c_mks;
        if (!monopoleHOMs_struct.Q)
          waveNumber.i = waveNumber.r/ (2.0 * resonator->Q/(resonator->deQFactor * monopoleHOMs_struct.deQFactorMultiplier));
        else
          waveNumber.i = waveNumber.r/ (2.0 * monopoleHOMs_struct.Q);
        ik = cmul(I, waveNumber);
        mik = cmul(mI, cconj(waveNumber));
        cArg1 = cmulr(I, synchrotronAngFrequency/ revolutionFrequency);
        cArg = cadd(cmulr(ik, circumference), cArg1);
        denom1Term = cexp_oag(cArg);
        cArg = cadd(cmulr(mik, circumference), cArg1);
        denom2Term = cexp_oag(cArg);
#pragma omp parallel num_threads(threads) default(none) private(i,j,numer1Term,numer2Term,numer3Term,numer4Term, t) firstprivate(bunches, ik, circumference, harmonicNumber, mik, one, denom1Term, denom2Term, cArg1, wakePot) shared(chi,bucketSelected, formFactor)
        {
          double r, exp_x[2], cr[2], ci[2], cos_x[2], sin_x[2];
          COMPLEX c, c2;
#pragma omp for
          for (i=0; i<bunches; i++) {
            for (j=0; j<bunches; j++) {
              if (i>j) { 
                //numer1Term = cmul( ik, cexp_oag(cmulr(ik, (bucketSelected[i]-bucketSelected[j])*circumference/ harmonicNumber)));
                //numer2Term = cmul( mik, cexp_oag(cmulr(mik, (bucketSelected[i]-bucketSelected[j])*circumference/ harmonicNumber)));
                //t = cadd(cdiv(numer1Term, csub(one, denom1Term)), cdiv(numer2Term, csub(one, denom2Term)));
                r = (bucketSelected[i]-bucketSelected[j])*circumference/ harmonicNumber;
                cr[0] = ik.r * r;
                ci[0] = ik.i * r;
                cr[1] = mik.r * r;
                ci[1] = mik.i * r;
                  
                vdExp(2,cr,exp_x);
                vdSinCos(2,ci,sin_x,cos_x);
                
                cr[0] = exp_x[0] * cos_x[0];
                ci[0] = exp_x[0] * sin_x[0];
                numer1Term.r = ik.r * cr[0] - ik.i * ci[0];
                numer1Term.i = ik.r * ci[0] + ik.i * cr[0];

                cr[1] = exp_x[1] * cos_x[1];
                ci[1] = exp_x[1] * sin_x[1];
                numer2Term.r = mik.r * cr[1] - mik.i * ci[1];
                numer2Term.i = mik.r * ci[1] + mik.i * cr[1];

                c.r = one.r - denom1Term.r;
                c.i = one.i - denom1Term.i;
                r = c.r * c.r + c.i * c.i;
                if (r == 0) {
                  bombClinchor("division by zero in cdiv()");
                }
                c2.r = (numer1Term.r * c.r + numer1Term.i * c.i)/r;
                c2.i = (-numer1Term.r * c.i + numer1Term.i * c.r)/r;

                c.r = one.r - denom2Term.r;
                c.i = one.i - denom2Term.i;
                r = c.r * c.r + c.i * c.i;
                if (r == 0) {
                  bombClinchor("division by zero in cdiv()");
                }
                t.r = (numer2Term.r * c.r + numer2Term.i * c.i)/r + c2.r;
                t.i = (-numer2Term.r * c.i + numer2Term.i * c.r)/r + c2.i;
              }
              else if (i<=j) {
                //numer3Term = cmul( ik, cexp_oag( cadd( cmulr( ik, circumference + (bucketSelected[i]-bucketSelected[j])*circumference/ harmonicNumber), cArg1)));
                //numer4Term = cmul( mik, cexp_oag( cadd( cmulr( mik, circumference + (bucketSelected[i]-bucketSelected[j])*circumference/ harmonicNumber), cArg1)));
                //t = cadd( cdiv( numer3Term, csub( one, denom1Term)), cdiv(numer4Term, csub( one, denom2Term)));
                r = circumference + (bucketSelected[i]-bucketSelected[j])*circumference/ harmonicNumber;
                cr[0] = ik.r * r + cArg1.r;
                ci[0] = ik.i * r + cArg1.i;
                cr[1] = mik.r * r + cArg1.r;
                ci[1] = mik.i * r + cArg1.i;

                vdExp(2,cr,exp_x);
                vdSinCos(2,ci,sin_x,cos_x);

                cr[0] = exp_x[0] * cos_x[0];
                ci[0] = exp_x[0] * sin_x[0];
                numer3Term.r = ik.r * cr[0] - ik.i * ci[0];
                numer3Term.i = ik.r * ci[0] + ik.i * cr[0];

                cr[1] = exp_x[1] * cos_x[1];
                ci[1] = exp_x[1] * sin_x[1];
                numer4Term.r = mik.r * cr[1] - mik.i * ci[1];
                numer4Term.i = mik.r * ci[1] + mik.i * cr[1];

                c.r = one.r - denom1Term.r;
                c.i = one.i - denom1Term.i;
                r = c.r * c.r + c.i * c.i;
                if (r == 0) {
                  bombClinchor("division by zero in cdiv()");
                }
                c2.r = (numer3Term.r * c.r + numer3Term.i * c.i)/r;
                c2.i = (-numer3Term.r * c.i + numer3Term.i * c.r)/r;

                c.r = one.r - denom2Term.r;
                c.i = one.i - denom2Term.i;
                r = c.r * c.r + c.i * c.i;
                if (r == 0) {
                  bombClinchor("division by zero in cdiv()");
                }
                t.r = (numer4Term.r * c.r + numer4Term.i * c.i)/r + c2.r;
                t.i = (-numer4Term.r * c.i + numer4Term.i * c.r)/r + c2.i;
              }
              chi[i*bunches+j].real += wakePot/2 * t.r * formFactor[j];
              chi[i*bunches+j].imag += wakePot/2 * t.i * formFactor[j];
            }
          }
        }
      }
  }
#pragma omp parallel num_threads(threads) default(none) shared(chi,bunchCurrent) firstprivate(bunches,momentumCompaction,energyGeV) private(i,j,tmp)
  {
#pragma omp for
    for (j=0; j<bunches; j++) {
      tmp = -momentumCompaction * bunchCurrent[j]/ 1e3 * c_mks / energyGeV/ 1e9;
      for (i=0; i<bunches; i++) {
        chi[i*bunches+j].real *= tmp;
        chi[i*bunches+j].imag *= tmp;
      }
    }
  }
  lapack_m = malloc(sizeof(lapack_complex_double) * bunches * bunches);
#pragma omp parallel num_threads(threads) default(none) shared(chi,lapack_m) firstprivate(bunches,synchrotronAngFrequency,dampingFactor) private(i,j)
  {
#pragma omp for
    for (i=0; i<bunches; i++) {
      for (j=0; j<bunches; j++) {
        lapack_m[i*bunches+j].real = - chi[i*bunches+j].real / 2.0/ synchrotronAngFrequency;
        lapack_m[i*bunches+j].imag = - chi[i*bunches+j].imag / 2.0/ synchrotronAngFrequency;
      }
      lapack_m[i*bunches+i].real += synchrotronAngFrequency;
      lapack_m[i*bunches+i].imag -= dampingFactor; // K. Thompson uses the term lambda/2 where lambda is twice the damping rate.
    }
  }

  free(chi);

  if ((doLongitudinalMotion_struct.eigenvectorFilename) && (eigenvectorFileWritten==0)) {
    lapack_eigenvectors = (lapack_complex_double*) tmalloc(sizeof(lapack_complex_double)*totalFilledBuckets*totalFilledBuckets);
    error_code = LAPACKE_zgeev( LAPACK_ROW_MAJOR, 'N', 'V', bunches, lapack_m, bunches, lapack_eigenvalues, NULL, totalFilledBuckets, lapack_eigenvectors, totalFilledBuckets);
  } else {
    error_code = LAPACKE_zgeev( LAPACK_ROW_MAJOR, 'N', 'N', bunches, lapack_m, bunches, lapack_eigenvalues, NULL, totalFilledBuckets, NULL, totalFilledBuckets);
  }
  if (error_code != 0) {
    fprintf(stderr,"frequencyOfMaxLongGrowthRate: Subroutine cg returned error code %ld.\n", error_code);
    bombClinchor(NULL);
  }

  if ((doLongitudinalMotion_struct.eigenvectorFilename) && (eigenvectorFileWritten==0)) {
    eigenvectorFileWritten=1;
    /* normalize the eigenvectors */
#pragma omp parallel num_threads(threads) default(none) shared(lapack_eigenvectors) firstprivate(totalFilledBuckets) private(i,j,sum)
    {
#pragma omp for
      for (j=0; j<totalFilledBuckets; j++) {
        sum = 0;
        for (i=0; i<totalFilledBuckets; i++) {
          sum += sqr(lapack_eigenvectors[i*totalFilledBuckets+j].real) + sqr(lapack_eigenvectors[i*totalFilledBuckets+j].imag);
        }
        sum = sqrt(sum/totalFilledBuckets);
        for (i=0; i<totalFilledBuckets; i++) {
          lapack_eigenvectors[i*totalFilledBuckets+j].real /= sum;
          lapack_eigenvectors[i*totalFilledBuckets+j].imag /= sum;
        }
      }
    }
    if( !SDDS_InitializeOutput(&eigenvectorPage, SDDS_BINARY, 1,
                               "CBM eigenvectors", "CBM eigenvectors",
                               doLongitudinalMotion_struct.eigenvectorFilename ))
      SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
    if(0>SDDS_DefineColumn(&eigenvectorPage, "bunchIndex", NULL, NULL, NULL, NULL, SDDS_LONG, 0))
      SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
    for (i=0; i<bunches; i++) {
      sprintf(modeName, "mode%04ldReal", i);
      if(0>SDDS_DefineColumn(&eigenvectorPage, modeName, NULL, NULL, NULL, NULL, SDDS_DOUBLE, 0))
        SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
      sprintf(modeName, "mode%04ldImag", i);
      if(0>SDDS_DefineColumn(&eigenvectorPage, modeName, NULL, NULL, NULL, NULL, SDDS_DOUBLE, 0))
        SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
    }
    if( !SDDS_WriteLayout(&eigenvectorPage) ||
        !SDDS_StartPage(&eigenvectorPage, bunches))
      SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
    for (i=0; i<bunches; i++) {
      if (!SDDS_SetRowValues(&eigenvectorPage, SDDS_BY_NAME|SDDS_PASS_BY_VALUE, i,
                             "bunchIndex", i, NULL))
        SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
    }
    
    /* to output a matrix one has to write columns of the transpose! */
    eigenvectors_real = (double*) malloc(sizeof(double)*totalFilledBuckets);
    eigenvectors_imag = (double*) malloc(sizeof(double)*totalFilledBuckets); 
    
    for (i=0; i<bunches; i++) {
      for (j=0; j<totalFilledBuckets; j++) {
        eigenvectors_real[j] = lapack_eigenvectors[j*bunches+i].real;
        eigenvectors_imag[j] = lapack_eigenvectors[j*bunches+i].imag;
      }
      sprintf(modeName, "mode%04ldReal", i);
      if (!SDDS_SetColumn(&eigenvectorPage, SDDS_SET_BY_NAME|SDDS_PASS_BY_REFERENCE,
                          eigenvectors_real, bunches, modeName))
        SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
      sprintf(modeName, "mode%04ldImag", i);
      if (!SDDS_SetColumn(&eigenvectorPage, SDDS_SET_BY_NAME|SDDS_PASS_BY_REFERENCE,
                          eigenvectors_imag, bunches, modeName))
        SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
    }
    if (!SDDS_WritePage(&eigenvectorPage) || !SDDS_Terminate(&eigenvectorPage))
      SDDS_PrintErrors(stderr, SDDS_EXIT_PrintErrors|SDDS_VERBOSE_PrintErrors);
    free(eigenvectors_real);
    free(eigenvectors_imag);
    free(lapack_eigenvectors);
  }
  free(formFactor);
  /* find the eigenvalue representing the fastest growing CBM */
  maxGrowthRate=-1e38;
  iMax = 0;
  for (iCBM=0; iCBM<bunches; iCBM++) {
    if (maxGrowthRate<lapack_eigenvalues[iCBM].imag) {
      maxGrowthRate=lapack_eigenvalues[iCBM].imag;
      iMax=iCBM;  
    }
  }
  result.r = lapack_eigenvalues[iMax].real;
  result.i = lapack_eigenvalues[iMax].imag;
  return(result);
}

COMPLEX frequencyOfMaxTransGrowthRate() {
  COMPLEX t;
  lapack_complex_double *chi;
  HOM_RESONATOR *resonator;
  double *formFactor, wakePot;
  COMPLEX waveNumber, cArg, cArg1;
  COMPLEX numer1Term, numer2Term, numer3Term, numer4Term, denom1Term, denom2Term;
  COMPLEX ik, mik;
  static COMPLEX I={0.0, 1.0};
  static COMPLEX mI={0.0, -1.0};
  static COMPLEX one={1.0, 0.0};
  long iRes, iCav, iCBM, iMax=0;
  long bunches;
  long i, j;
  double maxGrowthRate;
  char modeName[32];
  SDDS_TABLE eigenvectorPage;
  double sum;
  int error_code;
  double tmp;
  static int eigenvectorFileWritten=0;
  lapack_complex_double *lapack_m;
  COMPLEX result;
  lapack_complex_double *lapack_eigenvectors=NULL;
  double *eigenvectors_real, *eigenvectors_imag;
  /* shorter name for the same thing */
  bunches = totalFilledBuckets;
  chi = malloc(sizeof(lapack_complex_double) * bunches * bunches);
#pragma omp parallel num_threads(threads) default(none) shared(chi) firstprivate(bunches) private(i)
  {
#pragma omp for
    for (i=0; i<bunches*bunches; i++) {
      chi[i].real=0;
      chi[i].imag=0;
    }
  }
  formFactor = (double *) malloc(bunches*sizeof(double));
  /* I assume that chi is initialized to zero values */    
  for (iRes=0; iRes<dipoleHOMTypes; iRes++) {
    resonator = &dipoleHOMResonator[iRes];
    if (!dipoleHOMs_struct.Q&&(!resonator->Q||!resonator->RoQ))
      continue;
    if (resonator->cavities) {
      for (iCav=0; iCav<resonator->cavities; iCav++) {
        if (!bunchLengthPicoSecondData) {
          tmp = exp( - sqr(bunchLength/c_mks) * sqr(2 * PI * resonator->frequency[iCav]) );
#pragma omp parallel num_threads(threads) default(none) shared(formFactor) firstprivate(bunches,tmp) private(j)
          {
#pragma omp for
            for (j=0; j<bunches; j++) {
              /* double bunchLength has dimension of length (units m) */
              /* form factor doesn't have a 2 in denominator because
                 the voltage generated by a gaussian comes from the
                 charge squared */
              formFactor[j] = tmp;
            }
          }
        }
        else {
          tmp = sqr(2 * PI * resonator->frequency[iCav]);
#pragma omp parallel num_threads(threads) default(none) shared(formFactor,bunchDurationByBunch) firstprivate(bunches,tmp) private(j)
          {
#pragma omp for
            for (j=0; j<bunches; j++) {
              /* double array bunchDurationByBunch has dimension of time (units sec) */
              /* form factor doesn't have a 2 in denominator because
                 the voltage generated by a gaussian comes from the
                 charge squared */
              formFactor[j] = exp( - sqr(bunchDurationByBunch[j]) * tmp );
            }
          }
        }
        /*  K.Thompson's wake is defined as a positive number as is mine
            i.e. delta x' = W x with W>0 
            so I don't add a minus sign for the wake_pot variable here.    */
        wakePot = 2*PI*resonator->frequency[iCav] * resonator->RoQ ;
        waveNumber.r = 2*PI*resonator->frequency[iCav]/ c_mks;
        if (!dipoleHOMs_struct.Q)
          waveNumber.i = waveNumber.r/ 2.0/ 
            (resonator->Q/resonator->deQFactor/dipoleHOMs_struct.deQFactorMultiplier);
        else
          waveNumber.i =  waveNumber.r/ (2.0 * dipoleHOMs_struct.Q);
        ik = cmul(I, waveNumber);
        mik = cmul(mI, cconj(waveNumber));
        cArg1 = cmulr(I, betatronAngFrequency/ revolutionFrequency);
        cArg = cadd(cmulr(ik, circumference), cArg1);
        denom1Term = cexp_oag(cArg);
        cArg = cadd(cmulr(mik, circumference), cArg1);
        denom2Term = cexp_oag(cArg);
#pragma omp parallel num_threads(threads) default(none) private(i,j,numer1Term,numer2Term,numer3Term,numer4Term, t) firstprivate(bunches, I, ik, circumference, harmonicNumber, mI, mik, one, denom1Term, denom2Term, cArg1, wakePot) shared(chi,bucketSelected, formFactor)
        {
          double r, exp_x[2], cr[2], ci[2], cos_x[2], sin_x[2];
          COMPLEX c, c2;
#pragma omp for
          for (i=0; i<bunches; i++) {
            for (j=0; j<bunches; j++) {
              if (i>j) {
                
                //numer1Term = cmul( I, cexp_oag(cmulr(ik, (bucketSelected[i]-bucketSelected[j])*circumference/ harmonicNumber)));
                //numer2Term = cmul( mI, cexp_oag(cmulr(mik, (bucketSelected[i]-bucketSelected[j])*circumference/ harmonicNumber)));
                //t = cadd( cdiv( numer1Term, csub( one, denom1Term)), cdiv( numer2Term, csub( one, denom2Term)));
                
                r = (bucketSelected[i]-bucketSelected[j])*circumference/ harmonicNumber;
                cr[0] = ik.r * r;
                ci[0] = ik.i * r;
                cr[1] = mik.r * r;
                ci[1] = mik.i * r;
                  
                vdExp(2,cr,exp_x);
                vdSinCos(2,ci,sin_x,cos_x);

                cr[0] = exp_x[0] * cos_x[0];
                ci[0] = exp_x[0] * sin_x[0];

                numer1Term.r = I.r * cr[0] - I.i * ci[0];
                numer1Term.i = I.r * ci[0] + I.i * cr[0];

                cr[1] = exp_x[1] * cos_x[1];
                ci[1] = exp_x[1] * sin_x[1];
                numer2Term.r = mI.r * cr[1] - mI.i * ci[1];
                numer2Term.i = mI.r * ci[1] + mI.i * cr[1];

                c.r = one.r - denom1Term.r;
                c.i = one.i - denom1Term.i;
                r = c.r * c.r + c.i * c.i;
                if (r == 0) {
                  bombClinchor("division by zero in cdiv()");
                }
                c2.r = (numer1Term.r * c.r + numer1Term.i * c.i)/r;
                c2.i = (-numer1Term.r * c.i + numer1Term.i * c.r)/r;

                c.r = one.r - denom2Term.r;
                c.i = one.i - denom2Term.i;
                r = c.r * c.r + c.i * c.i;
                if (r == 0) {
                  bombClinchor("division by zero in cdiv()");
                }
                t.r = (numer2Term.r * c.r + numer2Term.i * c.i)/r + c2.r;
                t.i = (-numer2Term.r * c.i + numer2Term.i * c.r)/r + c2.i;
              }
              else if (i<=j) {
                
                //numer3Term = cmul( I, cexp_oag( cadd( cmulr( ik, circumference + (bucketSelected[i]-bucketSelected[j])*circumference/ harmonicNumber), cArg1)));
                //numer4Term = cmul( mI, cexp_oag( cadd( cmulr( mik, circumference + (bucketSelected[i]-bucketSelected[j])*circumference/ harmonicNumber), cArg1)));
                //t = cadd( cdiv( numer3Term, csub( one, denom1Term)), cdiv(numer4Term, csub( one, denom2Term)));                

                r = circumference + (bucketSelected[i]-bucketSelected[j])*circumference/ harmonicNumber;
                cr[0] = ik.r * r + cArg1.r;
                ci[0] = ik.i * r + cArg1.i;
                cr[1] = mik.r * r + cArg1.r;
                ci[1] = mik.i * r + cArg1.i;

                vdExp(2,cr,exp_x);
                vdSinCos(2,ci,sin_x,cos_x);

                cr[0] = exp_x[0] * cos_x[0];
                ci[0] = exp_x[0] * sin_x[0];
                numer3Term.r = I.r * cr[0] - I.i * ci[0];
                numer3Term.i = I.r * ci[0] + I.i * cr[0];

                cr[1] = exp_x[1] * cos_x[1];
                ci[1] = exp_x[1] * sin_x[1];
                numer4Term.r = mI.r * cr[1] - mI.i * ci[1];
                numer4Term.i = mI.r * ci[1] + mI.i * cr[1];

                c.r = one.r - denom1Term.r;
                c.i = one.i - denom1Term.i;
                r = c.r * c.r + c.i * c.i;
                if (r == 0) {
                  bombClinchor("division by zero in cdiv()");
                }
                c2.r = (numer3Term.r * c.r + numer3Term.i * c.i)/r;
                c2.i = (-numer3Term.r * c.i + numer3Term.i * c.r)/r;

                c.r = one.r - denom2Term.r;
                c.i = one.i - denom2Term.i;
                r = c.r * c.r + c.i * c.i;
                if (r == 0) {
                  bombClinchor("division by zero in cdiv()");
                }
                t.r = (numer4Term.r * c.r + numer4Term.i * c.i)/r + c2.r;
                t.i = (-numer4Term.r * c.i + numer4Term.i * c.r)/r + c2.i;
              }
              chi[i*bunches+j].real += wakePot/2 * t.r * formFactor[j];
              chi[i*bunches+j].imag += wakePot/2 * t.i * formFactor[j];
            }
          }
        }
      }
    }
  }
#pragma omp parallel num_threads(threads) default(none) shared(chi,bunchCurrent) firstprivate(bunches,betatronAngFrequency,betaAtRFCavs,energyGeV) private(i,j,tmp)
  {
#pragma omp for
    for (j=0; j<bunches; j++) {
      tmp = -bunchCurrent[j]/ 1e3 * betatronAngFrequency * betaAtRFCavs / energyGeV/ 1e9;
      for (i=0; i<bunches; i++) {
        chi[i*bunches+j].real *= tmp;
        chi[i*bunches+j].imag *= tmp;
      }
    }
  }
  
  lapack_m = malloc(sizeof(lapack_complex_double) * bunches * bunches);
#pragma omp parallel num_threads(threads) default(none) shared(chi,lapack_m) firstprivate(bunches,betatronAngFrequency,dampingFactor) private(i,j)
  {
#pragma omp for
    for (i=0; i<bunches; i++) {
      for (j=0; j<bunches; j++) {
        lapack_m[i*bunches+j].real = - chi[i*bunches+j].real / 2.0/ betatronAngFrequency;
        lapack_m[i*bunches+j].imag = - chi[i*bunches+j].imag / 2.0/ betatronAngFrequency;
      }
      lapack_m[i*bunches+i].real += betatronAngFrequency;
      lapack_m[i*bunches+i].imag -= dampingFactor; // K. Thompson uses the term lambda/2 where lambda is twice the damping rate.
    }
  }

  free(chi);

  if ((doTransverseMotion_struct.eigenvectorFilename) && (eigenvectorFileWritten==0)) {
    lapack_eigenvectors = (lapack_complex_double*) tmalloc(sizeof(lapack_complex_double)*totalFilledBuckets*totalFilledBuckets);
    error_code = LAPACKE_zgeev( LAPACK_ROW_MAJOR, 'N', 'V', bunches, lapack_m, bunches, lapack_eigenvalues, NULL, totalFilledBuckets, lapack_eigenvectors, totalFilledBuckets);
  } else {
    error_code = LAPACKE_zgeev( LAPACK_ROW_MAJOR, 'N', 'N', bunches, lapack_m, bunches, lapack_eigenvalues, NULL, totalFilledBuckets, NULL, totalFilledBuckets);
  }
  
  if (error_code != 0) {
    fprintf(stderr,"frequencyOfMaxTransGrowthRate: Subroutine zgeev returned error code %d.\n", error_code);
    for (iRes=0; iRes<dipoleHOMTypes; iRes++) {
      resonator = &dipoleHOMResonator[iRes];
      fprintf(stdout,"Resonator %ld\n",iRes);
      if (!dipoleHOMs_struct.Q&&(!resonator->Q||!resonator->RoQ))
        continue;
      if (resonator->cavities) 
        for (iCav=0; iCav<resonator->cavities; iCav++) {
          fprintf(stdout,"frequency[%ld] = %lf\n",iCav,resonator->frequency[iCav]);
        }
    }
    bombClinchor(NULL);
  }
  free(lapack_m);

  if ((doTransverseMotion_struct.eigenvectorFilename) && (eigenvectorFileWritten==0)) {
    eigenvectorFileWritten=1;
    /* normalize the eigenvectors */
#pragma omp parallel num_threads(threads) default(none) shared(lapack_eigenvectors) firstprivate(totalFilledBuckets) private(i,j,sum)
    {
#pragma omp for
      for (j=0; j<totalFilledBuckets; j++) {
        sum = 0;
        for (i=0; i<totalFilledBuckets; i++) {
          sum += sqr(lapack_eigenvectors[i*totalFilledBuckets+j].real) + sqr(lapack_eigenvectors[i*totalFilledBuckets+j].imag);
        }
        sum = sqrt(sum/totalFilledBuckets);
        for (i=0; i<totalFilledBuckets; i++) {
          lapack_eigenvectors[i*totalFilledBuckets+j].real /= sum;
          lapack_eigenvectors[i*totalFilledBuckets+j].imag /= sum;
        }
      }
    }
    if( !SDDS_InitializeOutput(&eigenvectorPage, SDDS_BINARY, 1,
                               "CBM eigenvectors", "CBM eigenvectors",
                               doTransverseMotion_struct.eigenvectorFilename ))
      SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
    if(0>SDDS_DefineColumn(&eigenvectorPage, "bunchIndex", NULL, NULL, NULL, NULL, SDDS_LONG, 0))
      SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
    for (i=0; i<bunches; i++) {
      sprintf(modeName, "mode%04ldReal", i);
      if(0>SDDS_DefineColumn(&eigenvectorPage, modeName, NULL, NULL, NULL, NULL, SDDS_DOUBLE, 0))
        SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
      sprintf(modeName, "mode%04ldImag", i);
      if(0>SDDS_DefineColumn(&eigenvectorPage, modeName, NULL, NULL, NULL, NULL, SDDS_DOUBLE, 0))
        SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
    }
    if( !SDDS_WriteLayout(&eigenvectorPage) ||
        !SDDS_StartPage(&eigenvectorPage, bunches))
      SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
    for (i=0; i<bunches; i++) {
      if (!SDDS_SetRowValues(&eigenvectorPage, SDDS_BY_NAME|SDDS_PASS_BY_VALUE, i,
                             "bunchIndex", i, NULL))
        SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
    }

    /* to output a matrix one has to write columns of the transpose! */
    eigenvectors_real = (double*) malloc(sizeof(double)*totalFilledBuckets);
    eigenvectors_imag = (double*) malloc(sizeof(double)*totalFilledBuckets); 
    for (i=0; i<bunches; i++) {
      for (j=0; j<totalFilledBuckets; j++) {
        eigenvectors_real[j] = lapack_eigenvectors[j*bunches+i].real;
        eigenvectors_imag[j] = lapack_eigenvectors[j*bunches+i].imag;
      }
      sprintf(modeName, "mode%04ldReal", i);
      if (!SDDS_SetColumn(&eigenvectorPage, SDDS_SET_BY_NAME|SDDS_PASS_BY_REFERENCE,
                          eigenvectors_real /*eigenvectorsT->ar[i]*/, bunches, modeName))
        SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
      sprintf(modeName, "mode%04ldImag", i);
      if (!SDDS_SetColumn(&eigenvectorPage, SDDS_SET_BY_NAME|SDDS_PASS_BY_REFERENCE,
                          eigenvectors_imag /*eigenvectorsT->ai[i]*/, bunches, modeName))
        SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
    }
    if (!SDDS_WritePage(&eigenvectorPage) || !SDDS_Terminate(&eigenvectorPage))
      SDDS_PrintErrors(stderr, SDDS_EXIT_PrintErrors|SDDS_VERBOSE_PrintErrors);
    free(eigenvectors_real);
    free(eigenvectors_imag);
    free(lapack_eigenvectors);
  }
  free(formFactor);
  /* find the eigenvalue representing the fastest growing CBM */
  maxGrowthRate=-1e38;
  iMax = 0;
  for (iCBM=0; iCBM<bunches; iCBM++) {
    if (maxGrowthRate<lapack_eigenvalues[iCBM].imag) {
      maxGrowthRate=lapack_eigenvalues[iCBM].imag;
      iMax=iCBM;  
    }
  }
  result.r = lapack_eigenvalues[iMax].real;
  result.i = lapack_eigenvalues[iMax].imag;
  return(result);
}

void loadDataFromTwissFile() 
/* load ringParameters data from the Twiss file, but only if there is not
 * a value already given. Column expected at Charge (nC) or Current (mA) and BunchLength (ps)
 */
{
  SDDS_DATASET SDDSin;
#define TWISSFILE_PARAMTERS 9
  char *parameterName[TWISSFILE_PARAMTERS] = {
    "pCentral", "U0", "alphac", "Sdelta0", "taudelta", "taux", "tauy", "nux", "nuy"
  };
  char *namelistQuantityName[9] = {
    "energyGeV", "energyLossPerTurnMeV", "momentumCompaction", "relativeEnergySpread", 
    "longDampingTime", "horTransDampingTime", "vertTransDampingTime", "horizontalTune", "verticalTune"
  };
  /* this is used to convert units */
  double factor[TWISSFILE_PARAMTERS] = {
    0.511e-3, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0
  };
  double longDampingTime=0, horTransDampingTime=0, vertTransDampingTime=0;;
  /* Uses data from file only if rate quantity is not specified in clinchor namelist. */
  if (longDampingRate!=0.0)
    longDampingTime = 1.0/longDampingRate;
  if (horTransDampingRate!=0.0 && longDampingTime)
    horTransDampingTime = 1.0/horTransDampingRate;
  if (vertTransDampingRate!=0.0 && longDampingTime)
    vertTransDampingTime = 1.0/vertTransDampingRate;
  double *target[TWISSFILE_PARAMTERS] = {
    &energyGeV, &energyLossPerTurnMeV, &momentumCompaction, &relativeEnergySpread, 
    &longDampingTime, &horTransDampingTime, &vertTransDampingTime, &horizontalTune, &verticalTune 
  };
  long i, rows, nRf;
  double betaxSum, betaySum;
  double *sData, *betax, *betay;
  char **elementType;
  
  if (!SDDS_InitializeInput(&SDDSin, twissFile) || SDDS_ReadPage(&SDDSin)!=1) {
    SDDS_SetError("Problem reading twiss file");
    SDDS_PrintErrors(stderr, SDDS_EXIT_PrintErrors|SDDS_VERBOSE_PrintErrors);
  }
  for (i=0; i<TWISSFILE_PARAMTERS; i++) {
    /* reads from file only if quantity is not specified in clinchor namelist. */
    if (*(target[i])==0) {
      if (!SDDS_GetParameterAsDouble(&SDDSin, parameterName[i], target[i])) {
        SDDS_SetError("Problem getting data from twiss file---are radiation-integral-related quantities included?");
        SDDS_PrintErrors(stderr, SDDS_EXIT_PrintErrors|SDDS_VERBOSE_PrintErrors);
      }
      *(target[i]) *= factor[i];
      fprintf(stdout, "Value of %s taken from file: %e\n", namelistQuantityName[i], *(target[i]));
    }
  }
  /* File actually stores damping time (taux) not rate. clinchor needs rate. */
  /* Uses data from file only if rate quantity is not specified in clinchor namelist. */
  if (longDampingRate==0.0 && longDampingTime)
    longDampingRate = 1.0/longDampingTime;
  if (horTransDampingRate==0.0 && longDampingTime)
    horTransDampingRate = 1.0/horTransDampingTime;
  if (vertTransDampingRate==0.0 && longDampingTime)
    vertTransDampingRate = 1.0/vertTransDampingTime;

  if (!(rows=SDDS_RowCount(&SDDSin)) || !(sData = SDDS_GetColumnInDoubles(&SDDSin, "s"))
      || !(betax = SDDS_GetColumnInDoubles(&SDDSin, "betax"))
      || !(betay = SDDS_GetColumnInDoubles(&SDDSin, "betay"))
      || !(elementType = SDDS_GetColumn(&SDDSin, "ElementType"))
      ) {
    SDDS_SetError("Problem getting column data from twiss file.");
    SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors);
    exit(1);
  }
  if (!circumference)
    circumference = sData[rows-1];
  free(sData);
  nRf=0;
  betaxSum=betaySum=0;
  for (i=0; i<rows; i++) {
    if (strcmp(elementType[i], "RFCA")==0 || 
        strcmp(elementType[i], "RFCW")==0) {
      nRf ++;
      betaxSum += betax[i];
      betaySum += betay[i];
    }
  }
  if (!betaxAtRFCavities) {
    if (!nRf)
      SDDS_Bomb("No RFCA or RFCW elements found in twiss file---you should supply betaxAtRFCavities and betayAtRFCavities");
    betaxAtRFCavities = betaxSum/nRf;
    fprintf(stdout, "Average value for betaxAtRFCavities from file: %e\n", betaxAtRFCavities);
  }
  if (!betayAtRFCavities) {
    if (!nRf)
      SDDS_Bomb("No RFCA or RFCW elements found in twiss file---you should supply betaxAtRFCavities and betayAtRFCavities");
    betayAtRFCavities = betaySum/nRf;
    fprintf(stdout, "Average value for betayAtRFCavities from file: %e\n", betayAtRFCavities);
  }
  
  SDDS_Terminate(&SDDSin);
}

void loadDataFromBunchLengthFile() 
/* load bunch length data from external file, which gives  
   more realistic bunch length factor which reduces
   effect of HOM of higher frequency. */

{
  SDDS_TABLE SDDSin;
  char *BunchCurrentColumnName, *BunchChargeColumnName, *BunchLengthColumnName;
  char *filename;
  
  BunchCurrentColumnName = BunchChargeColumnName = BunchLengthColumnName = NULL;
  filename = bunchLengthTableFile;

  if (!SDDS_InitializeInput(&SDDSin, filename) || SDDS_ReadPage(&SDDSin)!=1) {
    SDDS_SetError("Problem reading bunch length file");
    SDDS_PrintErrors(stderr, SDDS_EXIT_PrintErrors|SDDS_VERBOSE_PrintErrors);
  }
  bunchLengthValues = SDDS_CountRowsOfInterest(&SDDSin);
  BunchCurrentColumnName=SDDS_FindColumn(&SDDSin, FIND_NUMERIC_TYPE,
                                         "Current", "BunchCurrent", NULL);
  if (BunchCurrentColumnName){
    switch(SDDS_CheckColumn(&SDDSin, BunchCurrentColumnName, "mA", SDDS_DOUBLE, stderr)){
    case SDDS_CHECK_OKAY:
      break;
    case SDDS_CHECK_WRONGUNITS:
      switch(SDDS_CheckColumn(&SDDSin, BunchCurrentColumnName, NULL, SDDS_DOUBLE, stderr)){
      case SDDS_CHECK_OKAY:
        /* it is ok to forget to put in the units in the SDDS file. The units will be
           assumed to be mA */
        break;
      case SDDS_CHECK_WRONGUNITS:
        fprintf(stderr,"Invalid value for units of column %s in file %s. Permitted values are mA.\n", BunchCurrentColumnName, filename);
        bombClinchor(NULL);
        break;
      }
      break;
    default:
      printf("Something wrong with column %s in file %s.\n", BunchCurrentColumnName, filename);
      bombClinchor(NULL);
    }
  }
  BunchChargeColumnName=SDDS_FindColumn(&SDDSin, FIND_NUMERIC_TYPE,
                                        "Charge", "BunchCharge", NULL);
  if (BunchChargeColumnName){
    switch(SDDS_CheckColumn(&SDDSin, BunchChargeColumnName, "nC", SDDS_DOUBLE, stderr)){
    case SDDS_CHECK_OKAY:
      break;
    case SDDS_CHECK_WRONGUNITS:
      switch(SDDS_CheckColumn(&SDDSin, BunchChargeColumnName, NULL, SDDS_DOUBLE, stderr)){
      case SDDS_CHECK_OKAY:
        /* it is ok to forget to put in the units in the SDDS file. The units will be
           assumed to be nC */
        break;
      case SDDS_CHECK_WRONGUNITS:
        fprintf(stderr,"Invalid value for units of column %s in file %s. Permitted values are nC.\n", BunchChargeColumnName, filename);
        bombClinchor(NULL);
        break;
      }
      break;
    default:
      printf("Something wrong with column %s in file %s.\n", BunchChargeColumnName, filename);
      bombClinchor(NULL);
    }
  }
  BunchLengthColumnName=SDDS_FindColumn(&SDDSin, FIND_NUMERIC_TYPE,
                                        "Length", "BunchLength", NULL);
  if (BunchLengthColumnName){
    switch(SDDS_CheckColumn(&SDDSin, BunchLengthColumnName, "ps", SDDS_DOUBLE, stderr)){
    case SDDS_CHECK_OKAY:
      break;
    case SDDS_CHECK_WRONGUNITS:
      switch(SDDS_CheckColumn(&SDDSin, BunchLengthColumnName, NULL, SDDS_DOUBLE, stderr)){
      case SDDS_CHECK_OKAY:
        /* it is ok to forget to put in the units in the SDDS file. The units will be
           assumed to be ps */
        break;
      case SDDS_CHECK_WRONGUNITS:
        fprintf(stderr,"Invalid value for units of column %s in file %s. Permitted values are ps.\n", BunchLengthColumnName, filename);
        bombClinchor(NULL);
        break;
      }
      break;
    default:
      printf("Something wrong with column %s in file %s.\n", BunchLengthColumnName, filename);
      bombClinchor(NULL);
    }
  }

  /* if both columns for current current and charge exist, then current column will 
     take precedence. */
  if (!BunchCurrentColumnName && !BunchChargeColumnName) {
    printf("No valid columns for bunch current or bunch charge was given in file %s.\n", filename);
    bombClinchor(NULL);
  }
  if (!BunchLengthColumnName) {
    printf("No valid column for bunch length was given in file %s.\n", filename);
    bombClinchor(NULL);
  }
  
  /* double arrays bunchCurrent bunchCharge bunchLength are global */
  bunchLengthPicoSecondData = SDDS_GetColumnInDoubles(&SDDSin, BunchLengthColumnName);
  if (BunchCurrentColumnName)
    bunchCurrentData = SDDS_GetColumnInDoubles(&SDDSin, BunchCurrentColumnName);
  else
    bunchChargeData = SDDS_GetColumnInDoubles(&SDDSin, BunchChargeColumnName);
  
  SDDS_Terminate(&SDDSin);
}

void do_semaphore_setup(char **semaphoreFile, NAMELIST_TEXT *nltext, char *rootname)
{
  
  set_namelist_processing_flags(STICKY_NAMELIST_DEFAULTS);
  set_print_namelist_flags(0);
  if (processNamelist(&semaphores, nltext)==NAMELIST_ERROR)
    bomb("Problem processing semaphores command", NULL);
  print_namelist(stdout, &semaphores);
 
  if (started)
    SDDS_CopyString(&semaphoreFile[0], started);
  if (done)
    SDDS_CopyString(&semaphoreFile[1], done);
  if (failed)
    SDDS_CopyString(&semaphoreFile[2], failed);

  /* usually one would expect values for strings %s.started, %s.done, %s.failed
   */
  if (rootname) {
    semaphoreFile[0] = compose_filename(semaphoreFile[0], rootname);
    semaphoreFile[1] = compose_filename(semaphoreFile[1], rootname);
    semaphoreFile[2] = compose_filename(semaphoreFile[2], rootname);
  }
}

void createSemaphoreFile(char *filename)
{
  FILE *fp;
  if (filename) {
    fprintf(stdout, "Creating semaphore file %s\n", filename);
  } else 
    return;
  if (!(fp = fopen(filename, "w"))) {
    fprintf(stdout, "Problem creating semaphore file %s\n", filename);
    bombClinchor(NULL);
  } else { /* Put the CPU time in the file */
    fprintf(fp, "%8.2f\n", cpu_time()/100.0); 
  }
  fclose(fp);
}

void bombClinchor(char *error)
{
  if (error)
    fprintf(stdout, "error: %s\n", error);
  if (semaphoreFile[2]) 
    createSemaphoreFile(semaphoreFile[2]);
  exit(1);
}


char *compose_filename(char *template, char *root_name)
{
  char *ptr;

  if (str_in(template, "%s")) {
    ptr = tmalloc(sizeof(char)*(unsigned long)(strlen(template)+strlen(root_name)+1));
    sprintf(ptr, template, root_name);
    return(ptr);
  }
  else
    return(template);
}
