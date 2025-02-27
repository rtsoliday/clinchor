/*************************************************************************\
* Copyright (c) 2002 The University of Chicago, as Operator of Argonne
* National Laboratory.
* Copyright (c) 2002 The Regents of the University of California, as
* Operator of Los Alamos National Laboratory.
* This file is distributed subject to a Software License Agreement found
* in the file LICENSE that is included with this distribution. 
\*************************************************************************/

#namelist ringParameters static
STRING twissFile = NULL;
double energyGeV = 0.0;
double circumference = 0.0;
double energyLossPerTurnMeV = 0.0;
double rfVoltageMV = 0.0;
long harmonicNumber = 1;
double momentumCompaction = 0.0;
double relativeEnergySpread = 0.0;
double bunchLengtheningFactor = 1.0;
STRING bunchLengthTableFile = NULL;
STRING bunchLengthUsedFile = NULL;
double longDampingRate = 0.0;
double transDampingRate = 0.0;
double horTransDampingRate = 0.0;
double vertTransDampingRate = 0.0;
double horizontalTune = 0;
double verticalTune = 0;
double betaxAtRFCavities = 0;
double betayAtRFCavities = 0;
STRING parameters = NULL;
#end

#namelist monopoleHOMs struct
STRING filename = NULL;
long clearPreviousMonopoleHOMs = 0;
double deQFactorMultiplier = 1 ;
double Q = 0;
double staggeringStepMultiplier = 1;
long detuneFundamental = 0;
long keepFundamentalFixed = 0;
#end

#namelist doLongitudinalMotion struct
long normalModes = 1;
long doLaplace = 0;
STRING eigenvectorFilename = NULL;
STRING CBMFrequencyFilename = NULL;
#end

#namelist dipoleHOMs struct
STRING filename = NULL;
long clearPreviousDipoleHOMs = 0;
double deQFactorMultiplier = 1;
double Q = 0;
double staggeringStepMultiplier = 1;
long detuneFundamental = 0;
#end

#namelist doTransverseMotion struct
long normalModes = 1;
long doLaplace = 0;
long verticalDirection = 0;
long horizontalDirection = 0;
STRING CBMFrequencyFilename = NULL;
STRING eigenvectorFilename = NULL;
#end

#namelist sweepHOMFrequency struct
long resonatorIndex = 0;
double frequencyRange = 0.0;
long points = 100;
STRING filename = NULL;
#end

#namelist randomizeHOMFrequencies struct
double spread = 0.0;
long seed = -987654321;
long uniform = 1;
long samples = 100;
STRING CBMFrequencyFilename = NULL;
STRING HOMFrequencyFilename = NULL;
#end

/* Only does uniform distributions (not gaussian) for now to simulate variation
from topup operations */
#namelist randomizeBunchCurrent struct
long seed = -987654323;
long uniform = 1;
double absoluteSpreadMA = 0;
double relativeSpread = 0;
#end

#namelist semaphores
    STRING started = "%s.started";
    STRING done = "%s.done";
    STRING failed = "%s.failed";
#end
