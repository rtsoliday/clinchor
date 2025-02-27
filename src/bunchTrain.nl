/*************************************************************************\
* Copyright (c) 2002 The University of Chicago, as Operator of Argonne
* National Laboratory.
* Copyright (c) 2002 The Regents of the University of California, as
* Operator of Los Alamos National Laboratory.
* This file is distributed subject to a Software License Agreement found
* in the file LICENSE that is included with this distribution. 
\*************************************************************************/

#namelist bunchTrain static
long startBucket = 0;
long bucketInterval = 1;
long bunches = 1;
double currentPerBucketMA = 0.0;
double totalCurrentMA = 0.0;
long clearPreviousPatterns = 0;
#end
 
