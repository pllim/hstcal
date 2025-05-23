# include <stdio.h>
# include <string.h>        /* for strncmp, strcmp */

#include "hstcal.h"
# include "hstio.h"
# include "acs.h"
# include "acsinfo.h"
# include "hstcalerr.h"        /* defines error codes */

static int checkAtoD (Hdr *, ACSInfo *, int *, int *);
static int checkBias (Hdr *, ACSInfo *, int *, int *);
static int checkBlev (Hdr *, ACSInfo *, int *, int *);
static int checkCCD (Hdr *, ACSInfo *, int *);
static int checkDQI (Hdr *, ACSInfo *, int *, int *);
static int checkSink (Hdr *, ACSInfo *, int *, int *);


/* This routine gets the names of reference images and tables from the
 primary header and checks for dummy pedigree.

 Warren Hack, 1998 June 8:
 Original ACS version based on Phil Hodge's CALSTIS routine...
 **
 **    No Major modifications from STIS code...

 Warren Hack, 2001 Feb 7:
 Finished revisions for supporting Post-Flash
 processing.
 **

 Pey Lian Lim, 2012 Dec 12:
 Moved FLSHCORR to ACS2D.
 **

 Pey Lian Lim, 2012 Dec 19:
 Make sure check is skipped if input *acs flag is not PERFORM.
 Otherwise, acsccd.e flags do not work properly.
 **

 Pey Lian Lim, 2013 Aug 12:
 Separated PCTECORR from ACSCCD.
 **

 Pey Lian Lim, 2017 Feb 21:
 Added new SINKCORR.
 **

 Michele De La Pena, 2020 May 18:
 Added new SATUFILE: Full-well saturation image.
*/
int GetACSFlags (ACSInfo *acs, Hdr *phdr) {

    extern int status;

    int missing = 0;    /* true if any calibration file is missing */
    int nsteps = 0;     /* number of calibration steps to perform */

    int GetccdSw (ACSInfo *, Hdr *);

    /* Get the values for the Calibration Switches from the
    **    header for processing.
    */
    if (GetccdSw (acs, phdr) )
        return(status);

    /* Although this should never be called with MAMA data, we
       just want to be safe...
    */
    if (acs->detector != MAMA_DETECTOR) {
        if (checkCCD (phdr, acs, &missing))
            return (status);
    }

    /* Check each reference file that we need. */

    if (checkDQI (phdr, acs, &missing, &nsteps))
        return (status);

    if (checkAtoD (phdr, acs, &missing, &nsteps))
        return (status);

    if (checkBlev (phdr, acs, &missing, &nsteps))    /* no reference file */
        return (status);

    if (checkBias (phdr, acs, &missing, &nsteps))
        return (status);

    if (checkSink (phdr, acs, &missing, &nsteps))
        return (status);

    if (missing) {
        return (status = CAL_FILE_MISSING);
    } else if (nsteps < 1) {
        trlwarn ("No calibration switch was set to PERFORM,\n"
                 "            or all reference files had PEDIGREE = DUMMY.");
        return (status = NOTHING_TO_DO);
    } else {
        return (status);
    }
}


/* If this step is to be performed, check for the existence of the
   atod table.
*/
static int checkAtoD (Hdr *phdr, ACSInfo *acs, int *missing, int *nsteps) {

    /* arguments:
       Hdr *phdr        i: primary header
       ACSInfo *acs   i: switches, file names, etc
       int *missing     io: incremented if the file is missing
       int *nsteps      io: incremented if this step can be performed
    */

    extern int status;

    int calswitch;
    int GetSwitch (Hdr *, char *, int *);
    int GetTabRef (RefFileInfo *, Hdr *, char *, RefTab *, int *);
    void MissingFile (char *, char *, int *);

    /* Are we supposed to do this step? */
    if (acs->atodcorr == PERFORM) {

        if (GetSwitch (phdr, "ATODCORR", &calswitch))
            return (status);
        if (calswitch == COMPLETE) {
            acs->atodcorr = OMIT;
            return (status);
        }

        /* Get the table name, check that the file exists, and get
           pedigree and descrip.
        */
        if (GetTabRef (acs->refnames, phdr,
                       "ATODTAB", &acs->atod, &acs->atodcorr))
            return (status);
        if (acs->atod.exists != EXISTS_YES)
            MissingFile ("ATODTAB", acs->atod.name, missing);

        if (acs->atodcorr == PERFORM)
            (*nsteps)++;
    }

    return (status);
}


/* If this step is to be performed, check for the existence of the
   bias file.  If it exists, get the pedigree and descrip keyword values.
*/
static int checkBias (Hdr *phdr, ACSInfo *acs, int *missing, int *nsteps) {

  /* arguments:
   Hdr *phdr        i: primary header
   ACSInfo *acs   i: switches, file names, etc
   int *missing     io: incremented if the file is missing
   int *nsteps      io: incremented if this step can be performed
   */

    extern int status;

    int saveBiasCorr = GOOD_PEDIGREE;

    int calswitch;
    int GetSwitch (Hdr *, char *, int *);
    int GetImageRef (RefFileInfo *, Hdr *, char *, RefImage *, int *);
    void MissingFile (char *, char *, int *);

    /* Are we supposed to do this step? */
    if (acs->biascorr == PERFORM) {

        if (GetSwitch (phdr, "BIASCORR", &calswitch))
            return (status);
        if (calswitch == COMPLETE) {
            acs->biascorr = OMIT;
            return (status);
        }

        if (GetImageRef (acs->refnames, phdr,
                         "BIASFILE", &acs->bias, &acs->biascorr))
            return (status);

        if (acs->bias.exists != EXISTS_YES)
            MissingFile ("BIASFILE", acs->bias.name, missing);
        if (acs->biascorr == PERFORM)
            (*nsteps)++;

        /* Save the value for recovery */
        saveBiasCorr = acs->biascorr;

        /*
          Also check for the new full-well saturation image which is
          applied after BIASCORR, conversion to elections, and BLEVCORR
          are done. Since the reference file is not associated with its
          own "calibration step keyword" (e.g., SATUCORR), just using the
          BIASCORR key as a standin here - make sure the BIASCORR retains
          its value as set in the above code.

          This is a kludge.
       */
        if (GetImageRef (acs->refnames, phdr,
                         "SATUFILE", &acs->satmap, &acs->biascorr))
            return (status);

        /* Recover the biascorr setting */
        acs->biascorr = saveBiasCorr;

        if (acs->satmap.exists != EXISTS_YES)
            MissingFile ("SATUFILE", acs->satmap.name, missing);
    }

    return (status);
}


static int checkBlev (Hdr *phdr, ACSInfo *acs, int *missing, int *nsteps) {

    /* arguments:
       Hdr *phdr        i: primary header
       ACSInfo *acs   i: switches, file names, etc
       int *nsteps      io: incremented if this step can be performed
    */

    extern int status;

    int calswitch;
    int GetSwitch (Hdr *, char *, int *);
    int GetTabRef (RefFileInfo *, Hdr *, char *, RefTab *, int *);
    void MissingFile (char *, char *, int *);

    /* Are we supposed to do this step? */
    if (acs->blevcorr == PERFORM) {

        if (GetSwitch (phdr, "BLEVCORR", &calswitch))
            return (status);

        if (calswitch == COMPLETE) {
            acs->blevcorr = OMIT;
            return (status);
        }

        if (acs->blevcorr == PERFORM)
            (*nsteps)++;
    }

    return (status);
}


/* If the detector is the CCD, we need the CCD parameters table for
   BIASCORR, DARKCORR, and PHOTCORR.  This routine checks that the table
   exists.

   We also need the table for initializing the error array, but we
   don't have a flag for that step.  That's why we need this table
   regardless of which steps are to be performed.
*/
static int checkCCD (Hdr *phdr, ACSInfo *acs, int *missing) {

    /* arguments:
       Hdr *phdr        i: primary header
       ACSInfo *acs   i: switches, file names, etc
       int *missing     io: incremented if the table is missing
    */

    extern int status;
    int calswitch;            /* returned by GetTabRef and ignored */
    int GetTabRef (RefFileInfo *, Hdr *, char *, RefTab *, int *);
    void MissingFile (char *, char *, int *);
    int GotFileName (char *);

    if (acs->detector == MAMA_DETECTOR)
        return (status);

    if (GetTabRef (acs->refnames, phdr,
                   "CCDTAB", &acs->ccdpar, &calswitch))
        return (status);

    if (acs->ccdpar.exists != EXISTS_YES) {

        MissingFile ("CCDTAB", acs->ccdpar.name, missing);

    } else if (acs->ccdpar.goodPedigree != GOOD_PEDIGREE) {

        (*missing)++;
        trlerror("CCDTAB `%s' is a dummy table.", acs->ccdpar.name);
    }

    /* Get OSCNTAB here as it applies to all CCD processing as well.
       WJH 6 May 1999
    */
    if (GetTabRef (acs->refnames, phdr,
                   "OSCNTAB", &acs->oscn, &acs->blevcorr))
        return (status);

    if (acs->oscn.exists != EXISTS_YES) {
        if (GotFileName (acs->oscn.name)) {
            MissingFile ("OSCNTAB", acs->oscn.name, missing);
        }
    }

    return (status);
}


/* Check whether we should assign initial values to the data quality
   array.  There is a reference table but not an image for this step.

   For the CCD, there are two steps to DQICORR, checking for saturation
   and using the BPIXTAB to initialize the data quality array.  If no
   bad pixel table was specified (name is blank or "N/A"), or if the
   table is dummy, we can still do this calibration step, but it will
   just consist of checking and flagging saturation.  For the MAMAs,
   however, if this switch is set to PERFORM, the table must exist.

   Note that, unlike most other switches, dqicorr will not be reset to
   omit if the header value is "COMPLETE", since it would not cause any
   problem to perform this step more than once.  The user might do so
   deliberately in order to accumulate the flags from more than one table.
*/
static int checkDQI (Hdr *phdr, ACSInfo *acs, int *missing, int *nsteps) {

    /* arguments:
       Hdr *phdr        i: primary header
       ACSInfo *acs   i: switches, file names, etc
       int *missing     io: incremented if the table is missing
       int *nsteps      io: incremented if this step can be performed
    */

    extern int status;

    int GotFileName (char *);
    int GetTabRef (RefFileInfo *, Hdr *, char *, RefTab *, int *);
    void MissingFile (char *, char *, int *);

    if (acs->dqicorr == PERFORM) {

        if (GetTabRef (acs->refnames, phdr,
                       "BPIXTAB", &acs->bpix, &acs->dqicorr))
            return (status);

        if (acs->bpix.exists != EXISTS_YES) {

            if (acs->detector == MAMA_DETECTOR ||
                GotFileName (acs->bpix.name)) {

                MissingFile ("BPIXTAB", acs->bpix.name, missing);
            }
        }

        if (acs->dqicorr == PERFORM)
            (*nsteps)++;
    }

    return (status);
}


/* If this step is to be performed, check for the existence of the
   sink file.  If it exists, get the pedigree and descrip keyword values.
*/
static int checkSink (Hdr *phdr, ACSInfo *acs, int *missing, int *nsteps) {
    /* arguments:
       Hdr *phdr        i: primary header
       ACSInfo *acs     i: switches, file names, etc
       int *missing     io: incremented if the file is missing
       int *nsteps      io: incremented if this step can be performed
    */
    extern int status;

    int calswitch;
    int GetSwitch (Hdr *, char *, int *);
    int GetImageRef (RefFileInfo *, Hdr *, char *, RefImage *, int *);
    void MissingFile (char *, char *, int *);

    /* Are we supposed to do this step? */
    if (acs->sinkcorr == PERFORM) {

        if (GetSwitch (phdr, "SINKCORR", &calswitch))
            return (status);
        if (calswitch == COMPLETE) {
            acs->sinkcorr = OMIT;
            return (status);
        }


        if (GetImageRef (acs->refnames, phdr,
                         "SNKCFILE", &acs->sink, &acs->sinkcorr))
            return (status);

        if (acs->sink.exists != EXISTS_YES)
            MissingFile ("SNKCFILE", acs->sink.name, missing);
        if (acs->sinkcorr == PERFORM)
            (*nsteps)++;
    }

    return (status);
}
