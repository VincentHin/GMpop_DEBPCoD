/*
 NAME
 GMpop_DEBPCoD.c
 
 NOTES
 The problem specific definition file for the model describing the dynamics
 of a population of marine mammals that is feeding on a shared resource.
 When ENVIRON_DIM is set equal to 1 in GMpop_DEBPCoD.h, the mammals are assumed to grow
 exponentially, experiencing no density dependence.
 When ENVIRON_DIM is set equal to 2 in GMpop_DEBPCoD.h, the mammals are assumed to compete
 for the limiting resource, which is assumed to be replenished with a constant 
 productivity and turn-over rate, resulting in semi-chemostat growth.
 
 - The model tracks individual mammals. Reproduction and death are discrete
 events, while all other processes occur in continuous time.
 - Reproduction events occur, however, continuously throughout the year
 - Pregnancy starts occurs when females have accumulated sufficient reserves
 to cover all additional costs for pregnancy.
 
 HISTORY
 VH - Dec 3, 2019: Last Revision.
 */
#include "escbox.h"

#ifndef DEBUG
#define DEBUG                     0
#endif

// Rescaling is only useful in case of the density-independent model version
#if (ENVIRON_DIM == 2)
#define RESCALEINITPOP            0                                                 // Rescale initial population to RESCALEPOP
#define RESCALEPERIOD             -1                                                // Rescale population every so many years, use -1 for no rescaling
#else
#define RESCALEINITPOP            1                                               	// Rescale initial population to RESCALEPOP
#define RESCALEPERIOD             200                                               // Rescale population every so many years, use -1 for no rescaling
#define RESCALEPOP                100	                                            // If given in ISF file, scale initial population back to this number
#endif

#define INITPOP                   50                                                // Initial seeding population
#define	FIXED_PREGDELAY			  0													// Fix the pregnancy delay period to 445 days


/*
 *====================================================================================================================================
 *
 *  LABELLING ENVIRONMENT AND I-STATE VARIABLES
 *
 *====================================================================================================================================
 */
#define time                      (env[0])
#if (ENVIRON_DIM == 2)
#define resource                  (env[1])
#else
#define resource                  (RESOURCE)
#endif

//		number										// Column 1
#define age                       (i_state( 0))		// Column 2
#define reserves                  (i_state( 1))
#define totalmort                 (i_state( 2))

#define IDmale                    (i_const( 0))		// Column 5
#define IDfemale                  (i_const( 1))
#define IDwaiting				          (i_const( 2))		    					// Waiting for pregnancy (0 or 1)
#define IDpregnant                (i_const( 3))		// Column 8
#define IDlactating               (i_const( 4))
#define IDindividualID            (i_const( 5))
#define IDmotherID                (i_const( 6))		// Column 11
#define IDmotherIndex             (i_const( 7))

#define IDminloglifespan          (i_const( 8))		// Column 13
#define IDpregwait				        (i_const( 9))										// Future age at which female becomes pregnant
#define IDinseminated             (i_const(10))                   // Age at start of last pregnancy
#define IDcalved                  (i_const(11))   // Column 16    // Age of last calving
#define IDweaned                  (i_const(12))					          // Age of last successful weaning

#define IDtotalfemcalves          (i_const(13))                   // Total number of female calves produced during life (Using for R0 calculation)
#define IDtotalcalves          	  (i_const(14))										// Total number of calves produced during life
#define IDtotalweaned             (i_const(15))   // Column 20    // Total number of calves weaned during life
#define IDagefirstreceptive       (i_const(16))										// Age of becoming receptive for first time
#define IDagefirstbirth	          (i_const(17))										// Age at birth of first calf
#define IDagefirstweaning         (i_const(18))										// Age at weaning of first calf

#define IDlength                  (i_const(19))		// Column 24
#define IDbones                   (i_const(20))
#define IDweight                  (i_const(21))
#define IDweightM                 (i_const(22))		// Column 27
#define IDfatratio                (i_const(23))

#define IDingest                  (i_const(24))
#define IDringest                 (i_const(25))		// Column 30
#define IDmingest                 (i_const(26))

#define IDmaint                   (i_const(27))
#define IDgrowth                  (i_const(28))		// Column 33
#define IDpregcosts               (i_const(29))
#define IDlactcosts               (i_const(30))
#define IDnet_energy              (i_const(31))		// Column 36

#define IDmortality               (i_const(32))
#define IDbackground              (i_const(33))
#define IDstarvation              (i_const(34))		// Column 39
#define IDstarvdays	              (i_const(35))

#define IDfeedinglevel            (i_const(36))		// Column 41

#define ismale(n)                 (popIDcard[MAMMALS][(n)][IDmale]      > TINY)
#define isfemale(n)               (popIDcard[MAMMALS][(n)][IDfemale]    > TINY)
#define iscalve(n)                (pop[MAMMALS][(n)][age]               < Tl)
#define iswaiting(n)              (popIDcard[MAMMALS][(n)][IDwaiting]   > TINY)
#define ispregnant(n)             (popIDcard[MAMMALS][(n)][IDpregnant]  > TINY)
#define islactating(n)            (popIDcard[MAMMALS][(n)][IDlactating] > TINY)

#define isdead(n)                 (pop[MAMMALS][(n)][number] < 0.0)
#define isdying(n)                ((popIDcard[MAMMALS][(n)][IDfatratio] < MINFATRATIO) || (pop[MAMMALS][(n)][totalmort] >= popIDcard[MAMMALS][n][IDminloglifespan]) || (pop[MAMMALS][(n)][age] > 21900))

#define StdDev(ss, n)             sqrt((ss)/((double)(n - 1)))
#define isone(x)                  ((x > 1.0) || iszero(x - 1.0))


/*
 *====================================================================================================================================
 *
 *  DEFINING AND LABELLING CONSTANTS AND PARAMETERS
 *
 *====================================================================================================================================
 */

#define MAMMALS                   0
#define YEAR                      365.0                                             // Total year length
#define VOLUME                    1.0E6                                             // Volume for consumer versus resource
#define MINFATRATIO               0.005                                             // Minimum fat level below individuals die with certainty
#define TINY                      1.0E-3

#define Tg                        parameter[ 0]                                     // Gestation period
#define Tl                        parameter[ 1]                                     // Lactation period
#define Tf                        parameter[ 2]                                     // Age at independent foraging
#define Tr                        parameter[ 3]                                     // Age at 100% foraging capacity

#define Lb                        parameter[ 4]                                     // Length at birth
#define Linf                      parameter[ 5]                                     // Asymptotic length females
#define Linf_m                    parameter[ 6]                                     // Asymptotic length males
#define K                         parameter[ 7]                                     // Von Bertalanffy growth rate females
#define K_m		                  parameter[ 8]                                     // Von Bertalanffy growth rate males

#define OmegaS                    parameter[ 9]                                     // Weight-length scalar
#define OmegaE                    parameter[ 10]                                     // Weight-length exponent
#define OmegaM                    parameter[ 11]                                     // Relative maintenance reserves
#define OmegaC                    parameter[12]                                     // Relative maintenance foetus

#define Rho                       parameter[13]                                     // Target reserve-weight ratio
#define RhoS                      parameter[14]                                     // Starvation reserve-weight ratio

#define EtaF                      parameter[15]                                     // Steepness reserves-assimilation response
#define Gamma                     parameter[16]                                     // Shape parameter of resource assimilation-age response

#define PhiM                      parameter[17]                                     // Energy provisioning ratio through lactation

#define XiM                       parameter[18]                                     // Non-linearity in mother reserves-lactation rate
#define XiF                       parameter[19]                                     // Non-linearity in resource assimilation - calve age rate
#define XiC                       parameter[20]                                     // Non-linearity in milk assimilation-calve age rate

#define SigmaM                    parameter[21]                                     // Field metabolic rate scalar
#define SigmaG                    parameter[22]                                     // Structural mass growth costs
#define SigmaL                    parameter[23]                                     // Reserve-milk conversion ratio

#define Alpha1                    parameter[24]                                     // Background mortality parameter female
#define Beta1					  parameter[25]                                     // Background mortality parameter female
#define Alpha2                    parameter[26]                                     // Background mortality parameter female
#define Beta2                     parameter[27]                                     // Background mortality parameter female
#define Mu_male                   parameter[28]                                     // Background mortality parameter male

#define PREGCHANCE				  parameter[29]                                     // Daily chance of becoming pregnant

#define MuS                       parameter[30]                                     // Starvation mortality scalar

#define EpsPlus                   parameter[31]                                     // Anabolic conversion efficiency
#define EpsMin                    parameter[32]                                     // Catabolic conversion efficiency

#define RESOURCE                  parameter[33]                                     // Resource density

#define RESOURCETURNOVER          parameter[34]                                     // Resource turn-over rate

#define UNIDEVSEED                parameter[35]                                     // Seed value random generator

#define RESOURCEAMPLITUDE         parameter[36]                                     // Amplitude of annual resource fluctuations
#define DISTURBANCESTART          parameter[37]                                     // Start date of annual disturbance
#define DISTURBANCEPERIOD         parameter[38]                                     // Duration of annual disturbance

/*
 *====================================================================================================================================
 *
 *  DEFINING GLOBAL VARIABLES AND USER-IMPLEMENTED FUNCTIONS
 *
 *====================================================================================================================================
 */

extern void                       SievePop(void);

static void                       UpdateIDcards(double *env, population *pop);
static void                       YearlyDeaths(double *env, population *pop, FILE *fpnt);

static void                       UpdateStats(double value, double *mean, double *sum_sq, double *minval, double *maxval, long n, int countzero);
static double                     bexp(double x);
static double                     UniDev(void);

static long unsigned              UniDevSeed        = 0L;
static long unsigned              individualIDcount = 0L;
static long unsigned              initialIDcount = 0L;
static double                     SizeAtBirth, TotalNeonateCosts;

static long unsigned              Observations = 0L;
static long unsigned              ObsAgeFirstReceptive = 0L;
static long unsigned              ObsAgeFirstBirth = 0L;
static long unsigned              ObsAgeFirstWeaning = 0L;
static long unsigned              ObsStarvDays = 0L;
static long unsigned              ObsIBI = 0L;
static long unsigned              ObsIWI = 0L;
static double                     meanR0, ssR0, minR0, maxR0;
static double                     meanTotalCalves, ssTotalCalves, minTotalCalves, maxTotalCalves;
static double                     meanTotalWeaned, ssTotalWeaned, minTotalWeaned, maxTotalWeaned;
static double					  meanAgeFirstReceptive, ssAgeFirstReceptive, minAgeFirstReceptive, maxAgeFirstReceptive;
static double					  meanAgeFirstBirth, ssAgeFirstBirth, minAgeFirstBirth, maxAgeFirstBirth;
static double					  meanAgeFirstWeaning, ssAgeFirstWeaning, minAgeFirstWeaning, maxAgeFirstWeaning;
static double					  meanStarvDays, ssStarvDays, minStarvDays, maxStarvDays;
static double					  meanReproFemales, ssReproFemales, minReproFemales, maxReproFemales;
static double					  meanStarvingFemales, ssStarvingFemales, minStarvingFemales, maxStarvingFemales;
static double                     meanIBI, ssIBI, minIBI, maxIBI;
static double					  meanIWI, ssIWI, minIWI, maxIWI;

static char                       msg[256];

static FILE                       *deathfile;
static FILE 					  *R0file;

#define SEXCLASSES	              8
#define AGECLASSES	              6

static int                        Deaths[SEXCLASSES][AGECLASSES];

/*
 *====================================================================================================================================
 *
 *  USER INITIALIZATION ROUTINE ALLOWS OPERATIONS ON INITIAL POPULATIONS
 *
 *====================================================================================================================================
 */

void UserInit(int argc, char **argv, double *env, population *pop)
  
{
  int     i, j, n;
  //	int		k;
  double  tmp;
  double  NeonateWeightCosts, NeonateMaintCosts, NeonateLactCosts;
  char	  fn[1024];
  
  switch (argc)
  {
  case 4:
    RESOURCE    = atof(argv[3]);
  case 3:
    UNIDEVSEED  = atof(argv[2]);
  default:
    break;
  }
  
  // Initialize the random number generator
  if (UNIDEVSEED > 0.0) UniDevSeed = (long)(UNIDEVSEED + TINY);
  UniDev();
  
  SizeAtBirth         = OmegaS*pow(Lb, OmegaE);
  
  NeonateWeightCosts  = SigmaG*SizeAtBirth;                                         // Structural mass costs in MJ
  NeonateWeightCosts /= EpsMin;                                                     // Convert into kg stored reserves
  NeonateWeightCosts += RhoS*SizeAtBirth/(1-RhoS);                                  // Add the reserves transferred at birth to the calve
  
  NeonateMaintCosts   = SigmaM*(4/(3*OmegaE+4))*pow(OmegaC*OmegaS*pow(Lb, OmegaE), 0.75)*Tg;      // Maintenance costs in MJ
  NeonateMaintCosts  /= EpsMin;                                                     // Convert into kg stored reserves
  
  NeonateLactCosts    = 0.0;                                                        // DO NOT KNOW YET HOW MUCH THIS IS
  NeonateLactCosts   /= EpsMin;                                                     // Convert into kg stored reserves
  
  TotalNeonateCosts   = NeonateWeightCosts;
  
  // Initialize the average statistics
  Observations = ObsAgeFirstReceptive = ObsAgeFirstBirth = ObsAgeFirstWeaning = ObsStarvDays = ObsIBI = ObsIWI = 0.0;
  
  meanR0 = ssR0 = minR0 = maxR0 = 0.0;
  meanTotalCalves = ssTotalCalves = minTotalCalves = maxTotalCalves = 0.0;
  meanTotalWeaned = ssTotalWeaned = minTotalWeaned = maxTotalWeaned = 0.0;
  meanAgeFirstReceptive = ssAgeFirstReceptive = minAgeFirstReceptive = maxAgeFirstReceptive = 0.0;
  meanAgeFirstBirth = ssAgeFirstBirth = minAgeFirstBirth = maxAgeFirstBirth = 0.0;
  meanAgeFirstWeaning = ssAgeFirstWeaning = minAgeFirstWeaning = maxAgeFirstWeaning = 0.0;
  meanStarvDays = ssStarvDays = minStarvDays = maxStarvDays = 0.0;
  meanReproFemales = ssReproFemales = minReproFemales = maxReproFemales = 0.0;
  meanStarvingFemales = ssStarvingFemales = minStarvingFemales = maxStarvingFemales = 0.0;
  meanIBI = ssIBI = minIBI = maxIBI = 0.0;
  meanIWI = ssIWI = minIWI = maxIWI = 0.0;
  
  if (cohort_no[MAMMALS] < 1)
  {
    // No ISF file found. Initialize with small seeding population
    time      = 0.0;
#if (ENVIRON_DIM == 2)
    resource  = RESOURCE;
#endif
    
    AddCohorts(pop, MAMMALS, INITPOP);
    for (i = 0; i < cohort_no[MAMMALS]; i++)
    {
      // Introduce individuals with a randong age between 5 and 25 years
      pop[MAMMALS][i][number]  = 1.0;
      pop[MAMMALS][i][age]     = (5.0 + 20.0*UniDev())*YEAR;
      
      // Set all IDcards to zero
      for (n = 0; n < I_CONST_DIM; n++) popIDcard[MAMMALS][i][n] = 0.0;
      
      // Male (n=1) or female (n=0) calve produced
      n = (UniDev() < 0.5);
      popIDcard[MAMMALS][i][IDmale]          = n;
      popIDcard[MAMMALS][i][IDfemale]        = 1.0 - n;
      
      // That have their minimum reserves/weight ratio
      if(isfemale(i)){
        tmp = Linf - (Linf - Lb)*bexp(-K*pop[MAMMALS][i][age]);
      } else {
        tmp = Linf_m - (Linf_m - Lb)*bexp(-K_m*pop[MAMMALS][i][age]);
      }
      pop[MAMMALS][i][reserves]  = OmegaS*pow(tmp, OmegaE);
      pop[MAMMALS][i][reserves] *= RhoS/(1.0 - RhoS);
      
      // And only experienced background mortality (should still be adjusted for males)
      pop[MAMMALS][i][totalmort] = -Alpha1*bexp(-Beta1 * pop[MAMMALS][i][age]) / Beta1 + Alpha2*bexp(Beta2 * pop[MAMMALS][i][age]) / Beta2 + Alpha1 / Beta1 - Alpha2 / Beta2;
      
      popIDcard[MAMMALS][i][IDindividualID]  = i;
      popIDcard[MAMMALS][i][IDmotherID]      = -1;
      popIDcard[MAMMALS][i][IDmotherIndex]   = -1;
      
      // Generate a random number between 0 and 1. Store this as minus the logarithm of this value
      // The individual will die of natural mortality if integrated mortality rate exceeds this random deviate.
      do
      {
        tmp = UniDev();
        tmp = max(tmp, DBL_MIN);
        tmp = min(tmp, 1.0 - DBL_MIN);
        tmp = -log(tmp);
      }
      while (tmp <= pop[MAMMALS][i][totalmort]);
      popIDcard[MAMMALS][i][IDminloglifespan] = tmp;
    }
  }
  else
  {
    // ISF file found. Assume initial population is fully specified
    
#if (RESCALEINITPOP)
    // Kill all individuals
    for (i = 0; i < cohort_no[MAMMALS]; i++) pop[MAMMALS][i][number]  = -1.0;
    
    int females = 0;
    // Resurrect a specific number of individuals at random (aselective sampling without replacement) until RESCALEPOP females remain
    for (i = 0; i < cohort_no[MAMMALS]; i++)
    {
      while (1)
      {
        j = (int)floor(UniDev()*cohort_no[MAMMALS]);
        j = max(j, 0);
        j = min(j, cohort_no[MAMMALS] - 1);
        if (pop[MAMMALS][j][number] > 0) continue;                            // Already resurrected, select a new individual
        pop[MAMMALS][j][number] = 1.0;
        /*			  if (islactating(j))													// If female is lactating, also resurrect her calf
         {
         for (k = (cohort_no[MAMMALS] - 1); k > -1 ; k--)
         {
         if (iszero(popIDcard[MAMMALS][j][IDindividualID] - popIDcard[MAMMALS][k][IDmotherID]))
         {
         pop[MAMMALS][k][number] = 1.0;
         break;
         }
         }				
         }
         */
        if (isfemale(j) && !iscalve(j)) females++;							// Only count non-calf females
        break;
      }
      if (females >= RESCALEPOP) break;
    }
#endif
    
    SievePop();
    
    // Reset the IDmotherIndex of all individuals
    for (i = 0; i < cohort_no[MAMMALS]; i++)
    {
      popIDcard[MAMMALS][i][IDmotherIndex] = -1;
      if (iscalve(i))
      {
        for (j = 0; j < cohort_no[MAMMALS]; j++)
          if (iszero(popIDcard[MAMMALS][j][IDindividualID] - popIDcard[MAMMALS][i][IDmotherID]))
          {
            if (j >= i)
            {
              sprintf(msg, "UserInit(): Invalid value of IDmotherIndex %d: Cannot equal or exceed child's index %d!", j, i);
              ErrorExit(1, msg);
            }
            
            popIDcard[MAMMALS][i][IDmotherID]    = j;
            popIDcard[MAMMALS][i][IDmotherIndex] = j;
            
            popIDcard[MAMMALS][j][IDlactating]   = 1;
            break;
          }
      }
    }
    for (i = 0; i < cohort_no[MAMMALS]; i++) popIDcard[MAMMALS][i][IDindividualID] = i;
  }
  individualIDcount = cohort_no[MAMMALS];
  initialIDcount = individualIDcount;
  
  ReportNote("\n    %s\n    %s%s",
#if (BIFURCATION == 1)
    "Bifurcation of the individual-based marine mammal model",
#else
    "Dynamics of the individual-based marine mammal model",
#endif /* BIFURCATION */
#if (ENVIRON_DIM == 2)
    "Density-dependence through competition for limiting, dynamic resource",
#else
    "Density-independent population growth at constant resource density",
#endif /* BIFURCATION */
    "\n");
  
  ReportNote("%-67s:  %ld", "Seeding value of the random-number generator", UniDevSeed);
  ReportNote("%-67s:  S=%.1f kg (W=%.1f kg)", "Calve size at birth", SizeAtBirth, SizeAtBirth/(1.0 - RhoS));
  ReportNote("%-67s:  NeonateWeightCosts =%.2f kg", "Total costs (in kg fat) of neonate structure and reserves", NeonateWeightCosts);
  ReportNote("%-67s:  NeonateMaintCosts  =%.2f kg", "Cumulative costs (in kg fat) of foetal maintenance requirements", NeonateMaintCosts);
  //  ReportNote("%-67s:  NeonateMaintCosts  =%.2f kg", "Cumulative costs (in kg fat) of lactating during lactation period", NeonateLactCosts);
  ReportNote("%-67s:  TotalNeonateCosts  =%.2f kg", "Total costs of neonate determining pregnancy threshold", TotalNeonateCosts);
  ReportNote("%-67s:  Model version with separate lenght-age relationship and mortality rates for males. Commit: e7688e0 (July 2019)");
  
  UpdateIDcards(env, pop);
  memset(Deaths, 0, SEXCLASSES*AGECLASSES*sizeof(int));
  
  strcpy(fn, argv[1]);                                                                  // Open additional file for
  strcat(fn, "_deaths.out");		                                                    // death statistics
  deathfile = fopen(fn, "w");
  
  strcpy(fn, argv[1]);
  strcat(fn, "_R0.dat");
  R0file = fopen(fn, "a");
  
  return;
}


/*
 *====================================================================================================================================
 *
 *  SPECIFICATION OF THE NUMBER AND VALUES OF BOUNDARY POINTS
 *
 *====================================================================================================================================
 */

void SetBpointNo(double *env, population *pop, int *bpoint_no)
  
{
  bpoint_no[MAMMALS] = 0;
  
  return;
}


void SetBpoints(double *env, population *pop, population *bpoints) { return; }


/*
 *====================================================================================================================================
 *
 *  SPECIFICATION OF DERIVATIVES
 *
 *====================================================================================================================================
 */

void Gradient(double *env, population *pop, population *ofs, double *envgrad, population *popgrad, population *ofsgrad, population *bpoints)
  
{
  int     i;
#if (ENVIRON_DIM == 2)
  double  total_ringest = 0.0;
#endif
  
  UpdateIDcards(env, pop);
  
  for (i = 0; i < cohort_no[MAMMALS]; i++)
  {
    popgrad[MAMMALS][i][number]    = 0.0;
    popgrad[MAMMALS][i][age]       = 1.0;
    popgrad[MAMMALS][i][reserves]  = -1.0;
    popgrad[MAMMALS][i][totalmort] = popIDcard[MAMMALS][i][IDmortality];
    
    if (isdead(i) || isdying(i)) continue;
    
    popgrad[MAMMALS][i][reserves]  = popIDcard[MAMMALS][i][IDnet_energy];
    popgrad[MAMMALS][i][reserves] /= (popIDcard[MAMMALS][i][IDnet_energy] > 0) ? EpsPlus : EpsMin;
    
#if (ENVIRON_DIM == 2)
    total_ringest += popIDcard[MAMMALS][i][IDringest];
#endif
  }
  
  envgrad[0] = 1.0;
  
#if (ENVIRON_DIM == 2)
  envgrad[1]  = (RESOURCETURNOVER*RESOURCE)*(1.0 - RESOURCEAMPLITUDE*cos(2*M_PI*(time/YEAR)));
  envgrad[1] -= RESOURCETURNOVER*resource + total_ringest/VOLUME;
#endif
  
  return;
}


/*
 *====================================================================================================================================
 *
 *  SPECIFICATION OF EVENT LOCATION AND DYNAMIC COHORT CLOSURE
 *
 *====================================================================================================================================
 */

void EventLocation(double *env, population *pop, population *ofs, population *bpoints, double *events)
  
{
  int     i;
  double  eval;
  
  UpdateIDcards(env, pop);
  
  events[0] = -1.0;                                                                 // Event signalling giving birth
  events[1] = -1.0;                                                                 // Event signalling end of lactation
  events[2] = -1.0;                                                                 // Event signalling conception
  events[3] = -1.0;                                                                 // Event signalling entering the waiting period
  
  for (i = 0; i < cohort_no[MAMMALS]; i++)
  {
    if (ismale(i)) continue;
    
    if (ispregnant(i) && islactating(i))  										// A pregnant and lactating mother must first wean
    {
      eval      = pop[MAMMALS][i][age]/(popIDcard[MAMMALS][i][IDcalved] + Tl);
      events[1] = max(events[1], eval - 1.0);
    }
    else if (ispregnant(i))															// A pregnant, but non-lactating mother gives birth
    {
      eval      = pop[MAMMALS][i][age]/(popIDcard[MAMMALS][i][IDinseminated] + Tg);
      events[0] = max(events[0], eval - 1.0);
    }
    else																			// A non-pregnant mother can wean (if lactating) and become pregnant (or enter waiting period)
    {
      if (islactating(i))														// A lactating, non-pregnant mother weans
      {
        eval      = pop[MAMMALS][i][age]/(popIDcard[MAMMALS][i][IDcalved] + Tl);
        events[1] = max(events[1], eval - 1.0);
      }
      
      if (iswaiting(i))															// A waiting mother becomes pregnant when waiting time is due
      {
        eval	  = pop[MAMMALS][i][age]/(popIDcard[MAMMALS][i][IDpregwait]);
        events[2] = max(events[2], eval - 1.0);
      }
      else																		// Individual enters waiting period when pregnancy condition is met
      {
        eval      = min((pop[MAMMALS][i][reserves]/(RhoS*popIDcard[MAMMALS][i][IDweight] + TotalNeonateCosts)), (pop[MAMMALS][i][age] / (Tl + 365)));
        if (islactating(i)) eval = min(eval, pop[MAMMALS][i][age]/(popIDcard[MAMMALS][i][IDcalved] + Tl - Tg + 1.0)); // Enter waiting time 1 day after entering the last year of lactation
        events[3] = max(events[3], eval - 1.0);
      }
    }
  }
  return;
}


/*==================================================================================================================================*/

int ForceCohortEnd(double *env, population *pop, population *ofs, population *bpoints)
  
{
  int     i, j;
  double  eval;
  double  pregwait;
  
  // One or the other mother is coming to terms, so we have to end the cohort cycle to add a new cohort (individual)
  if (LocatedEvent == 0) return COHORT_END;
  
  UpdateIDcards(env, pop);
  
  // For other events (start of pregnancy or end of lactation period) we do not have to end the cohort cycle
  for (i = 0; i < cohort_no[MAMMALS]; i++)
  {
    if (ismale(i)) continue;
    
    // Reset the lactating flag for mothers that have calved more than Tl days ago
    eval = pop[MAMMALS][i][age]/(popIDcard[MAMMALS][i][IDcalved] + Tl);
    if (islactating(i) && isone(eval))
    {
      popIDcard[MAMMALS][i][IDlactating] = 0;
      // Calculate the R0 based on weaned calves
      popIDcard[MAMMALS][i][IDtotalweaned] += 1.0;
      if (popIDcard[MAMMALS][i][IDagefirstweaning] == 0) popIDcard[MAMMALS][i][IDagefirstweaning] = pop[MAMMALS][i][age];
      // Update the IWI statistic when weaning has happened before
      if (popIDcard[MAMMALS][i][IDweaned] > 0)
      {
        ObsIWI++;
        UpdateStats(pop[MAMMALS][i][age] - popIDcard[MAMMALS][i][IDweaned], &meanIWI, &ssIWI, &minIWI, &maxIWI, ObsIWI, 0);
      }
      // Set IDcalved after updating IWI
      popIDcard[MAMMALS][i][IDweaned] = pop[MAMMALS][i][age];
      
      // Find the calve of the mother with index i that stops lactating. Starting at the youngest cohort will lead to a hit sooner
      for (j = (cohort_no[MAMMALS] - 1); j > -1 ; j--)
      {
        if (iszero(popIDcard[MAMMALS][i][IDindividualID] - popIDcard[MAMMALS][j][IDmotherID]))
        {
          // If so, reset the mother's index of the calve and go to the next mother that stops lactating
          popIDcard[MAMMALS][j][IDmotherIndex] = -1;
          break;
        }
      }
    }
    
    // Check which mother's are ready for pregnancy. Do this after resetting lactation flags to allow for immediate restart of pregnancy after lactation
    
    // Waiting mother is ready to become pregnant
    eval = pop[MAMMALS][i][age]/(popIDcard[MAMMALS][i][IDpregwait]);
    if (!ispregnant(i) && iswaiting(i) && isone(eval))
    {
      popIDcard[MAMMALS][i][IDinseminated] = pop[MAMMALS][i][age];
      popIDcard[MAMMALS][i][IDpregnant]    = 1;
      popIDcard[MAMMALS][i][IDwaiting]     = 0;	
    }
    
    // A non-waiting mother is ready to enter the waiting period		
    eval = min((pop[MAMMALS][i][reserves]/(RhoS*popIDcard[MAMMALS][i][IDweight] + TotalNeonateCosts)), (pop[MAMMALS][i][age] / (Tl + 365)));
    if (islactating(i)) eval = min(eval, pop[MAMMALS][i][age]/(popIDcard[MAMMALS][i][IDcalved] + Tl - Tg + 1.0));
    
    if (!ispregnant(i) && !iswaiting(i) && isone(eval))
    {
      // Calculate the waiting time
#if (FIXED_PREGDELAY == 1)
      pregwait = 445;
#else
      double tmp;
      
      tmp = UniDev();
      tmp = max(tmp, DBL_MIN);
      tmp = min(tmp, 1.0 - DBL_MIN);
      pregwait = -log(tmp) / PREGCHANCE;
#endif			
      popIDcard[MAMMALS][i][IDpregwait]		= pop[MAMMALS][i][age] + pregwait;
      popIDcard[MAMMALS][i][IDwaiting]		= 1;
      if (popIDcard[MAMMALS][i][IDagefirstreceptive] == 0) popIDcard[MAMMALS][i][IDagefirstreceptive] = pop[MAMMALS][i][age];
    }
  }
  return NO_COHORT_END;
}


/*
 *====================================================================================================================================
 *
 *  SPECIFICATION OF BETWEEN COHORT CYCLE DYNAMICS
 *
 *====================================================================================================================================
 */

void InstantDynamics(double *env, population *pop, population *ofs)
  
{
  int           i, j, mother;
  //  int			k
  int           ExistingCohortNo, calve, n, sex, death = 0, warn;
  double        tmp, eval, pregwait;
  static double lasttime = -1.0;
#if (RESCALEPERIOD > 0)
  int           females;
  static double lastrescale = 0.0;
#endif
  
  // Check for double execution in a row
  if (time < (lasttime + TINY*cohort_limit)) return;
  lasttime = time;
  
  // Reinitialize statistics and observations during bifucation run every BifOutput Period
#if (BIFURCATION == 1)
  double dummy;
  
  dummy = fmod(time - BifPeriod + BifOutput, BifPeriod);
  
  if (dummy < TINY)
  {
    Observations = ObsAgeFirstReceptive = ObsAgeFirstBirth = ObsAgeFirstWeaning = ObsStarvDays = ObsIBI = ObsIWI = 0.0;
    
    meanR0 = ssR0 = minR0 = maxR0 = 0.0;
    meanTotalCalves = ssTotalCalves = minTotalCalves = maxTotalCalves = 0.0;
    meanTotalWeaned = ssTotalWeaned = minTotalWeaned = maxTotalWeaned = 0.0;
    meanAgeFirstReceptive = ssAgeFirstReceptive = minAgeFirstReceptive = maxAgeFirstReceptive = 0.0;
    meanAgeFirstBirth = ssAgeFirstBirth = minAgeFirstBirth = maxAgeFirstBirth = 0.0;
    meanAgeFirstWeaning = ssAgeFirstWeaning = minAgeFirstWeaning = maxAgeFirstWeaning = 0.0;
    meanStarvDays = ssStarvDays = minStarvDays = maxStarvDays = 0.0;
    meanReproFemales = ssReproFemales = 0.0;
    meanStarvingFemales = ssStarvingFemales = 0.0;
    meanIBI = ssIBI = minIBI = maxIBI = 0.0;
    meanIWI = ssIWI = minIWI = maxIWI = 0.0;
  }
#endif
  
#if (RESCALEPERIOD > 0)
  if (time >= (1.0 + TINY)*(lastrescale + RESCALEPERIOD*YEAR))
  {
    // Kill all individuals
    for (i = 0; i < cohort_no[MAMMALS]; i++) pop[MAMMALS][i][number]  = -1.0;
    
    // Resurrect a specific number of individuals at random (aselective sampling without replacement) until RESCALEPOP females remain
    for (i = 0, females = 0; i < cohort_no[MAMMALS]; i++)
    {
      while (1)
      {
        j = (int)floor(UniDev()*cohort_no[MAMMALS]);
        j = max(j, 0);
        j = min(j, cohort_no[MAMMALS] - 1);
        if (pop[MAMMALS][j][number] > 0) continue;                            // Already resurrected, select a new individual
        pop[MAMMALS][j][number] = 1.0;
        /*
        if (islactating(j))													// If female is lactating, also resurrect her calf
        {
          for (k = (cohort_no[MAMMALS] - 1); k > -1 ; k--)
          {
            if (iszero(popIDcard[MAMMALS][j][IDindividualID] - popIDcard[MAMMALS][k][IDmotherID]))
            {
              pop[MAMMALS][k][number] = 1.0;
              break;
            }
          }				
        }
       */
        if (isfemale(j) && !iscalve(j)) females++;							// Only count non-calf females
        break;
      }
      if (females >= RESCALEPOP) break;
    }
    lastrescale = time;
  }
#endif
  
  UpdateIDcards(env, pop);
  
  //================================================================================
  // Remove dead and totally starved individuals as well as individuals that are currently dying of starvation
  for (i = 0; i < cohort_no[MAMMALS]; i++)
  {
    
    // Death occurring because of old-age, natural death or starvation
    if (isdead(i) || isdying(i)) pop[MAMMALS][i][number] = -1.0;
    
#if (ENVIRON_DIM != 2)
    // Discard males after their first year of life for computational efficiency
    if (ismale(i) && (!iscalve(i))) pop[MAMMALS][i][number] = -1.0;
#endif
    
    // Skip the rest of this loop for survivors
    if (!isdead(i)) continue;
    
    death = 1;
    if (isfemale(i))
    {
      Observations++;
      if (popIDcard[MAMMALS][i][IDagefirstreceptive] > 0) ObsAgeFirstReceptive++;
      if (popIDcard[MAMMALS][i][IDagefirstbirth] > 0) ObsAgeFirstBirth++;
      if (popIDcard[MAMMALS][i][IDagefirstweaning] > 0) ObsAgeFirstWeaning++;
      if (popIDcard[MAMMALS][i][IDstarvdays] > 0) ObsStarvDays++;
      
      UpdateStats(popIDcard[MAMMALS][i][IDtotalfemcalves], &meanR0, &ssR0, &minR0, &maxR0, Observations, 1);
      UpdateStats(popIDcard[MAMMALS][i][IDtotalcalves], &meanTotalCalves, &ssTotalCalves, &minTotalCalves, &maxTotalCalves, Observations, 1);
      UpdateStats(popIDcard[MAMMALS][i][IDtotalweaned], &meanTotalWeaned, &ssTotalWeaned, &minTotalWeaned, &maxTotalWeaned, Observations, 1);
      UpdateStats(popIDcard[MAMMALS][i][IDagefirstreceptive], &meanAgeFirstReceptive, &ssAgeFirstReceptive, &minAgeFirstReceptive, &maxAgeFirstReceptive, ObsAgeFirstReceptive, 0);
      UpdateStats(popIDcard[MAMMALS][i][IDagefirstbirth], &meanAgeFirstBirth, &ssAgeFirstBirth, &minAgeFirstBirth, &maxAgeFirstBirth, ObsAgeFirstBirth, 0);
      UpdateStats(popIDcard[MAMMALS][i][IDagefirstweaning], &meanAgeFirstWeaning, &ssAgeFirstWeaning, &minAgeFirstWeaning, &maxAgeFirstWeaning, ObsAgeFirstWeaning, 0);
      UpdateStats(popIDcard[MAMMALS][i][IDstarvdays], &meanStarvDays, &ssStarvDays, &minStarvDays, &maxStarvDays, ObsStarvDays, 0);
      UpdateStats((double)(popIDcard[MAMMALS][i][IDtotalcalves] > 0), &meanReproFemales, &ssReproFemales, &minReproFemales, &maxReproFemales, Observations, 1);
      UpdateStats((double)(popIDcard[MAMMALS][i][IDstarvdays] > 0), &meanStarvingFemales, &ssStarvingFemales, &minStarvingFemales, &maxStarvingFemales, Observations, 1);
      
      fprintf(R0file, "%16.8G\t %16.8G\t %16.8G\t %16.8G\t %16.8G\t %16.8G\t %16.8G\t %16.8G\t %16.8G\t %16.8G\t %16.8G\t %16.8G\t %16.8G\t %16.8G\t %16.8G\t %16.8G\n", 
        DISTURBANCEPERIOD, RESOURCEAMPLITUDE, time, resource,
        popIDcard[MAMMALS][i][IDindividualID], pop[MAMMALS][i][age], pop[MAMMALS][i][reserves], exp(-pop[MAMMALS][i][totalmort]), exp(-popIDcard[MAMMALS][i][IDminloglifespan]),
        popIDcard[MAMMALS][i][IDtotalfemcalves], popIDcard[MAMMALS][i][IDtotalcalves], popIDcard[MAMMALS][i][IDtotalweaned],
        popIDcard[MAMMALS][i][IDagefirstreceptive], popIDcard[MAMMALS][i][IDagefirstbirth], popIDcard[MAMMALS][i][IDagefirstweaning],
        popIDcard[MAMMALS][i][IDstarvdays]);
      
      // Find the calve of the lactating mother with index i that is dying. Starting at the youngest cohort will lead to a hit sooner
      if (islactating(i))
      {
        for (j = (cohort_no[MAMMALS] - 1); j > -1 ; j--)
        {
          // Check whether the individual with index j is the calve of the dying individual with index i
          if (iszero(popIDcard[MAMMALS][i][IDindividualID] - popIDcard[MAMMALS][j][IDmotherID]))
          {
            // If so, reset the mother's index of the calve and go to the next dying individual
            popIDcard[MAMMALS][j][IDmotherIndex] = -1;
            break;
          }
        }
      }
    }
    
    if (iscalve(i) && ((popIDcard[MAMMALS][i][IDmotherIndex] + TINY) >= 0))       // Consider dead calves that still have a mother
    {
      // Find the mother of this calve and reset her lactation flag
      mother = (int)(popIDcard[MAMMALS][i][IDmotherIndex] + TINY);
      
      // Do some error checking: mother's can not have an index larger than their child as they are older
      if (mother >= i)
      {
        sprintf(msg, "InstantDynamics(): Invalid value of IDmotherIndex %d: Cannot equal or exceed child's index %d!", mother, i);
        ErrorExit(1, msg);
      }
      else if (!iszero(popIDcard[MAMMALS][i][IDmotherID] - popIDcard[MAMMALS][mother][IDindividualID]))
      {
        sprintf(msg, "InstantDynamics(): IDmotherID: %.0f IDmotherIndex: %.0f mother: %d Mother's IDindividualID: %.0f", popIDcard[MAMMALS][i][IDmotherID],
          popIDcard[MAMMALS][i][IDmotherIndex], mother, popIDcard[MAMMALS][mother][IDindividualID]);
        ErrorExit(1, msg);
      }
      else
      {
        popIDcard[MAMMALS][mother][IDlactating]   = 0;
        
        // Check if mother can enter waiting period again. Leave mothers already waiting alone.
        eval = min((pop[MAMMALS][mother][reserves]/(RhoS*popIDcard[MAMMALS][mother][IDweight] + TotalNeonateCosts)), (pop[MAMMALS][mother][age] / (Tl + 365)));
        if ((!ispregnant(mother)) && (!iswaiting(mother)) && isone(eval))
        {
#if (FIXED_PREGDELAY == 1)
          pregwait = 445;
#else
          tmp = UniDev();
          tmp = max(tmp, DBL_MIN);
          tmp = min(tmp, 1.0 - DBL_MIN);
          pregwait = -log(tmp) / PREGCHANCE;
#endif
          popIDcard[MAMMALS][mother][IDpregwait]  = pop[MAMMALS][mother][age] + pregwait;
          popIDcard[MAMMALS][mother][IDwaiting]   = 1;
        }
      }
    }
  }
  
  //================================================================================
  // Check the pregnancy and lactating flags for mothers and reset the mother index for calves
  ExistingCohortNo = cohort_no[MAMMALS];
  for (i = 0; i < ExistingCohortNo; i++)
  {
    // Skip this loop for males, calves and dead individuals as it is about pregnancy and lactation
    if (ismale(i) || isdead(i) || iscalve(i)) continue;
    
    warn = 1;
    // Reset the lactating flag for mothers that have calved more than Tl days ago
    eval = pop[MAMMALS][i][age]/(popIDcard[MAMMALS][i][IDcalved] + Tl);
    if (islactating(i) && isone(eval))
    {
      popIDcard[MAMMALS][i][IDlactating] = 0;
      popIDcard[MAMMALS][i][IDtotalweaned] += 1.0;
      if (popIDcard[MAMMALS][i][IDagefirstweaning] == 0) popIDcard[MAMMALS][i][IDagefirstweaning] = pop[MAMMALS][i][age];
      // Update the IWI statistic when weaning has happened before
      if (popIDcard[MAMMALS][i][IDweaned] > 0)
      {
        ObsIWI++;
        UpdateStats(pop[MAMMALS][i][age] - popIDcard[MAMMALS][i][IDweaned], &meanIWI, &ssIWI, &minIWI, &maxIWI, ObsIWI, 0);
      }
      // Set IDcalved after updating IWI
      popIDcard[MAMMALS][i][IDweaned] = pop[MAMMALS][i][age];
      
#if (!DEBUG)
      if (fabs(pop[MAMMALS][i][age]/(popIDcard[MAMMALS][i][IDcalved] + Tl) - 1.0) > TINY)
#endif
{
  fprintf(stderr,
    "T = %12.3f: InstantDynamics(): Missed end of lactation period.                               Individual ID %6.0f; Age     : %12.3f; Should have stopped at: %12.3f  (%6.4f)\n",
    time, popIDcard[MAMMALS][i][IDindividualID], pop[MAMMALS][i][age], (popIDcard[MAMMALS][i][IDcalved] + Tl), 
    (pop[MAMMALS][i][age]/(popIDcard[MAMMALS][i][IDcalved] + Tl)));
}
      // no need to find the calve of the mother with index i that stops lactating. Is taking care of in the next loop
      warn = 0;
    }
    
    eval = min((pop[MAMMALS][i][reserves]/(RhoS*popIDcard[MAMMALS][i][IDweight] + TotalNeonateCosts)), (pop[MAMMALS][i][age] / (Tl + 365)));
    if (islactating(i)) eval = min(eval, pop[MAMMALS][i][age]/(popIDcard[MAMMALS][i][IDcalved] + Tl - Tg + 1.0));
    
    if ((!iswaiting(i)) && (!ispregnant(i)) && isone(eval))	// Lactating mother becomes receptive again and enters waiting period
    {
#if (FIXED_PREGDELAY == 1)
      pregwait = 445;
#else
      tmp = UniDev();
      tmp = max(tmp, DBL_MIN);
      tmp = min(tmp, 1.0 - DBL_MIN);
      pregwait = -log(tmp) / PREGCHANCE;
#endif
      popIDcard[MAMMALS][i][IDpregwait]   = pop[MAMMALS][i][age] + pregwait;
      popIDcard[MAMMALS][i][IDwaiting]    = 1;
      if (popIDcard[MAMMALS][i][IDagefirstreceptive] == 0) popIDcard[MAMMALS][i][IDagefirstreceptive] = pop[MAMMALS][i][age];
#if (!DEBUG)
      eval = fabs(min((pop[MAMMALS][i][reserves]/(RhoS*popIDcard[MAMMALS][i][IDweight] + TotalNeonateCosts) - 1.0), (pop[MAMMALS][i][age] / (Tl + 365) - 1.0)));
      if (islactating(i)) eval = min(eval, fabs(pop[MAMMALS][i][age]/(popIDcard[MAMMALS][i][IDcalved] + Tl - Tg + 1.0) - 1.0));
      if (eval > TINY)
#endif
        if (warn)
        {
          if (islactating(i))
            fprintf(stderr,
              "T = %12.3f: InstantDynamics(): Missed start of waiting period of lactating individual.       Individual ID %6.0f; Reserves: %12.6f; Should have started at: %12.6f  (%6.4f)\n",
              time, popIDcard[MAMMALS][i][IDindividualID], pop[MAMMALS][i][reserves], (RhoS *popIDcard[MAMMALS][i][IDweight] + TotalNeonateCosts),
              (pop[MAMMALS][i][reserves]/(RhoS *popIDcard[MAMMALS][i][IDweight] + TotalNeonateCosts)));
          else
            fprintf(stderr,
              "T = %12.3f: InstantDynamics(): Missed start of waiting period of non-lactating individual.   Individual ID %6.0f; Reserves: %12.6f; Should have started at: %12.6f  (%6.4f)\n",
              time, popIDcard[MAMMALS][i][IDindividualID], pop[MAMMALS][i][reserves], (RhoS *popIDcard[MAMMALS][i][IDweight] + TotalNeonateCosts),
              (pop[MAMMALS][i][reserves]/(RhoS *popIDcard[MAMMALS][i][IDweight] + TotalNeonateCosts)));
        }
    }
    else if ((iswaiting(i)) && (!ispregnant(i)))
    {
      eval = pop[MAMMALS][i][age]/(popIDcard[MAMMALS][i][IDpregwait]);  // Mother waiting time is over
      if (isone(eval))
      {
        popIDcard[MAMMALS][i][IDinseminated] = pop[MAMMALS][i][age];
        popIDcard[MAMMALS][i][IDpregnant]    = 1;
        popIDcard[MAMMALS][i][IDwaiting]     = 0;
#if (!DEBUG)
        eval = fabs(pop[MAMMALS][i][age]/(popIDcard[MAMMALS][i][IDpregwait]) - 1.0);
        if (eval > TINY)
#endif
          if (warn)
          {
            if (islactating(i))
              fprintf(stderr,
                "T = %12.3f: InstantDynamics(): Missed start of pregnancy period of lactating individual.     Individual ID %6.0f; Age     : %12.3f; Should have stopped at: %12.3f  (%6.4f)\n",
                time, popIDcard[MAMMALS][i][IDindividualID], pop[MAMMALS][i][age], popIDcard[MAMMALS][i][IDpregwait],
                (pop[MAMMALS][i][age]/popIDcard[MAMMALS][i][IDpregwait]));
            else
              fprintf(stderr,
                "T = %12.3f: InstantDynamics(): Missed start of pregnancy period of non-lactating individual. Individual ID %6.0f; Age     : %12.3f; Should have stopped at: %12.3f  (%6.4f)\n",
                time, popIDcard[MAMMALS][i][IDindividualID], pop[MAMMALS][i][age], popIDcard[MAMMALS][i][IDpregwait],
                (pop[MAMMALS][i][age]/popIDcard[MAMMALS][i][IDpregwait]));
          }
      }
    }
    else if (ispregnant(i))
    {
      eval = pop[MAMMALS][i][age]/(popIDcard[MAMMALS][i][IDinseminated] + Tg);  // Mother is coming to term
      if (isone(eval))
      {
        // Give birth to a calve
        calve                           = AddCohorts(pop, MAMMALS, 1);
        pop[MAMMALS][calve][number]     = 1.0;
        pop[MAMMALS][calve][age]        = 0.0;
        pop[MAMMALS][calve][reserves]   = RhoS*SizeAtBirth/(1-RhoS);
        pop[MAMMALS][calve][totalmort]  = 0.0;
        
        for (n = 0; n < I_CONST_DIM; n++) popIDcard[MAMMALS][calve][n] = 0.0;
        
        sex = (UniDev() < 0.5);
        popIDcard[MAMMALS][calve][IDmale]         = sex;
        popIDcard[MAMMALS][calve][IDfemale]       = !sex;
        popIDcard[MAMMALS][calve][IDpregwait]     = -1.0;
        popIDcard[MAMMALS][calve][IDinseminated]  = -1.0;
        popIDcard[MAMMALS][calve][IDcalved]       = -1.0;
        popIDcard[MAMMALS][calve][IDindividualID] = (double)individualIDcount++;
        popIDcard[MAMMALS][calve][IDmotherID]     = popIDcard[MAMMALS][i][IDindividualID];
        popIDcard[MAMMALS][calve][IDmotherIndex]  = i;
        
        // Generate a random number between 0 and 1. Store this as minus the logarithm of this value
        // The individual will die of natural mortality if integrated mortality rate exceeds this random deviate.
        tmp = UniDev();
        tmp = max(tmp, DBL_MIN);
        tmp = min(tmp, 1.0 - DBL_MIN);
        popIDcard[MAMMALS][calve][IDminloglifespan] = -log(tmp);
        
        // Substract calve fat from mother's reserves
        pop[MAMMALS][i][reserves] -= pop[MAMMALS][calve][reserves];
        
        // Reset the pregnant flag and set the lactation flag
        popIDcard[MAMMALS][i][IDwaiting]      = 0.0;
        popIDcard[MAMMALS][i][IDpregnant]     = 0.0;
        popIDcard[MAMMALS][i][IDlactating]    = 1.0;
        if (popIDcard[MAMMALS][i][IDagefirstbirth] == 0) popIDcard[MAMMALS][i][IDagefirstbirth] = pop[MAMMALS][i][age];
        // For every birth, update the IBI when its not the first birth
        if (popIDcard[MAMMALS][i][IDcalved] > 0)
        {
          ObsIBI++;
          UpdateStats(pop[MAMMALS][i][age] - popIDcard[MAMMALS][i][IDcalved], &meanIBI, &ssIBI, &minIBI, &maxIBI, ObsIBI, 0);
        }
        // IDcalved should be set after updating IBI!
        popIDcard[MAMMALS][i][IDcalved]       = pop[MAMMALS][i][age];
        popIDcard[MAMMALS][i][IDtotalcalves] += 1.0;                         	   // Count all offspring
        popIDcard[MAMMALS][i][IDtotalfemcalves] += !sex;                         // Count only female offspring
      }
    }
  }
  
  //================================================================================
  YearlyDeaths(env, pop, deathfile);                                                // Classify the number of deaths
  
  // Some individuals have died, this will mean that the cohort index of mothers has changed.
  if (death)
  {
    SievePop();
    
    // Now reset the IDmotherIndex of all individuals
    for (i = 0; i < cohort_no[MAMMALS]; i++)
    {
      popIDcard[MAMMALS][i][IDmotherIndex] = -1;
      // Only for calves that are still suckling set the mother's index to something else than -1
      if (iscalve(i))
      {
        // Check whether the individual with index j is the mother of the suckling calve with index i
        for (j = 0; j < cohort_no[MAMMALS]; j++)
          if (islactating(j) && iszero(popIDcard[MAMMALS][j][IDindividualID] - popIDcard[MAMMALS][i][IDmotherID]))
          {
            // If so set the mother index of the suckling calve with index i equal to
            // the mother's index j and skip to the next individual
            popIDcard[MAMMALS][i][IDmotherIndex] = j;
            break;
          }
      }
    }
  }
  
  return;
}


/*
 *====================================================================================================================================
 *
 *  SPECIFICATION OF OUTPUT VARIABLES
 *
 *====================================================================================================================================
 */

void DefineOutput(double *env, population *pop, double *output)
  
  /*
   *  DefineOutput  - Routine defines the output variables and writes numbers and fatratio's 
   *                  to additional output files
   *
   *      Array output[]:
   *
   *      Column Index    Variable
   *
   *           1     -    Time
   *           2     0    Time in years
   *
   *           3     1    Total number of all individuals
   *           4     2    Total number of males
   *           5     3    Total number of all females
   *           6     4    Total number of resting females
   *           7     5    Total number of pregnant females
   *           8     6    Total number of lactating females
   *           9     7    Total number of pregnant & lactating females
   *          10     8    Total number of waiting females
   *          11     9    Total number of waiting & lacatating females
   *
   *          12    10    Total number of all individuals    0 -  1 year old
   *          13    11    Total number of males              0 -  1 year old
   *          14    12    Total number of all females        0 -  1 year old
   *
   *          15    13    Total number of all individuals    1 -  Tl year old
   *          16    14    Total number of males              1 -  Tl year old
   *          17    15    Total number of all females        1 -  Tl year old
   *
   *          18    16    Total number of all individuals    Tl - 8 year old
   *          19    17    Total number of males              Tl - 8 year old
   *          20    18    Total number of all females        Tl - 8 year old
   *
   *          21    19    Total number of all individuals    8 - 15 year old
   *          22    20    Total number of males              8 - 15 year old
   *          23    21    Total number of all females        8 - 15 year old
   *
   *          24    22    Total number of all individuals   15 - 25 year old
   *          25    23    Total number of males             15 - 25 year old
   *          26    24    Total number of all females       15 - 25 year old
   *
   *          27    25    Total number of all individuals      > 25 year old
   *          28    26    Total number of males                > 25 year old
   *          29    27    Total number of all females          > 25 year old
   *
   *          30    28    Mean R0
   *          31    29    Std. in R0
   *          32    30    Minimum in R0
   *          33    31    Maximum in R0
   *          34    32    Resource density
   *
   */
{
  int   i;
  
  UpdateIDcards(env, pop);
  
  output[0] = time/YEAR;
  for (i = 0; i < cohort_no[MAMMALS]; i++)
  {
    output[ 1] += 1.0;
    output[ 2] += popIDcard[MAMMALS][i][IDmale];
    output[ 3] += popIDcard[MAMMALS][i][IDfemale];
    if (isfemale(i) && !iswaiting(i) && !ispregnant(i) && !islactating(i))  output[ 4] += 1.0; // Resting females 
    if (isfemale(i) && !iswaiting(i) &&  ispregnant(i) && !islactating(i))  output[ 5] += 1.0; // Pregnant females 
    if (isfemale(i) && !iswaiting(i) && !ispregnant(i) &&  islactating(i))  output[ 6] += 1.0; // Lactating females 
    if (isfemale(i) && !iswaiting(i) &&  ispregnant(i) &&  islactating(i))  output[ 7] += 1.0; // Pregnant & Lactating females 
    if (isfemale(i) &&  iswaiting(i) && !ispregnant(i) && !islactating(i))  output[ 8] += 1.0; // Waiting females
    if (isfemale(i) &&  iswaiting(i) && !ispregnant(i) &&  islactating(i))  output[ 9] += 1.0; // Waiting & Lactating females 
    
    if (pop[MAMMALS][i][age] < (YEAR - TINY))                                     // 0-1 year
    {
      output[10] += 1.0;
      output[11] += popIDcard[MAMMALS][i][IDmale];
      output[12] += popIDcard[MAMMALS][i][IDfemale];
    }
    else if (pop[MAMMALS][i][age] < (Tl - TINY))                            		// 1-Tl year
    {
      output[13] += 1.0;
      output[14] += popIDcard[MAMMALS][i][IDmale];
      output[15] += popIDcard[MAMMALS][i][IDfemale];
    }
    else if (pop[MAMMALS][i][age] < (8.0*YEAR - TINY))                          	// Tl-8 year
    {
      output[16]  += 1.0;
      output[17]  += popIDcard[MAMMALS][i][IDmale];
      output[18]  += popIDcard[MAMMALS][i][IDfemale];
    }
    else if (pop[MAMMALS][i][age] < (15.0*YEAR - TINY))                           // 8-15 year
    {
      output[19]  += 1.0;
      output[20]  += popIDcard[MAMMALS][i][IDmale];
      output[21]  += popIDcard[MAMMALS][i][IDfemale];
    }
    else if (pop[MAMMALS][i][age] < (25*YEAR - TINY))                           	// 15-25 year
    {
      output[22]  += 1.0;
      output[23]  += popIDcard[MAMMALS][i][IDmale];
      output[24]  += popIDcard[MAMMALS][i][IDfemale];
    }
    else                                                                          // > 25 year
    {
      output[25]  += 1.0;
      output[26]  += popIDcard[MAMMALS][i][IDmale];
      output[27]  += popIDcard[MAMMALS][i][IDfemale];
    }
  } 
  
  output[28] = meanR0;
  output[29] = (Observations > 1) ? StdDev(ssR0, Observations) : 0;
  output[30] = minR0;
  output[31] = maxR0;
  
  output[32] = meanTotalCalves;
  output[33] = (Observations > 1) ? StdDev(ssTotalCalves, Observations) : 0;
  output[34] = minTotalCalves;
  output[35] = maxTotalCalves;
  
  output[36] = meanTotalWeaned;
  output[37] = (Observations > 1) ? StdDev(ssTotalWeaned, Observations) : 0;
  output[38] = minTotalWeaned;
  output[39] = maxTotalWeaned;
  
  output[40] = meanAgeFirstReceptive;
  output[41] = (ObsAgeFirstReceptive > 1) ? StdDev(ssAgeFirstReceptive, ObsAgeFirstReceptive) : 0;
  output[42] = minAgeFirstReceptive;
  output[43] = maxAgeFirstReceptive;
  
  output[44] = meanAgeFirstBirth;
  output[45] = (ObsAgeFirstBirth > 1) ? StdDev(ssAgeFirstBirth, ObsAgeFirstBirth) : 0;
  output[46] = minAgeFirstBirth;
  output[47] = maxAgeFirstBirth;
  
  output[48] = meanAgeFirstWeaning;
  output[49] = (ObsAgeFirstWeaning > 1) ? StdDev(ssAgeFirstWeaning, ObsAgeFirstWeaning) : 0;
  output[50] = minAgeFirstWeaning;
  output[51] = maxAgeFirstWeaning;
  
  output[52] = meanIBI;
  output[53] = (ObsIBI > 1) ? StdDev(ssIBI, ObsIBI) : 0;
  output[54] = minIBI;
  output[55] = maxIBI;
  
  output[56] = meanIWI;
  output[57] = (ObsIWI > 1) ? StdDev(ssIWI, ObsIWI) : 0;
  output[58] = minIWI;
  output[59] = maxIWI;
  
  output[60] = meanStarvDays;
  output[61] = (ObsStarvDays > 1) ? StdDev(ssStarvDays, ObsStarvDays) : 0;
  output[62] = minStarvDays;
  output[63] = maxStarvDays;
  
  output[64] = meanReproFemales;
  output[65] = (Observations > 1) ? StdDev(ssReproFemales, Observations) : 0;
  
  output[66] = meanStarvingFemales;
  output[67] = (Observations > 1) ? StdDev(ssStarvingFemales, Observations) : 0;
  
  output[68] = (double)Observations;
  output[69] = (double)ObsAgeFirstReceptive;
  output[70] = (double)ObsAgeFirstBirth;
  output[71] = (double)ObsAgeFirstWeaning;
  output[72] = (double)ObsStarvDays;
  output[73] = (double)ObsIBI;
  output[74] = (double)ObsIWI;
  
#if (ENVIRON_DIM == 2)
  output[75] = resource;
#endif
  
  return;
}


/*
 *====================================================================================================================================
 *
 *  ADDITIONAL, PROBLEM-SPECIFIC FUNCTIONS
 *
 *====================================================================================================================================
 */

static void UpdateIDcards(double *env, population *pop)
  
  /*
   * UpdateIDcards - Updates the values of all IDconst variables for the current time and environmental conditions
   */
  
{
  int           i, mother, notdisturbed;
  double        len, bones, fat, weight, weightM, fatratio;
  double        tau, foetuslen, foetusbones;
  double        ringest, mingest, milkage, tmp, feeding;
  double        maint, growth, pregcosts;
  double        background, starvation;
  double        timeinyear;
  static double lasttimez = 0.0;
  
  timeinyear   = time - floor(time/YEAR)*YEAR;
  notdisturbed = ((timeinyear >= DISTURBANCESTART) && (timeinyear <= (DISTURBANCESTART + DISTURBANCEPERIOD))) ? 0 : 1;
  
  for (i = 0; i < cohort_no[MAMMALS]; i++)
  {
    if (isdead(i)) continue;
    
    if (isfemale(i)) {
      len       = Linf - (Linf - Lb)*bexp(-K*pop[MAMMALS][i][age]);
    } else {
      len       = Linf_m - (Linf_m - Lb)*bexp(-K_m*pop[MAMMALS][i][age]);
    }
    bones     = OmegaS*pow(len, OmegaE);
    fat       = pop[MAMMALS][i][reserves];
    weight    = bones + fat;
    weightM   = bones + OmegaM*fat;
    
    // Add foetus' weight to mother's weight
    pregcosts = 0.0;
    if (isfemale(i) && ispregnant(i))
    {
      tau           = pop[MAMMALS][i][age] - popIDcard[MAMMALS][i][IDinseminated];
      foetuslen     = Lb*tau/Tg;
      foetusbones   = OmegaS*pow(foetuslen, OmegaE);
      weight       += foetusbones;
      weightM      += OmegaC*foetusbones;
      pregcosts     = SigmaG*OmegaS*OmegaE*pow(Lb/Tg, OmegaE)*pow(tau, OmegaE-1);
    }
    
    fatratio  = fat/weight;
    feeding   = 1.0/(1 + bexp(-EtaF*(Rho/fatratio - 1.0)));
    
    // Determine the mortality rates
    if (isfemale(i) || iscalve(i)) {
      background = Alpha1 * bexp(-Beta1*pop[MAMMALS][i][age]) + Alpha2 * bexp(Beta2*pop[MAMMALS][i][age]);
    } else {
      background = Mu_male;
    }
    
    starvation = (fatratio < MINFATRATIO) ? max(MuS*(RhoS/MINFATRATIO - 1.0), 0) : max(MuS*(RhoS/fatratio - 1.0), 0);
    
    // Initialize the IDconst to be reset with default empty values
    popIDcard[MAMMALS][i][IDlength]     = len;
    popIDcard[MAMMALS][i][IDbones]      = bones;
    popIDcard[MAMMALS][i][IDweight]     = weight;
    popIDcard[MAMMALS][i][IDweightM]    = weightM;
    popIDcard[MAMMALS][i][IDfatratio]   = fatratio;
    
    popIDcard[MAMMALS][i][IDingest]     = 0.0;
    popIDcard[MAMMALS][i][IDringest]    = 0.0;
    popIDcard[MAMMALS][i][IDmingest]    = 0.0;
    
    popIDcard[MAMMALS][i][IDmaint]      = 0.0;
    popIDcard[MAMMALS][i][IDgrowth]     = 0.0;
    popIDcard[MAMMALS][i][IDpregcosts]  = 0.0;
    popIDcard[MAMMALS][i][IDlactcosts]  = 0.0;
    popIDcard[MAMMALS][i][IDnet_energy] = 0.0;
    
    popIDcard[MAMMALS][i][IDmortality]    = background + starvation;
    popIDcard[MAMMALS][i][IDbackground]   = background;
    popIDcard[MAMMALS][i][IDstarvation]   = starvation;
    popIDcard[MAMMALS][i][IDfeedinglevel] = feeding;
    
    popIDcard[MAMMALS][i][IDstarvdays]   += (starvation > 0.0) ? max((time - lasttimez),0.0) : 0.0;
    
    // Skip energteics of cohorts that are dying through fat-depletion or that have reached their maximum lifespan
    if (isdying(i)) continue;
    
    // Calculate resource ingestion
    ringest   = notdisturbed*feeding*resource*pow(bones, 2.0/3.0);
    ringest  *= pow(pop[MAMMALS][i][age], Gamma) / (pow(Tr, Gamma) + pow(pop[MAMMALS][i][age], Gamma));
    
    // Calculate milk ingestion if mother is not dead
    mingest   = 0.0;
    if (iscalve(i) && ((popIDcard[MAMMALS][i][IDmotherIndex] + TINY) >= 0))
    {
      // Find the mother's index
      mother = (int)(popIDcard[MAMMALS][i][IDmotherIndex] + TINY);
      
      // Do some error checking: mother's can not have an index larger than their child as they are older
      if (mother >= i)
      {
        sprintf(msg, "UpdateIDcards(): Invalid value of IDmotherIndex %d: Cannot equal or exceed child's index %d!", mother, i);
        ErrorExit(1, msg);
      }
      if (!iszero(popIDcard[MAMMALS][i][IDmotherID] - popIDcard[MAMMALS][mother][IDindividualID]))
      {
        sprintf(msg, "UpdateIDcards(): IDmotherID: %.0f IDmotherIndex: %.0f mother: %d Mother's IDindividualID: %.0f", popIDcard[MAMMALS][i][IDmotherID],
          popIDcard[MAMMALS][i][IDmotherIndex], mother, popIDcard[MAMMALS][mother][IDindividualID]);
        ErrorExit(1, msg);
      }
      
      mingest  = PhiM*pow(bones, 2.0/3.0);
      
      // Adjust lactation to mother's condition
      tmp  = max((1 - XiM)*(pop[MAMMALS][mother][reserves] - RhoS*popIDcard[MAMMALS][mother][IDweight]), 0.0);
      tmp /= (Rho - RhoS)*popIDcard[MAMMALS][mother][IDweight] - XiM*(pop[MAMMALS][mother][reserves] - RhoS*popIDcard[MAMMALS][mother][IDweight]);
      mingest *= max(tmp, 0.0);
      
      // Adjust milk intake to calve's condition
      mingest *= feeding;
      milkage  = max(0, (1 - (pop[MAMMALS][i][age] - Tf)/(Tl-Tf))/(1 - XiC*(pop[MAMMALS][i][age]-Tf)/(Tl-Tf)));
      mingest *= min(1, milkage);
      
      // Add milk costs to the mother's IDcard, both to lactation costs and substract it from net-energy production
      popIDcard[MAMMALS][mother][IDlactcosts]   = max(mingest / SigmaL, 0.0);
      popIDcard[MAMMALS][mother][IDnet_energy] -= popIDcard[MAMMALS][mother][IDlactcosts];
    }
    
    // Compute the maintenance costs
    maint = SigmaM*pow(weightM, 0.75);
    
    // Compute the growth costs
    if (isfemale(i)) {
      growth = OmegaE*SigmaG*K*(Linf - len)*OmegaS*pow(len,OmegaE-1);
    } else {
      growth = OmegaE*SigmaG*K_m*(Linf_m - len)*OmegaS*pow(len,OmegaE-1);
    }
    
    // Initialize the IDconst to be reset with default empty values
    popIDcard[MAMMALS][i][IDingest]     = ringest + mingest;
    popIDcard[MAMMALS][i][IDringest]    = ringest;
    popIDcard[MAMMALS][i][IDmingest]    = mingest;
    
    popIDcard[MAMMALS][i][IDmaint]      = maint;
    popIDcard[MAMMALS][i][IDgrowth]     = growth;
    popIDcard[MAMMALS][i][IDpregcosts]  = pregcosts;
    popIDcard[MAMMALS][i][IDlactcosts]  = 0.0;
    popIDcard[MAMMALS][i][IDnet_energy] = ringest + mingest - maint - growth - pregcosts;
  }
  
  lasttimez = time;
  
  return;
}


/*
 *====================================================================================================================================
 *
 *  USER-DEFINED SUPPORTING ROUTINES
 *
 *====================================================================================================================================
 */

/* SEXCLASSES
 0: All males
 1: All females
 2: Non-waiting, non-pregnant, non-lactating females (resting)
 3: Non-waiting, pregnant, non-lactating females (pregnant)
 4: Non-waiting, non-pregnant, lactating females (lactating)
 5: Non-waiting, pregnant, lactating females (preglact)
 6: Waiting, non-pregnant, non-lactating (waiting)
 7: Waiting, non-pregnant, lactating females (waitlact)
 
 AGECLASSES
 0:  0 -  1 years old (yearling)
 1:  1 - Tl years old (suckling)
 2: Tl -  8 years old (juvenile)
 3:  8 - 15 years old (growing)
 4: 15 - 25 years old (mature)
 5:    > 25 years old (scenescent)
 */

static void   YearlyDeaths(double *env, population *pop, FILE *fpnt)
  
{
  register int  i, j, sex, ageclass;
  double        dummy;
  
  for (i=0; i<cohort_no[MAMMALS]; i++)
  {
    if (isdead(i))                                                                // Classify all individuas currently dying
    {
      if      (ismale(i))                                        sex = 0;
      else if (isfemale(i) && !iswaiting(i) && !ispregnant(i) && !islactating(i)) sex = 2;
      else if (isfemale(i) && !iswaiting(i) &&  ispregnant(i) && !islactating(i)) sex = 3;
      else if (isfemale(i) && !iswaiting(i) && !ispregnant(i) &&  islactating(i)) sex = 4;
      else if (isfemale(i) && !iswaiting(i) &&  ispregnant(i) &&  islactating(i)) sex = 5;
      else if (isfemale(i) &&  iswaiting(i) && !ispregnant(i) && !islactating(i)) sex = 6;
      else sex = 7;
      
      if      (pop[MAMMALS][i][age] < (YEAR + TINY))    ageclass = 0;
      else if (pop[MAMMALS][i][age] < (Tl + TINY))      ageclass = 1;
      else if (pop[MAMMALS][i][age] < (8.0*YEAR+TINY))  ageclass = 2;
      else if (pop[MAMMALS][i][age] < (15.0*YEAR+TINY)) ageclass = 3;
      else if (pop[MAMMALS][i][age] < (25.0*YEAR+TINY)) ageclass = 4;
      else ageclass = 5;
      
      Deaths[sex][ageclass]++;
    }
  }
  
  // Is it new year's eve?
  dummy = fmod(time, YEAR);
  if (fpnt && ((dummy < TINY) || (dummy > YEAR - TINY)))                            // Write out and reset
  {
    for (i=2; i<SEXCLASSES; i++)                                                  // Add up all females
      for (j=0; j<AGECLASSES; j++)
        Deaths[1][j] += Deaths[i][j];
    
    (void)fprintf(fpnt, "%d", (int)(time/YEAR));
    
    for (i=0; i<SEXCLASSES; i++)
      for (j=0; j<AGECLASSES; j++)
      {
        (void)fprintf(fpnt, "\t%5d", Deaths[i][j]);
        Deaths[i][j] = 0;
      }
      (void)fprintf(fpnt, "\n");
    (void)fflush(fpnt);
  }
  
  return;
}


static void UpdateStats(double value, double *mean, double *sum_sq, double *minval, double *maxval, long n, int countzero)
  
  /*
   *  UpdateStats - Routine that updates the value of the mean, the sum of the
   *                squared deviations from the mean, the maximum and minimum
   *                for a series of number of which "value" is the next element.
   *
   *  Arguments   - value   : The current element in the series of numbers.
   *                mean    : The mean value of all preceding number in the series.
   *                sum_sq  : The sum of the squared deviations of all preceding numbers.
   *                minval  : The minimum value of all preceding number in the series.
   *                maxval  : The maximum value of all preceding number in the series.
   *                n       : The order number of the current value in the series.
   *                          It hence equals the number of observations on which
   *                          the resulting (!!) statistics are based.
   */
  
{
  double  deviation;
  
  // Very clever!
  // Writing out the expressions for mean and variance shows that it is indeed
  // correct.
  
  if(value == 0 && countzero == 0) return;
  
  deviation  = value -*mean;
  *mean     += (deviation/n);
  *sum_sq   += (deviation*(value -*mean));
  
  if ((n == 1) || (value >*maxval)) *maxval = value;
  if ((n == 1) || (value <*minval)) *minval = value;
  
  return;
}


/*==================================================================================================================================*/

#define MAX_EXP 50.0

static double bexp(double x)
  
{
  double  pw = x;
  
  pw = max(pw, -MAX_EXP);
  pw = min(pw, MAX_EXP);
  
  return exp(pw);
}


// From here onwards we have to undefine time as an alias for env[0]

#undef time
#include <sys/time.h>
#include <time.h>

#define MAX_STORED                100
#ifndef RAND_MAX
#define RAND_MAX                  32767
#endif

#define M1                        259200
#define IA1                       7141
#define IC1                       54773
#define RM1                       (1.0/M1)
#define M2                        134456
#define IA2                       8121
#define IC2                       28411
#define RM2                       (1.0/M2)
#define M3                        243000
#define IA3                       4561
#define IC3                       51349

double UniDev(void)
  
  /*
   * UniDev - Routine is adapted from the routine ran1() in "Numerical
   *          Recipes in C". It returns a uniform random deviate between
   *          0.0 and 1.0 and is entirely portable.
   *
   * Return - randev : the random number.
   */
  
{
  static long       ix1, ix2, ix3;
  static double     prevrans[MAX_STORED];
  static int        iff = 0;
  int               j;
  double            randev;
  struct timeval    tv;
  struct timezone   tz;
  
  if (iff == 0)
  {
    // time(NULL) returns seconds since some day
    if (!UniDevSeed)
    {
      gettimeofday(&tv, &tz);
      UniDevSeed = (unsigned)(tv.tv_usec % RAND_MAX);
      srand(UniDevSeed);
      UniDevSeed = (unsigned)(rand() % RAND_MAX);
    }
    ix1 = (IC1 - UniDevSeed) % M1;
    ix1 = (IA1*ix1 + IC1) % M1;
    ix2 = ix1 % M2;
    ix1 = (IA1*ix1 + IC1) % M1;
    ix3 = ix1 % M3;
    for (j = 0; j < MAX_STORED; j++)
    {
      ix1         = (IA1*ix1 + IC1) % M1;
      ix2         = (IA2*ix2 + IC2) % M2;
      prevrans[j] = (ix1 + ix2*RM2)*RM1;
    }
    iff = 1;
  }
  
  ix1 = (IA1*ix1 + IC1) % M1;
  ix2 = (IA2*ix2 + IC2) % M2;
  ix3 = (IA3*ix3 + IC3) % M3;
  j   = (int)((MAX_STORED*ix3)/(double)M3);
  
  if ((j < 0) || (j >= MAX_STORED)) ErrorExit(1, "UniDev: Invalid index.");
  
  randev      = prevrans[j];
  prevrans[j] = (ix1 + ix2*RM2)*RM1;
  
  return randev;
}

#undef M1
#undef IA1
#undef IC1
#undef RM1
#undef M2
#undef IA2
#undef IC2
#undef RM2
#undef M3
#undef IA3
#undef IC3


/*==================================================================================================================================*/
