/* 
	Title: modeltest
	Programmer: David Posada
	Date started: 03/05/98
	Purpose: Evaluate the fit of several models of evolution to a given data and unrooted tree.
	Compare different models of DNA substitution using likelihood ratio test (and a chi-square distribution)
	or/and the AIC criterion
	
		
	COPYRIGHT
	--------
	Copyright (C) 1998-2005 David Posada
	Facultad de Biolog�a, Universidad de Vigo, Campus Universitario 36310 Vigo, Spain	
	dposada@uvigo.es

	This program is free software; you can redistribute it and/or
	modify it under the terms of the GNU General Public License
	as published by the Free Software Foundation; either version 2
	of the License, or (at your option) any later version.

	This program is distributed in the hope that it will be useful,
	but WITHOUT ANY WARRANTY; without even the implied warranty of
	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
	GNU General Public License for more details.

	You should have received a copy of the GNU General Public License
	along with this program; if not, write to the Free Software
	Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.

	
	HISTORY
	-------
	Version 2.0 (June 99):    		Models TrN, K81 and its equivalents with equal base frequencies added

	Version 2.1 (October 99): 		Changed (-1) #free parameters (eg. JC has 0 instead of 1).

	Version 3.0beta1 (November 99):	Models TIM, TIMef, TVM and TVMef added. Now we have 56 models.

	Version 3.0beta3 (January 00):	Now PAUP* outputs base frequencies estimates in the likelihood scores file, 
									so changes were done to read this new file

	Version 3.0 (February 00):		Several aesthetic changes

	Version 3.01(March 00):			The program ouputs now a block of PAUP commands
									to implement the likelihood settings for the best-fit model				

	Version 3.02 (June 00):			The frequencies of C and G were interchanged in the output. Fixed. Thanks to Carlos Lopez-Vaamonde

	Version 3.03 (June 00):			The mixed chi-square distribution is added for I and G tests

	Version 3.04 (July 00):			Several cosmetic changes

	Version 3.05 (Feb 01):			In the windows version, the AIC[55] gave an AIC of 0 to the GTRIG. Now AIC[56]. (Juan Suarez)
									TIM+G reported invariable sites instead of gamma shape (Cymon Cox)

	Version 3.06 (Apr 01):			Print likelihood scores by default
									In the windows version there was a bug by which the file scores.txt was always
									the standard input (Andy Vierstraete)
									Using GNU licencese (I should have done this a long time ago) (thanks to Naoki Takebayashi)

	Version 3.07 (Apr 01):			Mistake in the print likelihoods corrected
									Some model names were not displayed complete in the screen (e.g TrNef+I+G)

	Version 3.1 (March 02):			Calculate delta AIC and Akaike weights
									Removed option to turn off the use of mixed chi square for I and G LRTs

	Version 3.2 (March 03):			Aesthetic changes
									TrN+I had 5 df instead of 6

	Version 3.3 (Nov 03):			Option to include branch length estimates as parameters
									Option to calculate AICc
									Changing some option letters accordingly

	Version 3.4 (March 04):			There was a typo printing the Rd value for K81uf+I.
									It was printing p-inv instead. (thanks to Michael Sorenson)

	Version 3.5 (May 04):			This is a minor update that does not affect the calculations. AIC weights were sorted by 
									their value, but because these can be almost zero (zero for the computer) for several models, 
									their order would not make sense in the light of the AIC values. 
									Now the program order the AIC weights by the AIC scores. 

	Version 3.6  (November 04):		The program includes now model averaged estimates
									The program includes now variable importance calculations
									New option (-w) to define confidence interval of models used for model-averaged estimates. By
										default this interval is 1.0, so all models are included in model-averaged estimates
									Use double precision in the AIC calculator (thanks to Renee Park)
									Likelihoods and number of parameters are now printed in the AIC table
									Aesthetic changes
									Argument for specifying sample size is now -n (it was -c)

	Version 3.7 (December 04):		Averaged estimate for alpha is now alpha(G) and alpha(IG), 
										instead of alpha(G) and alpha(G+IG)  (thanks to Roman Biek, now this
										is congruent with Posada and Buckley (2004))
									Now report confidence level for hLRTs
									Minor aesthetic changes

				(March 05):			Fixed cosmetic bug: In the comments of the AAC PAUP* block it was printing the name 
									of the hlRT model instead of the name of the AIC model

				(June 05):			Implemented BIC (option "-b").Changed some function and array names accordingly (AIC -> IC)
									Removed AICcalc function (option "-i)
									Removed AICfile	function(option "-f)
									
*/


#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <ctype.h>
#include <string.h>
#include <unistd.h>

#ifdef macintosh
	#include <unix.h>
	#include <sioux.h>
	#include <console.h>
	#include <io.h>
#endif


/* Constants */
#define	BIGX            20.0                                 /* max value to represent exp (x) */
#define	LOG_SQRT_PI     0.5723649429247000870717135          /* log (sqrt (pi)) */
#define	I_SQRT_PI       0.5641895835477562869480795          /* 1 / sqrt (pi) */
#define	Z_MAX           6.0                                  /* maximum meaningful z value */ 
#define	ex(x)           (((x) < -BIGX) ? 0.0 : exp (x))   
#define	MAX_PROB        0.999999   
#define	MIN_PROB        0.000001
#define	PROGRAM_NAME   "Modeltest"
#define	VERSION_NUMBER "3.7"
#define	SUCCESS         1
#define	FAILURE         0
#define	YES				1
#define	NO				0
#define	BIGNUMBER       9999999
#define	NUM_MODELS		56
#define	NA				-99999
	
/* Structures */
typedef struct {
float	ln;
int		parameters;        
char	*name;
}		ModelSt;

/* Prototypes */
static void		ParseArgs(int, char**);
static int		RecognizeInputFormat(); // change prototype by Minh
static int		ReadPaupScores(); // change prototype by Minh
static void		Initialize();
static int 	ReadScores(); // change prototype by Minh
static void 	PrintTitle(FILE *fp);
static void 	PrintDate(FILE *fp);
static void 	PrintOS (FILE *fp);
static void 	CheckExpirationDate (int month, int year, int different);
static void 	RatioCalc();
static void	 	hLRT();
static void 	CalculateIC();
static void 	InformationWeights();
static void 	SetModel(char[]);
static void 	Allocate();
/*static void 	Free();*/
void 	Free();
static void 	PrintUsage();
static void 	Output(char *selection, float value);
static void		WritePaupBlock(int ishLRT);
static double	LRT(ModelSt *model0, ModelSt *model1);
static double	LRTmix(ModelSt *model0, ModelSt *model1);
static void 	PrintRunSettings();
static void		ModelAveraging();
static void		AverageEstimates (char *parameter, int numModels, int modelIndex[], int estimateIndex[], 
						double *importance, double *averagedEstimate, double minWeightToAverage);
static double	FindMinWeightToAverage ();
static char 	*CheckNA (double value);

#ifdef macintosh
	static void Sioux();
#endif

float	ChiSquare (float x, int);
float	Normalz (float);                          
float 	TestEqualBaseFrequencies(ModelSt *, ModelSt *);
float	TestTiequalsTv(ModelSt *, ModelSt *);
float	TestEqualTiRates(ModelSt *, ModelSt *);
float	TestEqualTvRates(ModelSt *, ModelSt *);
float	TestTwoTvRates(ModelSt *, ModelSt *);
float	TestEqualSiteRates(ModelSt *, ModelSt *);
float	TestInvariableSites(ModelSt *, ModelSt *);

/* Global variables */
ModelSt *JC, *F81;
ModelSt *JC, *JCI, *JCG, *JCIG, *F81, *F81I, *F81G, *F81IG, *K80, *K80I, *K80G, *K80IG; 
ModelSt	*HKY, *HKYI, *HKYG, *HKYIG, *SYM, *SYMI, *SYMG, *SYMIG, *GTR, *GTRI, *GTRG, *GTRIG;
ModelSt	*TrN, *TrNI, *TrNG, *TrNIG, *TrNef, *TrNefI, *TrNefG, *TrNefIG, *K81, *K81I, *K81G, *K81IG,*K81uf, *K81ufI, *K81ufG, *K81ufIG;
ModelSt	*TVM, *TVMI, *TVMG, *TVMIG, *TVMef, *TVMefI,*TVMefG, *TVMefIG, *TIM, *TIMI, *TIMG, *TIMIG, *TIMef, *TIMefI,*TIMefG, *TIMefIG;
float 	score[488];
float  	ln[NUM_MODELS];
float	IC[NUM_MODELS];
float	wIC[NUM_MODELS];
int		orderedIC[NUM_MODELS]; /* keeps the index of the ordered (small to big) ICs */
float	alpha, minIC;        
int		format, file_id;                    
int		print_scores;
int		DEBUGLEVEL; 
char	*modelhLRT; 
char	*modelIC;   
char 	*nexttok1;
char 	*nexttok2;
char	filename[80];
char	infilename[80];
ModelSt	*model; 
ModelSt	*order;
FILE	*fpin;
int		mixchi;
int		numTaxa, numBL, sampleSize;
int		useAICc, useBL, useBIC;
int		lastModelConfidence;
float 	averagingConfidenceInterval;	
double	cumConfidenceWeight;


/* Parameter estimates for the selected model */
float fA, fC, fG, fT;
float TiTv;
float Ra, Rb, Rc, Rd, Re, Rf;
float shape;
float pinv; 
float theln;


/****************************** MAIN ***********************************/
int run_modeltest(int argc, char **argv)
{
	float start, secs;

	Allocate();
	alpha= 0.01;            /* default level of significance (aprox Bonferroni)  */
	mixchi = YES;			/* by default use mixed chi-square distribution */
	print_scores = NO;		/* by default print the likelihood scores for all models */
	numTaxa = 0;
	sampleSize = 0;
	useBL = NO; 			/* by default do not include branch length estimates as parameters */
	useAICc = NO;			/* by default do not use the AICc correction */
	useBIC = NO;			/* by default use AIC */
	averagingConfidenceInterval = 1.0; /* by default include all models in model-averaged estimates */
	
	/* CheckExpirationDate (int month, int year, int different);*/

	#ifdef macintosh
		_fcreator = 'R*ch';
		_ftype = 'TEXT';
		argc=ccommand(&argv);
		Sioux();
	#endif

	start = clock();
	PrintTitle(stdout);
	PrintDate(stdout);
	PrintOS(stdout);
	ParseArgs(argc, argv);
	
	file_id = isatty(fileno(stdin));
	if (file_id)
		{
		fprintf(stderr, "\n\nNo input file\n\n");
		PrintUsage();
		exit(1);
		}

	if (!RecognizeInputFormat()) return EXIT_FAILURE;
	PrintRunSettings();
	
	/* Do hLRTs */
	printf("\n\n\n\n---------------------------------------------------------------");
	printf("\n*                                                             *");
	printf("\n*         HIERARCHICAL LIKELIHOD RATIO TESTS (hLRTs)          *");
	printf("\n*                                                             *");
	printf("\n---------------------------------------------------------------\n");

	printf("\nConfidence level = %4.2f\n", alpha);

	hLRT();
	SetModel(modelhLRT);

	if (format == 0)
		{	
		if (shape > 999) /* alpha shape = infinity */
			{
			printf("\n\nWARNING: Although the model %s was initially selected, gamma (G) was removed ",modelhLRT);
			printf("because the estimated shape equals infinity, which implies equal rates among sites.");

			/* removing +G */
			nexttok1 = strtok(modelhLRT, "+");
			strcpy(modelhLRT, nexttok1);
			if (!strcmp(strtok(NULL, "+"),"I"))
				strcat(modelhLRT,"+I");
			}		
		Output (modelhLRT,0);
		WritePaupBlock(YES);
		}	
	else
		printf("\n hLRT model = %s", modelhLRT);

	/* Do IC */
	
	if (useBIC == NO)
		{
		if (useAICc == YES)
			{
			printf("\n\n\n\n\n---------------------------------------------------------------");
			printf("\n*                                                             *");
			printf("\n*        CORRECTED AKAIKE INFORMATION CRITERION (AICc)        *");
			printf("\n*                                                             *");
			printf("\n---------------------------------------------------------------\n");
			}
		else
			{
			printf("\n\n\n\n\n---------------------------------------------------------------");
			printf("\n*                                                             *");
			printf("\n*             AKAIKE INFORMATION CRITERION (AIC)              *");
			printf("\n*                                                             *");
			printf("\n---------------------------------------------------------------\n");
			}
		}
	else
		{
			printf("\n\n\n\n\n---------------------------------------------------------------");
			printf("\n*                                                             *");
			printf("\n*             BAYESIAN INFORMATION CRITERION (BIC)            *");
			printf("\n*                                                             *");
			printf("\n---------------------------------------------------------------\n");
		}	
		
	CalculateIC();
	SetModel(modelIC);

	if (format == 0)
		{	
		Output (modelIC,minIC);
		WritePaupBlock(NO);
		}
	else if (useAICc == YES)
		printf("\n AICc model = %s", modelIC);
	else	
		printf("\n AIC model = %s", modelIC);

	InformationWeights();
	ModelAveraging();
	/*Free(); remove this by Minh */

	secs = (double)(clock() - start) / CLOCKS_PER_SEC;
	printf("\n\n_________________________________________________________________");
	printf("\nProgram is done.");
	printf("\nTime processing: %G seconds", secs);
	printf("\nIf you need help type '-?' in the command line of the program.\n\n");

	return EXIT_SUCCESS;
}


/******************** PrintRunSettings **************************/
static void PrintRunSettings()
{
	int	i;
	
	/* Check settings */
	if (useBIC == YES && sampleSize == 0)
		{
		fprintf(stderr,"\n\n\nERROR: If you want to use the BIC (-b), you need to specify also sample size (-n)\n\n");
	    PrintUsage();
	    exit(-1);
		}

	printf("\n\nRun settings:");
	
	if (useBIC == YES)
		{
		if (sampleSize <= 0)
			{
			fprintf(stderr,"\n\nIf you want to use the BIC (-b), you need to specify also sample size (-n)\n");
		    PrintUsage();
		    exit(-1);
			}
		if (sampleSize <= model[NUM_MODELS-1].parameters)
			{
			fprintf (stdout, "\n\nYou have more parameters than data for some models!");
			fprintf (stdout, "\nExiting the program ...\n");
			exit(1);
			}

		fprintf (stdout, "\n  Using the BIC");
		fprintf (stdout, "\n    sample size = %d", sampleSize);
		}
	else
		{
		if (sampleSize > 0)
			{
			useAICc = YES;	
			/* Check the data is large enough*/
			if (sampleSize <= model[NUM_MODELS-1].parameters)
				{
				fprintf (stdout, "\n\nYou have more parameters than data for some models!");
				fprintf (stdout, "\nCalculations cannot be performed. Exiting the program ...\n");
				exit(1);
				}
			fprintf (stdout, "\n  Using the AICc correction");
			fprintf (stdout, "\n    sample size = %d", sampleSize);
			}
		else
			{
			useAICc = NO;	
			fprintf (stdout, "\n  Using the standard AIC (not the AICc)");
			}
		}
		
	if (numTaxa > 0)
		{
		useBL = YES;
		fprintf (stdout, "\n  Using branch lengths as parameters");
		fprintf (stdout, "\n    number of taxa = %d (%d branch lengths) have been estimated", numTaxa, numBL);

		for (i=0; i<NUM_MODELS; i++)
			{
			model[i].parameters += numBL;
			order[i].parameters += numBL;
			}
		}
	else
		{
		useBL = NO;
		fprintf (stdout, "\n  Not using branch lengths as parameters");
		}

	if (averagingConfidenceInterval == 1)
		fprintf (stdout, "\n  Including all models in model-averaging calculations");
	else 
		fprintf (stdout, "\n  Including only models within the %4.2f confidence interval", averagingConfidenceInterval);
			

		
}

/******************** ParseArgs **************************/
static void ParseArgs(int argc,char **argv)
{
	int i;
	char flag; 
 
	for (i=1; i<argc; i++)
	{
		argv[i]++;
		flag=*argv[i];
		argv[i]++;
		 
	    switch (flag)
	      {
	      case 'd':
	        DEBUGLEVEL = atoi(argv[i]);
	        if (DEBUGLEVEL >= 2)
	          fprintf(stderr,"DEBUGLEVEL set to %d\n",
	            DEBUGLEVEL);
	        break;

	      case 'a':
	        alpha = atof(argv[i]);
	        fprintf(stderr, "alpha level set to %.6f\n", 
	          alpha);
	        break;

	      case 'n':
	        sampleSize = atoi(argv[i]);
	        /*fprintf(stderr, "number of characters set to %d\n", sampleSize);*/
	        break;

	      case 't':
	        numTaxa = atoi(argv[i]);
	        numBL = 2 * numTaxa - 3;
	        /*fprintf(stderr, "number of taxa set to %d\n", numTaxa);*/
	        break;

	      case 'w':
	        averagingConfidenceInterval = atof(argv[i]);
	        if (averagingConfidenceInterval > 1)
	        	{
	        	fprintf (stderr, "\nError: confidence interval cannot be > 1");
	        	exit (1);
	        	}
	        else if (averagingConfidenceInterval <= 0)
	        	{
	        	fprintf (stderr, "\nError: confidence interval cannot be <= 0");
	        	exit (1);
	        	}
	        break;
   
	      case 'l':
	         printf("\n LRT CALCULATOR MODE \n");
	         RatioCalc();
	         break;
	
	      case 'b':
	        useBIC = YES;
	        break;

	      case '?': 
	        PrintUsage();
	        exit(1);
	        break;

	      case 'h': 
	        PrintUsage();
	        exit(1);
	        break;
	
	      default:
	        fprintf(stderr,"Unknown argument on the command line '%c'\n",flag);
	        PrintUsage();
	        exit(1);
	        break;
		}  
	} 
}


/********************* RecognizeInputFormat ***********************/
static int RecognizeInputFormat()
{
	int iochar;
		
	if (stdin == NULL)   
		{
		fprintf(stderr,"Error opening the input file");
		exit(0);
		}
		
	iochar=getc(stdin);
	
	if (iochar == (int)'T')    /* In the Paup matrix, in the first line there is the word 'Tree'*/
		{
		ungetc(iochar,stdin);
		printf("\nInput format: PAUP* scores file \n");
		format=0;
		if (!ReadPaupScores()) return 0;
		}
	else 
		{
		ungetc(iochar,stdin);
		printf("\nInput format: raw log likelihood scores \n");
		format=1;
		if (!ReadScores()) return 0;
		}
	return 1;
}


/***************************** ReadPaupScores ********************************/
static int ReadPaupScores()
{
	int 	iochar;
	int 	i,j,k;
	char 	string [120];
	i=0;

	while (!feof(stdin))
		{
		iochar=	getc(stdin);

		if (isdigit(iochar))    
			{
			ungetc(iochar,stdin);
			scanf("%f", &score[i]);
			if (DEBUGLEVEL >= 2)
				fprintf(stdout,"\nINFO:   Storing %f in score[%d]", score[i],i);
			i++;
			}
				
		if (isalpha (iochar))    
			{
			ungetc(iochar,stdin);
			scanf("%s",string);
			if (DEBUGLEVEL >= 2)
				fprintf(stdout,"\nINFO:   Reading string %s", string);			
			if (strcmp(string,"infinity") == 0)
				{
				score[i]= 999.999;
				if (DEBUGLEVEL >=2)
				fprintf(stdout,"INFO:   Storing %f in score[%d]", score[i],i);
				i++;
				}
			}
		}		

	if (ferror(stdin))
		{ 
		perror ("Modeltest");
		clearerr(stdin);
		}

	Initialize();

	if(print_scores == YES)
		{
		printf("\n\n  ** Log Likelihood scores **");
		printf("\n%-12.12s\t\t\t\t+I\t\t\t+G\t\t\t+I+G", " ");
		for(k=0; k<NUM_MODELS; k+=4)
			printf("\n%-10.10s =\t%9.4f\t%9.4f\t%9.4f\t%9.4f", model[k].name, model[k].ln, model[k+1].ln, model[k+2].ln, model[k+3].ln);
		printf("\n\n");
		}

		/*for(k=0; k<NUM_MODELS; k++)
			printf("\n%s", model[k].name);*/
		/*for(k=0; k<=NUM_MODELS; k++)
			printf("\n%12.8f", model[k].ln);*/
		/*for(k=0; k<=NUM_MODELS; k++)
			printf("\n%d", model[k].parameters);*/

	for (j=0; j<NUM_MODELS; j++)
		{	
		if (model[j].ln == 0 || i<487)
			{
			printf("\n\nThe input file is incomplete or incorrect. \nAre you using the most updated block of PAUP* commands?. ");
			printf("See that beginning in version 3.x Modeltest reads 56 instead of 40 (v2.x) or 24 (v1.x) models, so you should use the appropriate PAUP* block\n");
			printf("\nThis version of Modeltest is not compatible with versions of PAUP* older than PAUP*4.0beta3");
			printf("\nCheck the Modeltest and PAUP* web pages");
			return 0;
			}
		}	
	return 1;
}



/************** Initialize **********************/
void Initialize()
{
	
	JC=			model;		model[0].ln = order[0].ln = score[1];
	JCI= 		model + 1;	model[1].ln = order[1].ln = score[3];
	JCG= 		model + 2;	model[2].ln = order[2].ln = score[6];
	JCIG=		model + 3;	model[3].ln = order[3].ln = score[9];
	F81=		model + 4;	model[4].ln = order[4].ln = score[13];
	F81I=		model + 5;	model[5].ln = order[5].ln = score[19];
	F81G=		model + 6;	model[6].ln = order[6].ln = score[26];
	F81IG=		model + 7;	model[7].ln = order[7].ln = score[33];
	K80=		model + 8;	model[8].ln = order[8].ln = score[41];
	K80I=		model + 9;	model[9].ln = order[9].ln = score[44];
	K80G=		model + 10;	model[10].ln= order[10].ln= score[48];
	K80IG=		model + 11;	model[11].ln= order[11].ln= score[52];
	HKY=		model + 12;	model[12].ln= order[12].ln= score[57];
	HKYI=		model + 13;	model[13].ln= order[13].ln= score[64];
	HKYG=		model + 14;	model[14].ln= order[14].ln= score[72];
	HKYIG=		model + 15;	model[15].ln= order[15].ln= score[80];
	TrNef=		model + 16;	model[16].ln= order[16].ln= score[89];
	TrNefI=		model + 17;	model[17].ln= order[17].ln= score[96];
	TrNefG=		model + 18;	model[18].ln= order[18].ln= score[104];
	TrNefIG=	model + 19;	model[19].ln= order[19].ln= score[112];	
	TrN=		model + 20;	model[20].ln= order[20].ln= score[121];
	TrNI=		model + 21;	model[21].ln= order[21].ln= score[132];
	TrNG=		model + 22;	model[22].ln= order[22].ln= score[144];
	TrNIG=		model + 23;	model[23].ln= order[23].ln= score[156];	
	K81=		model + 24;	model[24].ln= order[24].ln= score[169];
	K81I=		model + 25;	model[25].ln= order[25].ln= score[176];
	K81G=		model + 26;	model[26].ln= order[26].ln= score[184];
	K81IG=		model + 27;	model[27].ln= order[27].ln= score[192];
	K81uf=		model + 28;	model[28].ln= order[28].ln= score[201];
	K81ufI=		model + 29;	model[29].ln= order[29].ln= score[212];
	K81ufG=		model + 30;	model[30].ln= order[30].ln= score[224];
	K81ufIG=	model + 31;	model[31].ln= order[31].ln= score[236];
	TIMef=		model + 32;	model[32].ln= order[32].ln= score[249];
	TIMefI=		model + 33;	model[33].ln= order[33].ln= score[256];
	TIMefG=		model + 34;	model[34].ln= order[34].ln= score[264];
	TIMefIG=	model + 35;	model[35].ln= order[35].ln= score[272];
	TIM=		model + 36;	model[36].ln= order[36].ln= score[281];
	TIMI=		model + 37;	model[37].ln= order[37].ln= score[292];
	TIMG=		model + 38;	model[38].ln= order[38].ln= score[304];
	TIMIG=		model + 39;	model[39].ln= order[39].ln= score[316];
	TVMef=		model + 40;	model[40].ln= order[40].ln= score[329];
	TVMefI=		model + 41;	model[41].ln= order[41].ln= score[336];
	TVMefG=		model + 42;	model[42].ln= order[42].ln= score[344];
	TVMefIG=	model + 43;	model[43].ln= order[43].ln= score[352];
	TVM=		model + 44;	model[44].ln= order[44].ln= score[361];
	TVMI=		model + 45;	model[45].ln= order[45].ln= score[372];
	TVMG=		model + 46;	model[46].ln= order[46].ln= score[384];
	TVMIG=		model + 47;	model[47].ln= order[47].ln= score[396];
	SYM=		model + 48;	model[48].ln= order[48].ln= score[409];
	SYMI=		model + 49;	model[49].ln= order[49].ln= score[416];
	SYMG=		model + 50;	model[50].ln= order[50].ln= score[424];
	SYMIG=		model + 51;	model[51].ln= order[51].ln= score[432];
	GTR=		model + 52;	model[52].ln= order[52].ln= score[441];
	GTRI=		model + 53;	model[53].ln= order[53].ln= score[452];
	GTRG=		model + 54;	model[54].ln= order[54].ln= score[464];
	GTRIG=		model + 55;	model[55].ln= order[55].ln= score[476];

   
	/* free parameters */
	
	/* JC */
	model[0].parameters= order[0].parameters= 0;
	model[1].parameters= order[1].parameters= 1;
	model[2].parameters= order[2].parameters= 1;
	model[3].parameters= order[3].parameters= 2;

	/* F81*/
	model[4].parameters= order[4].parameters= 3;
	model[5].parameters= order[5].parameters= 4;
	model[6].parameters= order[6].parameters= 4;
	model[7].parameters= order[7].parameters= 5;

	/* K80 */
	model[8].parameters= order[8].parameters= 1;
	model[9].parameters= order[9].parameters= 2;
	model[10].parameters= order[10].parameters= 2;
	model[11].parameters= order[11].parameters= 3;

	/* HKY */
	model[12].parameters= order[12].parameters= 4;
	model[13].parameters= order[13].parameters= 5;
	model[14].parameters= order[14].parameters= 5;
	model[15].parameters= order[15].parameters= 6;

	/* TrNef */
	model[16].parameters= order[16].parameters= 2;
	model[17].parameters= order[17].parameters= 3;
	model[18].parameters= order[18].parameters= 3;
	model[19].parameters= order[19].parameters= 4;

	/* TrN */
	model[20].parameters= order[20].parameters= 5;
	model[21].parameters= order[21].parameters= 6;
	model[22].parameters= order[22].parameters= 6;
	model[23].parameters= order[23].parameters= 7;
	
	/* K81 */
	model[24].parameters= order[24].parameters= 2;
	model[25].parameters= order[25].parameters= 3;
	model[26].parameters= order[26].parameters= 3;
	model[27].parameters= order[27].parameters= 4;

	/* K81uf*/
	model[28].parameters= order[28].parameters= 5;
	model[29].parameters= order[29].parameters= 6;
	model[30].parameters= order[30].parameters= 6;
	model[31].parameters= order[31].parameters= 7;

	/* TIMef */
	model[32].parameters= order[32].parameters= 3;
	model[33].parameters= order[33].parameters= 4;
	model[34].parameters= order[34].parameters= 4;
	model[35].parameters= order[35].parameters= 5;

	/* TIM */
	model[36].parameters= order[36].parameters= 6;
	model[37].parameters= order[37].parameters= 7;
	model[38].parameters= order[38].parameters= 7;
	model[39].parameters= order[39].parameters= 8;

	/* TVMef */
	model[40].parameters= order[40].parameters= 4;
	model[41].parameters= order[41].parameters= 5;
	model[42].parameters= order[42].parameters= 5;
	model[43].parameters= order[43].parameters= 6;

	/* TVM */
	model[44].parameters= order[44].parameters= 7;
	model[45].parameters= order[45].parameters= 8;
	model[46].parameters= order[46].parameters= 8;
	model[47].parameters= order[47].parameters= 9;

	/* SYM */
	model[48].parameters= order[48].parameters= 5;
	model[49].parameters= order[49].parameters= 6;
	model[50].parameters= order[50].parameters= 6;
	model[51].parameters= order[51].parameters= 7;

	/* GTR */
	model[52].parameters= order[52].parameters= 8;
	model[53].parameters= order[53].parameters= 9;
	model[54].parameters= order[54].parameters= 9;
	model[55].parameters= order[55].parameters= 10;


	order[0].name=  "JC";
	order[1].name=  "JC+I";
	order[2].name=  "JC+G";
	order[3].name=  "JC+I+G";
	order[4].name=  "F81";
	order[5].name=  "F81+I";
	order[6].name=  "F81+G";
	order[7].name=  "F81+I+G";
	order[8].name=  "K80";
	order[9].name=  "K80+I";
	order[10].name= "K80+G";
	order[11].name= "K80+I+G";
	order[12].name= "HKY";
	order[13].name= "HKY+I";
	order[14].name= "HKY+G";
	order[15].name= "HKY+I+G";
	order[16].name= "TrNef";
	order[17].name= "TrNef+I";
	order[18].name= "TrNef+G";
	order[19].name= "TrNef+I+G";
	order[20].name= "TrN";
	order[21].name= "TrN+I";
	order[22].name= "TrN+G";
	order[23].name= "TrN+I+G";
	order[24].name= "K81";
	order[25].name= "K81+I";
	order[26].name= "K81+G";
	order[27].name= "K81+I+G";
	order[28].name= "K81uf";
	order[29].name= "K81uf+I";
	order[30].name= "K81uf+G";
	order[31].name= "K81uf+I+G";
	order[32].name= "TIMef";
	order[33].name= "TIMef+I";
	order[34].name= "TIMef+G";
	order[35].name= "TIMef+I+G";
	order[36].name= "TIM";
	order[37].name= "TIM+I";
	order[38].name= "TIM+G";
	order[39].name= "TIM+I+G";
	order[40].name= "TVMef";
	order[41].name= "TVMef+I";
	order[42].name= "TVMef+G";
	order[43].name= "TVMef+I+G";
	order[44].name= "TVM";
	order[45].name= "TVM+I";
	order[46].name= "TVM+G";
	order[47].name= "TVM+I+G";
	order[48].name= "SYM";
	order[49].name= "SYM+I";
	order[50].name= "SYM+G";
	order[51].name= "SYM+I+G";
	order[52].name= "GTR";
	order[53].name= "GTR+I";
	order[54].name= "GTR+G";
	order[55].name= "GTR+I+G";

	model[0].name=  "JC";
	model[1].name=  "JC+I";
	model[2].name=  "JC+G";
	model[3].name=  "JC+I+G";
	model[4].name=  "F81";
	model[5].name=  "F81+I";
	model[6].name=  "F81+G";
	model[7].name=  "F81+I+G";
	model[8].name=  "K80";
	model[9].name=  "K80+I";
	model[10].name= "K80+G";
	model[11].name= "K80+I+G";
	model[12].name= "HKY";
	model[13].name= "HKY+I";
	model[14].name= "HKY+G";
	model[15].name= "HKY+I+G";
	model[16].name= "TrNef";
	model[17].name= "TrNef+I";
	model[18].name= "TrNef+G";
	model[19].name= "TrNef+I+G";
	model[20].name= "TrN";
	model[21].name= "TrN+I";
	model[22].name= "TrN+G";
	model[23].name= "TrN+I+G";
	model[24].name= "K81";
	model[25].name= "K81+I";
	model[26].name= "K81+G";
	model[27].name= "K81+I+G";
	model[28].name= "K81uf";
	model[29].name= "K81uf+I";
	model[30].name= "K81uf+G";
	model[31].name= "K81uf+I+G";
	model[32].name= "TIMef";
	model[33].name= "TIMef+I";
	model[34].name= "TIMef+G";
	model[35].name= "TIMef+I+G";
	model[36].name= "TIM";
	model[37].name= "TIM+I";
	model[38].name= "TIM+G";
	model[39].name= "TIM+I+G";
	model[40].name= "TVMef";
	model[41].name= "TVMef+I";
	model[42].name= "TVMef+G";
	model[43].name= "TVMef+I+G";
	model[44].name= "TVM";
	model[45].name= "TVM+I";
	model[46].name= "TVM+G";
	model[47].name= "TVM+I+G";
	model[48].name= "SYM";
	model[49].name= "SYM+I";
	model[50].name= "SYM+G";
	model[51].name= "SYM+I+G";
	model[52].name= "GTR";
	model[53].name= "GTR+I";
	model[54].name= "GTR+G";
	model[55].name= "GTR+I+G";
}




/******************* ReadScores ************************/
static int ReadScores()
{
	float score[NUM_MODELS];
	int i,j;

	i=0;
	score[NUM_MODELS-1]= 0;

	while (!feof(stdin) && i < NUM_MODELS)
		{
		int items = scanf("%f", &score[i]);
		if (items < 1 || items == EOF) return 0;
		i++;
		}
	
	Initialize();



	for (j=0; j<NUM_MODELS; j++)
		{	
		if (model[j].ln == 0 || i<487)
			{
		printf("\n\nThe input file is incomplete or incorrect. \nAre you using the most updated block of PAUP* commands?. ");
		printf("See that beginning in version 3.x Modeltest reads 56 instead of 40 (v2.x) or 24 (v1.x) models, so you should use the appropriate PAUP* block\n");
		return 0;
			}
		}
		
	if (i>NUM_MODELS+1)
		{
		printf("\n\nThe input file has more than %d scores", NUM_MODELS);
		return 0;
		}
	return 1;
}




/******************* LRT ******************************/
/* peforms a likelihood ratio test */

static double LRT(ModelSt *model0, ModelSt *model1)
{
	double delta;
	double prob;
	int df;
	
	delta = 2 * (model0->ln - model1->ln);
	df = model1->parameters - model0->parameters;
	
	if (delta == 0)
		prob = 1.0;
	else
		prob= ChiSquare(delta, df);

	printf("\n   Null model = %-9.9s\t\t\t  -lnL0 = %.4f", model0->name, model0->ln);
	printf("\n   Alternative model = %-9.9s\t  -lnL1 = %.4f", model1->name, model1->ln);
	printf("\n   2(lnL1-lnL0) = %9.4f\t\t      df = %d ", delta, df);


	if (prob == 1.0)
		printf("\n   P-value = >%f", MAX_PROB);
	else if (prob < 0.000001)
		printf("\n   P-value = <%f", MIN_PROB);
	else 
		printf("\n   P-value =  %f", prob);

		return prob;
}


/******************* LRTmix ******************************/
/* peforms a likelihood ratio test and uses mixed chi2*/

static double LRTmix(ModelSt *model0, ModelSt *model1)
{
	double delta, prob;
	int df;
	
	delta = 2 * (model0->ln - model1->ln);
	df = model1->parameters - model0->parameters;

	if (delta == 0)
		prob = 1.0;
	else
		{
		if (df == 1)
			prob = ChiSquare(delta,df)/2;
		else	
			prob= (ChiSquare(delta,df-1) + ChiSquare(delta,df)) / 2;
		}
	printf("\n   Null model = %-9.9s\t\t\t  -lnL0 = %.4f", model0->name, model0->ln);
	printf("\n   Alternative model = %-9.9s\t  -lnL1 = %.4f", model1->name, model1->ln);
	printf("\n   2(lnL1-lnL0) = %9.4f\t\t      df = %d ", delta, df);
	printf("\n   Using mixed chi-square distribution");

	if (prob == 1.0)
		printf("\n   P-value = >%f", MAX_PROB);
	else if (prob < 0.000001)
		printf("\n   P-value = <%f", MIN_PROB);
	else 
		printf("\n   P-value =  %f", prob);

		return prob;
}


/******************* TestEqualBaseFrequencies ****************/
float TestEqualBaseFrequencies(ModelSt *model0, ModelSt *model1)
{
	float P;
	printf("\n Equal base frequencies");
	P = LRT(model0, model1);
	return P;
}


/*******************  TestTiequalsTv  **********************/
float TestTiequalsTv(ModelSt *model0, ModelSt *model1)
{
	float P;
	printf("\n Ti=Tv");
	P = LRT(model0, model1);
	return P;
	
}

/******************* TestEqualTiRates *******************/
float TestEqualTiRates (ModelSt *model0, ModelSt *model1)
{
	float P;
	printf("\n Equal Ti rates");
	P = LRT(model0, model1);
	return P;
	
}


/******************* TestEqualTvRates *******************/
float TestEqualTvRates (ModelSt *model0, ModelSt *model1)
{
	float P;
	printf("\n Equal Tv rates");
	P = LRT(model0, model1);
	return P;
}


/******************* TestTwoTvRates *******************/
float TestTwoTvRates (ModelSt *model0, ModelSt *model1)
{
	float P;
	printf("\n Only two Tv rates");
	P = LRT(model0, model1);
	return P;
}


/********************* TestEqualSiteRates **********************/
float TestEqualSiteRates (ModelSt *model0, ModelSt *model1)
{
	float P;
	printf("\n Equal rates among sites");
	if (mixchi)
		P = LRTmix(model0, model1);
	else
		P = LRT(model0, model1);
	return P;
}


/******************** TestInvariableSites **********************/
float TestInvariableSites (ModelSt *model0, ModelSt *model1)
{
	float P;
	printf("\n No Invariable sites");
	if (mixchi)
		P = LRTmix(model0, model1);
	else
		P = LRT(model0, model1);
	return P;	
}



/**************  ChiSquare: probability of chi square value *************/
/*ALGORITHM Compute probability of chi square value.
Adapted from: 	Hill, I. D. and Pike, M. C.  Algorithm 299.Collected Algorithms for the CACM 1967 p. 243
Updated for rounding errors based on remark inACM TOMS June 1985, page 185. Found in Perlman.lib*/
	
float ChiSquare (float x, int df)  /* x: obtained chi-square value,  df: degrees of freedom */
{
	float	a, y, s;
	float	e, c, z;
	int 	even;         /* true if df is an even number */
	
	if (x <= 0.0 || df < 1)
		return (1.0);
	
	y= 1;
	
	a = 0.5 * x;
	even = (2*(df/2)) == df;
	if (df > 1)
		y = ex (-a);
	s = (even ? y : (2.0 * Normalz (-sqrt(x))));
	if (df > 2)
		{
		x = 0.5 * (df - 1.0);
		z = (even ? 1.0 : 0.5);
		if (a > BIGX)
			{
   		e = (even ? 0.0 : LOG_SQRT_PI);
			c = log (a);
			while (z <= x)
				{
				e = log (z) + e;
				s += ex (c*z-a-e);
				z += 1.0;
				}
			return (s);
			}
		else
			{
			e = (even ? 1.0 : (I_SQRT_PI / sqrt (a)));
			c = 0.0;
			while (z <= x)
				{
				e = e * (a / z);
				c = c + e;
				z += 1.0;
				}
			return (c * y + s);
		}
	}
	else
		return (s);
}

	
/************** Normalz: probability of normal z value *********************/
/*
ALGORITHM:	Adapted from a polynomial approximation in:
			Ibbetson D, Algorithm 209
			Collected Algorithms of the CACM 1963 p. 616
		Note:
			This routine has six digit accuracy, so it is only useful for absolute
			z values < 6.  For z values >= to 6.0, Normalz() returns 0.0.
*/

float Normalz (float z)        /*VAR returns cumulative probability from -oo to z VAR normal z value */
{	float 	y, x, w;
	
	if (z == 0.0)
		x = 0.0;
	else
		{
		y = 0.5 * fabs (z);
		if (y >= (Z_MAX * 0.5))
			x = 1.0;
		else if (y < 1.0)
			{
			w = y*y;
			x = ((((((((0.000124818987 * w
				-0.001075204047) * w +0.005198775019) * w
				-0.019198292004) * w +0.059054035642) * w
				-0.151968751364) * w +0.319152932694) * w
				-0.531923007300) * w +0.797884560593) * y * 2.0;
			}
		else
			{
			y -= 2.0;
			x = (((((((((((((-0.000045255659 * y
				+0.000152529290) * y -0.000019538132) * y
				-0.000676904986) * y +0.001390604284) * y
				-0.000794620820) * y -0.002034254874) * y
				+0.006549791214) * y -0.010557625006) * y
				+0.011630447319) * y -0.009279453341) * y
				+0.005353579108) * y -0.002141268741) * y
					+0.000535310849) * y +0.999936657524;
				}
			}
	return (z > 0.0 ? ((x + 1.0) * 0.5) : ((1.0 - x) * 0.5));	
}


/************** RatioCalc *********************/
/*
	Interactive calculation of LRTs
*/
static void RatioCalc()
{
	float score1, score2, ratio, prob;

	int df;

	printf("\nPlease, input the POSITIVE log likelihood score corresponding to \nthe null model> ");
	scanf("%f", &score1);

	while (score1 < 0)
		{
		printf("\nBad Input: the program doesn't accept negative likelihood scores");
		printf("\n\nPlease, input the POSITIVE log likelihood score corresponding to \nthe null model> ");
		scanf("%f", &score1);
		}

	printf("\nPlease, input the POSITIVE log likelihood score corresponding to \nthe alternative model> ");
	scanf("%f", &score2);

	while (score2 < 0)
		{
		printf("Bad Input: the program doesn't accept negative likelihood scores");
		printf("\nPlease, input the POSITIVE log likelihood score corresponding to \nthe alternative model> ");
		scanf("%f", &score2);
		}

	if (score1 < score2) 
		{
		printf("\n\nIncorrect input: the positive likelihood of the null model cannot be smaller than the positive likelihood of the alternative model.");
		printf("\nYou should enter the correct positive log likelihood scores again�\n");
		exit(0);
		}

	printf("\nPlease, input the number of degrees of freedom> ");
	scanf("%d", &df);

	while (df < 1)
		{
		printf("\nThe number of degrees of freedom should be at least 1");
		printf("\nPlease, input the number of degrees of freedom> ");
		scanf("%d", &df);
		}

	ratio = 2*(score1-score2);
	prob = ChiSquare(ratio,df);

	printf("\n\n_________________________ Results of Ratio Calculator _______________________\n");
	printf("\nThe ratio is %f ", ratio);
	printf("\n\nThe probability of observing this ratio likelihood test statistic under a correct null model is %f\n", prob); 
	if (prob < alpha)
		printf ("\nThis is significant at the alpha level of %.4f", alpha);
		else
		printf ("\nThis is not significant at the alpha level of %.4f\n", alpha);

	exit(0);
}



/*********************** CalculateIC ***************************/
/* Calculates the AIC or AICc value for each likelihood score */

void CalculateIC()
{
	float smallerIC;
	int i, K, n;

	n = sampleSize;
	
	for (i=0; i<NUM_MODELS; i++)
		{
		K = model[i].parameters;
		
		if (useBIC == YES)
			IC[i] = 2 * model[i].ln + K * log(n);	
		else
			{
			IC[i] =  2 * (model[i].ln  +  K);
			if (useAICc == YES)
				IC[i] += 2 * K * (K+1) / (double) (n-K-1);
			}
		}
	
	smallerIC = IC[0];
	for (i=1; i< NUM_MODELS; i++)
		if (IC[i] < smallerIC)
			smallerIC = IC[i];
			
	minIC = BIGNUMBER;

	for (i=NUM_MODELS-1; i>=0; i--)
		if (IC[i] == smallerIC)
			{
			strcpy(modelIC, model[i].name);
			minIC = IC[i];
			}
	
}

/*********************** InformationWeights ****************************/
/* Calculates delta IC and Information weights (w[i]) */ 

void InformationWeights ()
{
	int		i, j, sorted, pass;
	float	deltaIC[NUM_MODELS];
	float	ord[NUM_MODELS];
	float	sumExp, temp;
	float	cumWeight;
	
	sumExp = 0;
	
	for (i=0; i<NUM_MODELS; i++)
		{
		deltaIC[i] = IC[i] - minIC;
		sumExp += exp(-0.5*deltaIC[i]);
		}

	for (i=0; i<NUM_MODELS; i++)
		{
		wIC[i] = exp(-0.5*deltaIC[i]) / sumExp;
		ord[i] = IC[i];
		orderedIC[i] = i;
		}
		
/* Sort by IC score to print weights in order */
	sorted=NO; 
	pass=1; 
	while (sorted==NO)
		{
		sorted = YES;
		for (i=0; i < (NUM_MODELS-pass); i++)
			{ 
	 		if (ord[i] > ord [i+1])
				{
				temp = ord[i+1];
				ord[i+1]= ord[i];
				ord[i] = temp;

				temp = orderedIC[i+1];
				orderedIC[i+1] = orderedIC[i];
				orderedIC[i] = temp;

				sorted = NO; 
				}
			}	
		pass++;	
		} 

	if (useBIC == YES)
		{
		printf ("\n\n\n * MODEL SELECTION UNCERTAINTY : BIC Weights (approximate posterior probabilities)");
		printf ("\n\nModel             -lnL   K         BIC      delta     weight  cumWeight");
		}
	
	else
		{
		printf ("\n\n\n * MODEL SELECTION UNCERTAINTY : Akaike Weights");
		if (useAICc == NO)
			printf ("\n\nModel             -lnL   K         AIC      delta     weight  cumWeight");
		else
			printf ("\n\nModel             -lnL   K        AICc      delta     weight  cumWeight");
		}
		
	printf ("\n-----------------------------------------------------------------------");
	cumWeight = 0;
	for (i=0; i<NUM_MODELS; i++)
		{
		j=orderedIC[i];
		cumWeight += wIC[j];
		if (wIC[j] > 0.0001)
			printf("\n%-10s\t%10.4f\t%2d\t%10.4f\t%9.4f\t%8.4f\t%7.4f", model[j].name, model[j].ln, model[j].parameters, IC[j], deltaIC[j], wIC[j], cumWeight);
		else
			printf("\n%-10s\t%10.4f\t%2d\t%10.4f\t%9.4f\t%4.2e\t%7.4f", model[j].name, model[j].ln, model[j].parameters, IC[j], deltaIC[j], wIC[j], cumWeight);
		
		/*printf("\n%-10s\t%8.4f\t%8.4f\t\t%10.8f\t%10.8f", model[j].name, IC[j], deltaIC[j], wIC[j], cumWeight);*/
		}
	printf ("\n-----------------------------------------------------------------------");
	printf ("\n-lnL:\t\tNegative log likelihod");
	printf ("\n K:\t\t\tNumber of estimated parameters");
	printf ("\n IC:\t\tInformation Criterion");
	printf ("\n delta:\t\tInformation difference");
	printf ("\n weight:\tInformation weight");
	printf ("\n cumWeight:\tCumulative information weight");


}



/************** ModelAveraging **********************/
/*	Calculates the importance for different parameters 
	of the models (it is simply the sum of the Akaike
	weights for those models that include such parameter) 	
	and model averaged estimates
	
	This method is completely brute force (See Java version)

	Assumes TrN y TIM estimate only Rb, Re
		K81 estimates no R parameter
		TVM estimates only Ra, Rc, Rd
		GTR y SIM estimate ra, Rb, Rc, Rd, Re
*/

void ModelAveraging()
{
	double	minWeightToAverage;
	double	ifA, ifC, ifG, ifT, ititv, iRa, iRb, iRc, iRd, iRe, ipinvI, ialphaG, ipinvIG, ialphaIG;
	double	wfA, wfC, wfG, wfT, wtitv, wRa, wRb, wRc, wRd, wRe, wpinvI, walphaG, wpinvIG, walphaIG;

	/* which index for scores */
	int 	efA[] = {14,20,27,34,58,65,73,81,122,133,145,157,202,213,225,237,282,293,305,317,362,373,385,397,442,453,465,477};
	int 	efC[] = {15,21,28,35,59,66,74,82,123,134,146,158,203,214,226,238,283,294,306,318,363,374,386,398,443,454,466,478};
	int 	efG[] = {16,22,29,36,60,67,75,83,124,135,147,159,204,215,227,239,284,295,307,319,364,375,387,399,444,455,467,479};
	int 	efT[] = {17,23,30,37,61,68,76,84,125,136,148,160,205,216,228,240,285,296,308,320,365,376,388,400,445,456,468,480};
	int 	etitv[] = {42,45,49,53,62,69,77,85};
	int		eRa[] = {330,337,345,353,366,377,389,401,410,417,425,433,446,457,469,481};
	int		eRb[] = {91,98,106,114,127,138,150,162,251,258,266,274,287,298,310,322,411,418,426,434,447,458,470,482};
	int		eRc[] = {332,339,347,355,368,379,391,403,412,419,427,435,448,459,471,483};
	int		eRd[] = {333,340,348,356,369,380,392,404,413,420,428,436,449,460,472,484};
	int		eRe[] = {94,101,109,117,130,141,153,165,254,261,269,277,290,301,313,325,414,421,429,437,450,461,473,485};
	int		epinvI[] = {4,24,46,70,102,142,182,222,262,302,342,382,422,462};	
	int		ealphaG[] = {7,31,50,78,110,154,190,234,270,314,350,394,430,474};
/*	int 	epinvIG[] = {4,10,24,38,46,54,70,86,102,118,142,166,182,198,222,246,262,278,302,326,342,358,382,406,422,438,462,486};
	int		ealphaIG[] = {7,11,31,39,50,55,78,87,110,119,154,167,190,199,234,247,270,279,314,327,350,359,394,407,430,439,474,487};*/
	int 	epinvIG[] = {10,38,54,86,118,166,198,246,278,326,358,406,438,486};
	int		ealphaIG[] = {11,39,55,87,119,167,199,247,279,327,359,407,439,487};

	/* which index for models containing the parameter  */
	int 	mfA[] = {4,5,6,7,12,13,14,15,20,21,22,23,28,29,30,31,36,27,38,39,44,45,46,47,52,53,54,55}; 
	int 	mfC[] = {4,5,6,7,12,13,14,15,20,21,22,23,28,29,30,31,36,27,38,39,44,45,46,47,52,53,54,55};
	int 	mfG[] = {4,5,6,7,12,13,14,15,20,21,22,23,28,29,30,31,36,27,38,39,44,45,46,47,52,53,54,55};
	int 	mfT[] = {4,5,6,7,12,13,14,15,20,21,22,23,28,29,30,31,36,27,38,39,44,45,46,47,52,53,54,55};
	int 	mtitv[] = {8,9,10,11,12,13,14,15};
	int		mRa[] = {40,41,42,43,44,45,46,47,48,49,50,51,52,53,54,55};
	int		mRb[] = {16,17,18,19,20,21,22,23,32,33,34,35,36,37,38,39,48,49,50,51,52,53,54,55};
	int		mRc[] = {40,41,42,43,44,45,46,47,48,49,50,51,52,53,54,55};
	int		mRd[] = {40,41,42,43,44,45,46,47,48,49,50,51,52,53,54,55};
	int		mRe[] = {16,17,18,19,20,21,22,23,32,33,34,35,36,37,38,39,48,49,50,51,52,53,54,55};
	int		mpinvI[] = {1,5,9,13,17,21,25,29,33,37,41,45,49,53};	
	int		malphaG[] = {2,6,10,14,18,22,26,30,34,38,42,46,50,54};
/*	int 	mpinvIG[] = {1,3,5,7,9,11,13,15,17,19,21,23,25,27,29,31,33,35,37,39,41,43,45,47,49,51,53,55};
	int		malphaIG[] = {2,3,6,7,10,11,14,15,18,19,22,23,26,27,30,31,34,35,38,39,42,43,46,47,50,51,54,55};*/
	int 	mpinvIG[] =  {3,7,11,15,19,23,27,31,35,39,43,47,51,55};
	int		malphaIG[] = {3,7,11,15,19,23,27,31,35,39,43,47,51,55};

	ifA = ifC = ifG = ifT = ititv = iRa = iRb = iRc = iRd = iRe = ipinvI = ialphaG = ipinvIG = ialphaIG = 0;
	wfA = wfC = wfG = wfT = wtitv = wRa = wRb = wRc = wRd = wRe = wpinvI = walphaG = wpinvIG = walphaIG = 0;
	
	if (averagingConfidenceInterval < 1)
		minWeightToAverage = FindMinWeightToAverage ();
	else
		{
		minWeightToAverage = 0.0;
		cumConfidenceWeight = 1.0;
		}
	
	/* calculate importances and model-averaged estimates */
	AverageEstimates ("fA",28, mfA, efA, &ifA, &wfA, minWeightToAverage);
	AverageEstimates ("fC",28, mfC, efC, &ifC, &wfC, minWeightToAverage);
	AverageEstimates ("fG",28, mfG, efG, &ifG, &wfG, minWeightToAverage);
	AverageEstimates ("fT",28, mfT, efT, &ifT, &wfT, minWeightToAverage);
	AverageEstimates ("titv", 8, mtitv, etitv, &ititv, &wtitv, minWeightToAverage);
	AverageEstimates ("Ra",16, mRa, eRa, &iRa, &wRa, minWeightToAverage);
	AverageEstimates ("Rb",24, mRb, eRb, &iRb, &wRb, minWeightToAverage);
	AverageEstimates ("Rc",16, mRc, eRc, &iRc, &wRc, minWeightToAverage);
	AverageEstimates ("Rd",16, mRd, eRd, &iRd, &wRd, minWeightToAverage);
	AverageEstimates ("Re",24, mRe, eRe, &iRe, &wRe, minWeightToAverage);
	AverageEstimates ("pinv(I)",14, mpinvI, epinvI, &ipinvI, &wpinvI, minWeightToAverage);
	AverageEstimates ("alpha(G)",14, malphaG, ealphaG, &ialphaG, &walphaG, minWeightToAverage);
	AverageEstimates ("pinv(IG)",14, mpinvIG, epinvIG, &ipinvIG, &wpinvIG, minWeightToAverage);
	AverageEstimates ("alpha(IG)",14, malphaIG, ealphaIG, &ialphaIG, &walphaIG, minWeightToAverage);
	
	/* print results */

	if (useBIC == YES)
		printf ("\n\n\n\n* MODEL AVERAGING AND PARAMETER IMPORTANCE (using BIC weights)");
	else if (useAICc == YES)
		printf ("\n\n\n\n* MODEL AVERAGING AND PARAMETER IMPORTANCE (using AICc weights)");
	else
		printf ("\n\n\n\n* MODEL AVERAGING AND PARAMETER IMPORTANCE (using Akaike Weights)");


	if (averagingConfidenceInterval == 1)
		fprintf (stdout, "\n  Including all %d models", NUM_MODELS);
	else 
		{
		fprintf (stdout, "\n    Including only the best %d models within the aproximate %4.2f (%6.4f)\n    confidence interval", 
					lastModelConfidence+1, averagingConfidenceInterval, cumConfidenceWeight);
		fprintf (stdout, "\n      minimum weight to average is %6.4f", minWeightToAverage);
		fprintf (stdout, "\n      weights are reescaled by the interval cumulative weight (%6.4f)", cumConfidenceWeight);
		}

	printf ("\n\n                         Model-averaged");
	printf ("\nParameter   Importance        estimates");
	printf ("\n---------------------------------------");
	printf ("\nfA\t\t\t\t%6.4f\t\t%11s",ifA, CheckNA(wfA)); 
	printf ("\nfC\t\t\t\t%6.4f\t\t%11s",ifC, CheckNA(wfC)); 
	printf ("\nfG\t\t\t\t%6.4f\t\t%11s",ifG, CheckNA(wfG)); 
	printf ("\nfT\t\t\t\t%6.4f\t\t%11s",ifT, CheckNA(wfT)); 
	printf ("\nTiTv\t\t\t%6.4f\t\t%11s",ititv, CheckNA(wtitv)); 
	printf ("\nrAC\t\t\t\t%6.4f\t\t%11s",iRa, CheckNA(wRa)); 
	printf ("\nrAG\t\t\t\t%6.4f\t\t%11s",iRb, CheckNA(wRb)); 
	printf ("\nrAT\t\t\t\t%6.4f\t\t%11s",iRc, CheckNA(wRc)); 
	printf ("\nrCG\t\t\t\t%6.4f\t\t%11s",iRd, CheckNA(wRd)); 
	printf ("\nrCT\t\t\t\t%6.4f\t\t%11s",iRe, CheckNA(wRe)); 
	printf ("\npinv(I)\t\t\t%6.4f\t\t%11s",ipinvI, CheckNA(wpinvI)); 
	printf ("\nalpha(G)\t\t%6.4f\t\t%11s",ialphaG, CheckNA(walphaG)); 
	printf ("\npinv(IG)\t\t%6.4f\t\t%11s",ipinvIG, CheckNA(wpinvIG)); 
	printf ("\nalpha(IG)\t\t%6.4f\t\t%11s",ialphaIG, CheckNA(walphaIG)); 
	printf ("\n---------------------------------------");

	printf ("\n Values have been rounded.");
/*	printf ("\n fA:\t\tfrecuency of adenine");
	printf ("\n fC:\t\tfrecuency of cytosine");
	printf ("\n fG:\t\tfrecuency of guanine");
	printf ("\n fT:\t\tfrecuency of thymine");
	printf ("\n TiTv:\t\ttransition / transversion ratio");
	printf ("\n rAC:\t\tA<>C relative substitution rate");
	printf ("\n rAG:\t\tA<>G relative substitution rate");
	printf ("\n rAT:\t\tA<>T relative substitution rate");
	printf ("\n rCG:\t\tG<>G relative substitution rate");
	printf ("\n rCT:\t\tC<>T relative substitution rate");*/
	printf ("\n (I):\t\taveraged using only +I models.");
	printf ("\n (G):\t\taveraged using only +G models.");
/*	printf ("\n (I+IG):\taveraged using both +I and +I+G models");
	printf ("\n (G+IG):\taveraged using both +G and +I+G models");*/
	printf ("\n (IG):\t\taveraged using only +I+G models.");

}


/************** AverageEstimates **********************/
/*	
	Calculates parameter importance and averaged estimates
*/

void AverageEstimates (char *whichParameter, int numModels, int *modelIndex, int *estimateIndex, double *importance, 
				double *averagedEstimate, double minWeightToAverage)
{
	int	i;

	whichParameter = whichParameter; /* just to avoid warnings */
	
	/* printf ("\n\nDoing model averaging for %s", whichParameter); */
	
	for (i=0; i<numModels; i++)
		{
		if (wIC[modelIndex[i]] < minWeightToAverage)
			continue;
		*importance += wIC[modelIndex[i]];
		*averagedEstimate += wIC[modelIndex[i]] * score[estimateIndex[i]]; 
		
		/*printf ("\nsum  model = %10s\t\tweigth = %10.8f\t\testimate = %10.8f",  
					model[modelIndex[i]].name, wIC[modelIndex[i]], score[estimateIndex[i]]);*/
		}

	/* rescale importance to the total weight of the models included in the confidence interval */
	if (*importance  > 0)
		{
		*averagedEstimate /= *importance;
		*importance /= cumConfidenceWeight; 
		}
	else
		{
		*averagedEstimate = NA;
		}

		/*printf ("\n\t\t\t\t\t\t\t\t\t\t\taveraged estimate = %10.8f",*averagedEstimate);*/

}


/*********************** FindMinWeightToAverage ****************************/
/*  
	Finds the minimum weight we want to average, that is the weight corresponding
	to the last model in the confidence interval specified
*/ 

double FindMinWeightToAverage ()
{
	int		i;
	double	minWeight,cumWeight;
	
	cumWeight = 0;
	for (i=0; i<NUM_MODELS; i++)
		{
		cumWeight += wIC[orderedIC[i]];
		if (cumWeight > averagingConfidenceInterval)
			{
			minWeight = wIC[orderedIC[i]];
			lastModelConfidence = i;
			cumConfidenceWeight = cumWeight;
			break;
			}
		}

	return minWeight;
}





/********************* PrintPaupBlock ************************/
/* Prints a block of PAUP* commands for appending to the data file */
static void WritePaupBlock (int ishLRT)
{
	fT = 1 - (fA + fC + fG);
	
	printf("\n\n\n--\n\nPAUP* Commands Block:");
	printf(" If you want to implement the previous estimates as likelihod settings in PAUP*,");
	printf(" attach the next block of commands after the data in your PAUP file:\n\n");

	if (ishLRT == YES)
		{
		printf("\n[!\nLikelihood settings from best-fit model (%s) selected by hLRT in %s %s on ", modelhLRT, PROGRAM_NAME, VERSION_NUMBER);
		PrintDate (stdout);
		printf ("]");
		}
	else if (useBIC == YES)
		{
		printf("\n[!\nLikelihood settings from best-fit model (%s) selected by BIC in %s %s on ", modelIC, PROGRAM_NAME, VERSION_NUMBER);
		PrintDate (stdout);
		printf ("]");
		}
	else if (useAICc == YES)
		{
		printf("\n[!\nLikelihood settings from best-fit model (%s) selected by AICc in %s %s on ", modelIC, PROGRAM_NAME, VERSION_NUMBER);
		PrintDate (stdout);
		printf ("]");
		}
	else
		{
		printf("\n[!\nLikelihood settings from best-fit model (%s) selected by AIC in %s %s on ", modelIC, PROGRAM_NAME, VERSION_NUMBER);
		PrintDate (stdout);
		printf ("]");
		}
	
	printf("\n\nBEGIN PAUP;");
	printf("\nLset");

	/* Base frequencies */
	printf("  Base=");
	if (fA == fC && fA ==fG && fA == fT)
		printf("equal");
	else
		printf("(%.4f %.4f %.4f)",fA,fC,fG);

	/* Substitution rates */
	if (Ra == Rb && Ra == Rc && Ra == Rd && Ra == Re && Ra == Rf && TiTv == 0)
		printf("  Nst=1");
	else if (TiTv != 0)
		printf("  Nst=2  TRatio=%.4f", TiTv);
	else	
		printf("  Nst=6  Rmat=(%.4f %.4f %.4f %.4f %.4f)", Ra, Rb, Rc, Rd, Re);
	
	/* Rate variation */
	printf("  Rates=");
	if (shape == 0 || shape > 999)
		printf("equal");
	else
		printf("gamma  Shape=%.4f", shape);

	/* Invariable sites */
	printf("  Pinvar=");
	if (pinv == 0)
		printf("0");
	else
		printf("%.4f", pinv);

	printf(";\nEND;");
	printf("\n\n--");

}

/********************* Output ************************/
/* Prints the results of Modeltest  */
static void Output(char *selection, float value)
{
	int i, numK;
	
	fT = 1 - (fA + fG + fC);

	for (i=0; i<NUM_MODELS; i++)
		if (!strcmp (selection, model[i].name))  
			{
			theln = model[i].ln;
			numK = model[i].parameters;
			}

	printf("\n\n Model selected: %s", selection);	
	printf("\n  -lnL  =\t%7.4f", theln); 
	printf("\n   K    =\t%d", numK); 

	if (value > 0)
		{
		if (useBIC == YES)
			printf("\n   BIC  =\t%7.4f\n", value);
		else
			{
			if (useAICc == YES)
				printf("\n   AICc =\t%7.4f\n", value);
			else	
				printf("\n   AIC  =\t%7.4f\n", value);	
			}
		}
		
	printf("\n   Base frequencies: ");		
	if (fA == fC && fA ==fG && fA == fT)
		printf("\n     Equal frequencies");	
	else
	{
		printf("\n     freqA = \t%7.4f", fA);
		printf("\n     freqC = \t%7.4f", fC);
		printf("\n     freqG = \t%7.4f", fG);
		printf("\n     freqT = \t%7.4f", fT);
	}

	printf("\n   Substitution model: ");
	if (Ra == Rb && Ra == Rc && Ra == Rd && Ra == Re && Ra == Rf && TiTv == 0)
		printf("\n     All rates equal");
	else if (TiTv != 0)
		printf("\n    Ti/tv ratio =\t%7.4f", TiTv);	
	else
	{
		printf("\n     Rate matrix");
		printf("\n     R(a) [A-C] = \t%7.4f", Ra); 
		printf("\n     R(b) [A-G] = \t%7.4f", Rb); 
		printf("\n     R(c) [A-T] = \t%7.4f", Rc); 
		printf("\n     R(d) [C-G] = \t%7.4f", Rd); 
		printf("\n     R(e) [C-T] = \t%7.4f", Re); 
		printf("\n     R(f) [G-T] = \t%7.4f", 1.0); 
	}

		printf("\n   Among-site rate variation");
		if (pinv == 0)
			printf("\n     Proportion of invariable sites = 0");
		else
			printf("\n     Proportion of invariable sites (I) = \t%.4f", pinv);	
		
		printf("\n     Variable sites (G)");
		
		if (shape == 0)
			printf("\n      Equal rates for all sites");
		else if (shape > 999) /* shape is infinity */
			printf("\n      Equal rates for all sites (shape parameter = infinity)");	
		else
			printf("\n      Gamma distribution shape parameter = \t%.4f", shape);	
	
}


/********************* SetModel ************************/
/* Sets the parameter estimates for the selected model */
static void SetModel(char *selection)
{
	/* Default parameter estimates for the selected model (JC)*/
	fA = fC = fG = fT = 0.25;
	TiTv = 0;
	Ra = Rb = Rc = Rd = Re = Rf = 1.0;
	shape = 0.0;
	pinv = 0.0; 

	if (!strcmp (selection, "JCI"))   
		pinv = score[4]; 

	else if (!strcmp (selection, "JC+G"))   
		shape = score[7];

	else if (!strcmp (selection, "JC+I+G"))  	     
    {
		pinv = score[10];
		shape = score[11];
	}

	else if (!strcmp (selection, "F81"))       
	{
		fA = score[14];
		fC = score[15];
		fG = score[16];
		fT = score[17];
	}

	else if (!strcmp (selection, "F81+I"))       
	{
		fA = score[20];
		fC = score[21];
		fG = score[22];
		fT = score[23];
		pinv = score[24];
	}
	     
	else if (!strcmp (selection, "F81+G"))   
	{
		fA = score[27];
		fC = score[28];
		fG = score[29];
		fT = score[30];
		shape = score[31];
	}

	else if (!strcmp (selection, "F81+I+G"))       
    {
		fA = score[34];
		fC = score[35];
		fG = score[36];
		fT = score[37];
		pinv = score[38];
		shape = score[39];
	}

	else if (!strcmp (selection, "K80"))       
	{
		TiTv = score[42];
	}

	else if (!strcmp (selection, "K80+I"))       
    {
		TiTv = score[45];
		pinv = score[46];
	}   

	else if (!strcmp (selection, "K80+G"))   
    {
		TiTv = score[49];
		shape = score[50];
    }

	else if (!strcmp (selection, "K80+I+G"))       
	{
		TiTv = score[53];
		pinv = score[54];
		shape = score[55];
	}

	else if (!strcmp (selection, "HKY"))  
    {
		fA = score[58];
		fC = score[59];
		fG = score[60];
		fT = score[61];
		TiTv = score[62];
    }

	else if (!strcmp (selection, "HKY+I"))      
	{
		fA = score[65];
		fC = score[66];
		fG = score[67];
		fT = score[68];
		TiTv = score[69];
		pinv = score[70];
	}
	     
	else if (!strcmp (selection, "HKY+G"))   
	{
		fA = score[73];
		fC = score[74];
		fG = score[75];
		fT = score[76];
		TiTv = score[77];
		shape = score[78];
    }

	else if (!strcmp (selection, "HKY+I+G"))       
    {
		fA = score[81];
		fC = score[82];
		fG = score[83];
		fT = score[84];
		TiTv = score[85];
		pinv = score[86];
		shape = score[87];
    }


	else if (!strcmp (selection, "TrNef"))          
	{
		Ra = score[90];
		Rb = score[91];
		Rc = score[92];
		Rd = score[93];
		Re = score[94];
	}
	 
	else if (!strcmp (selection, "TrNef+I"))  
	{
		Ra = score[97];
		Rb = score[98];
		Rc = score[99];
		Rd = score[100];
		Re = score[101];
		pinv = score[102];

	}
	     
	else if (!strcmp (selection, "TrNef+G"))   
	{
		Ra = score[105];
		Rb = score[106];
		Rc = score[107];
		Rd = score[108];
		Re = score[109];
		shape = score[110];
    }

	else if (!strcmp (selection, "TrNef+I+G"))       
	{
		Ra = score[113];
		Rb = score[114];
		Rc = score[115];
		Rd = score[116];
		Re = score[117];
		pinv = score[118];
		shape = score[119];
    }

	else if (!strcmp (selection, "TrN"))          
	{
		fA = score[122];
		fC = score[123];
		fG = score[124];
		fT = score[125];
		Ra = score[126];
		Rb = score[127];
		Rc = score[128];
		Rd = score[129];
		Re = score[130];
	}
	 
	else if (!strcmp (selection, "TrN+I"))  
	{
		fA = score[133];
		fC = score[134];
		fG = score[135];
		fT = score[136];
		Ra = score[137];
		Rb = score[138];
		Rc = score[139];
		Rd = score[140];
		Re = score[141];
		pinv = score[142];
	}
	     
	else if (!strcmp (selection, "TrN+G"))   
	{
		fA = score[145];
		fC = score[146];
		fG = score[147];
		fT = score[148];
		Ra = score[149];
		Rb = score[150];
		Rc = score[151];
		Rd = score[152];
		Re = score[153];
		shape = score[154];
    }

	else if (!strcmp (selection, "TrN+I+G"))       
	{
		fA = score[157];
		fC = score[158];
		fG = score[159];
		fT = score[160];
		Ra = score[161];
		Rb = score[162];
		Rc = score[163];
		Rd = score[164];
		Re = score[165];
		pinv = score[166];
		shape = score[167];
    }

	else if (!strcmp (selection, "K81"))          
	{
		Ra = score[170];
		Rb = score[171];
		Rc = score[172];
		Rd = score[173];
		Re = score[174];
	}
	 
	else if (!strcmp (selection, "K81+I"))  
	{
		Ra = score[177];
		Rb = score[178];
		Rc = score[179];
		Rd = score[180];
		Re = score[181];
		pinv = score[182];
	}
	     

	else if (!strcmp (selection, "K81+G"))   
	{
		Ra = score[185];
		Rb = score[186];
		Rc = score[187];
		Rd = score[188];
		Re = score[189];
		shape = score[190];
    }

	else if (!strcmp (selection, "K81+I+G"))       
	{
		Ra = score[193];
		Rb = score[194];
		Rc = score[195];
		Rd = score[196];
		Re = score[197];
		pinv = score[198];
		shape = score[199];
    }

	else if (!strcmp (selection, "K81uf"))          
	{
		fA = score[202];
		fC = score[203];
		fG = score[204];
		fT = score[205];
		Ra = score[206];
		Rb = score[207];
		Rc = score[208];
		Rd = score[209];
		Re = score[210];
	}
	 
	else if (!strcmp (selection, "K81uf+I"))  
	{
		fA = score[213];
		fC = score[214];
		fG = score[215];
		fT = score[216];
		Ra = score[217];
		Rb = score[218];
		Rc = score[219];
		Rd = score[220];
		Re = score[221];
		pinv = score[222];
	}
	     
	else if (!strcmp (selection, "K81uf+G"))   
	{
		fA = score[225];
		fC = score[226];
		fG = score[227];
		fT = score[228];
		Ra = score[229];
		Rb = score[230];
		Rc = score[231];
		Rd = score[232];
		Re = score[233];
		shape = score[234];
    }

	else if (!strcmp (selection, "K81uf+I+G"))       
	{
		fA = score[237];
		fC = score[238];
		fG = score[239];
		fT = score[240];
		Ra = score[241];
		Rb = score[242];
		Rc = score[243];
		Rd = score[244];
		Re = score[245];
		pinv = score[246];
		shape = score[247];
    }

	else if (!strcmp (selection, "TIMef"))          
	{
		Ra = score[250];
		Rb = score[251];
		Rc = score[252];
		Rd = score[253];
		Re = score[254];
	}	 
	
	else if (!strcmp (selection, "TIMef+I"))  
	{
		Ra = score[257];
		Rb = score[258];
		Rc = score[259];
		Rd = score[260];
		Re = score[261];
		pinv = score[262];
	}
	     

	else if (!strcmp (selection, "TIMef+G"))   
	{
		Ra = score[265];
		Rb = score[266];
		Rc = score[267];
		Rd = score[268];
		Re = score[269];
		shape = score[270];
    }

	else if (!strcmp (selection, "TIMef+I+G"))       
	{
		Ra = score[273];
		Rb = score[274];
		Rc = score[275];
		Rd = score[276];
		Re = score[277];
		pinv = score[278];
		shape = score[279];
    }

	else if (!strcmp (selection, "TIM"))          
	{
		fA = score[282];
		fC = score[283];
		fG = score[284];
		fT = score[285];
		Ra = score[286];
		Rb = score[287];
		Rc = score[288];
		Rd = score[289];
		Re = score[290];
	}
	 
	else if (!strcmp (selection, "TIM+I"))  
	{
		fA = score[293];
		fC = score[294];
		fG = score[295];
		fT = score[296];
		Ra = score[297];
		Rb = score[298];
		Rc = score[299];
		Rd = score[300];
		Re = score[301];
		pinv = score[302];
	}
	     

	else if (!strcmp (selection, "TIM+G"))   
	{
		fA = score[305];
		fC = score[306];
		fG = score[307];
		fT = score[308];
		Ra = score[309];
		Rb = score[310];
		Rc = score[311];
		Rd = score[312];
		Re = score[313];
		shape = score[314];
    }

	else if (!strcmp (selection, "TIM+I+G"))       
	{
		fA = score[317];
		fC = score[318];
		fG = score[319];
		fT = score[320];
		Ra = score[321];
		Rb = score[322];
		Rc = score[323];
		Rd = score[324];
		Re = score[325];
		pinv = score[326];
		shape = score[327];
    }

	else if (!strcmp (selection, "TVMef"))          
	{
		Ra = score[330];
		Rb = score[331];
		Rc = score[332];
		Rd = score[333];
		Re = score[334];
	}
	 
	else if (!strcmp (selection, "TVMef+I"))  
	{
		Ra = score[337];
		Rb = score[338];
		Rc = score[339];
		Rd = score[340];
		Re = score[341];
		pinv = score[342];
	}
	     
	else if (!strcmp (selection, "TVMef+G"))   
	{
		Ra = score[345];
		Rb = score[346];
		Rc = score[347];
		Rd = score[348];
		Re = score[349];
		shape = score[350];
    }

	else if (!strcmp (selection, "TVMef+I+G"))       
	{
		Ra = score[353];
		Rb = score[354];
		Rc = score[355];
		Rd = score[356];
		Re = score[357];
		pinv = score[358];
		shape = score[359];
    }


	else if (!strcmp (selection, "TVM"))          
	{
		fA = score[362];
		fC = score[363];
		fG = score[364];
		fT = score[365];
		Ra = score[366];
		Rb = score[367];
		Rc = score[368];
		Rd = score[369];
		Re = score[370];
	}
	 
	else if (!strcmp (selection, "TVM+I"))  
	{
		fA = score[373];
		fC = score[374];
		fG = score[375];
		fT = score[376];
		Ra = score[377];
		Rb = score[378];
		Rc = score[379];
		Rd = score[380];
		Re = score[381];
		pinv = score[382];
	}
	     
	else if (!strcmp (selection, "TVM+G"))   
	{
		fA = score[385];
		fC = score[386];
		fG = score[387];
		fT = score[388];
		Ra = score[389];
		Rb = score[390];
		Rc = score[391];
		Rd = score[392];
		Re = score[393];
		shape = score[394];
    }

	else if (!strcmp (selection, "TVM+I+G"))       
	{
		fA = score[397];
		fC = score[398];
		fG = score[399];
		fT = score[400];
		Ra = score[401];
		Rb = score[402];
		Rc = score[403];
		Rd = score[404];
		Re = score[405];
		pinv = score[406];
		shape = score[407];
    }

	else if (!strcmp (selection, "SYM"))          
	{
		Ra = score[410];
		Rb = score[411];
		Rc = score[412];
		Rd = score[413];
		Re = score[414];
	 }
	 
	 else if (!strcmp (selection, "SYM+I"))  
	 {
		Ra = score[417];
		Rb = score[418];
		Rc = score[419];
		Rd = score[420];
		Re = score[421];
		pinv = score[422];
	}
	     
	 else if (!strcmp (selection, "SYM+G"))   
	 {
		Ra = score[425];
		Rb = score[426];
		Rc = score[427];
		Rd = score[428];
		Re = score[429];
		shape = score[430];
    }

	else if (!strcmp (selection, "SYM+I+G"))       
	 {
		Ra = score[433];
		Rb = score[434];
		Rc = score[435];
		Rd = score[436];
		Re = score[437];
		pinv = score[438];
		shape = score[439];
    }


	else if (!strcmp (selection, "GTR"))       
    {
		fA = score[442];
		fC = score[443];
		fG = score[444];
		fT = score[445];
		Ra = score[446];
		Rb = score[447];
		Rc = score[448];
		Rd = score[449];
		Re = score[450];
   }

	else if (!strcmp (selection, "GTR+I"))       
 	{
		fA = score[453];
		fC = score[454];
		fG = score[455];
		fT = score[456];
		Ra = score[457];
		Rb = score[458];
		Rc = score[459];
		Rd = score[460];
		Re = score[461];
		pinv = score[462];
    }
	     
	else if (!strcmp (selection, "GTR+G"))   
	{
		fA = score[465];
		fC = score[466];
		fG = score[467];
		fT = score[468];
		Ra = score[469];
		Rb = score[470];
		Rc = score[471];
		Rd = score[472];
		Re = score[473];
		shape = score[474];
	}

	else if (!strcmp (selection, "GTR+I+G"))       
	{
		fA = score[477];
		fC = score[478];
		fG = score[479];
		fT = score[480];
		Ra = score[481];
		Rb = score[482];
		Rc = score[483];
		Rd = score[484];
		Re = score[485];
		pinv = score[486];
		shape = score[487];
	} 
}

/********************* hLRT ************************/
/* Performs the hypothesis testing and print the results */

/*TestTwoTvRates*/

void hLRT()
{
	if (TestEqualBaseFrequencies (JC, F81) < alpha)
	{
		if (TestTiequalsTv (F81, HKY) < alpha)
		{ 	
			if (TestEqualTiRates (HKY,TrN) < alpha)	
			{
				if (TestEqualTvRates (TrN,TIM) < alpha)	
				{
					if (TestTwoTvRates (TIM,GTR) < alpha)
					{
						if (TestEqualSiteRates (GTR, GTRG) < alpha)
						{
							if (TestInvariableSites (GTRG, GTRIG) < alpha)
								strcpy(modelhLRT,"GTR+I+G");
							else
								strcpy(modelhLRT,"GTR+G");
						}
						else
						{
							if (TestInvariableSites (GTR, GTRI) < alpha)
								strcpy(modelhLRT,"GTR+I");	
							else
								strcpy(modelhLRT,"GTR");							
						}
					}
					else
					{
						if (TestEqualSiteRates (TIM, TIMG) < alpha)
						{			
							if (TestInvariableSites (TIMG, TIMIG) < alpha)
								strcpy(modelhLRT,"TIM+I+G");
								
							else
								strcpy(modelhLRT,"TIM+G");
						}
				
						else
						{
							if (TestInvariableSites (TIM, TIMI) < alpha)
								strcpy(modelhLRT,"TIM+I");
							else
								strcpy(modelhLRT,"TIM");
						}				
					}

				}
				else
				{
					if (TestEqualSiteRates (TrN, TrNG) < alpha)
					{
						if (TestInvariableSites (TrNG, TrNIG) < alpha)
							strcpy(modelhLRT,"TrN+I+G");
						else
							strcpy(modelhLRT,"TrN+G");
					}
			
					else
					{
						if (TestInvariableSites (TrN, TrNI) < alpha)
							strcpy(modelhLRT,"TrN+I");
						else
							strcpy(modelhLRT,"TrN");
					}	
				}
			}
			else		
			{
				if (TestEqualTvRates(HKY, K81uf) < alpha)
				{
					if (TestTwoTvRates(K81uf, TVM) < alpha)
					{
						if (TestEqualSiteRates (TVM, TVMG) < alpha)
						{
							if (TestInvariableSites (TVMG, TVMIG) < alpha)
								strcpy(modelhLRT,"TVM+I+G");
							else
								strcpy(modelhLRT,"TVM+G");
						}
				
						else
						{
							if (TestInvariableSites (TVM, TVMI) < alpha)
								strcpy(modelhLRT,"TVM+I");
							else
								strcpy(modelhLRT,"TVM");			
						}					

					}
					else
					{
						if (TestEqualSiteRates (K81uf, K81ufG) < alpha)
						{
							if (TestInvariableSites (K81ufG, K81ufIG) < alpha)
								strcpy(modelhLRT,"K81uf+I+G");
							else
								strcpy(modelhLRT,"K81uf+G");
						}
				
						else
						{
							if (TestInvariableSites (K81uf, K81ufI) < alpha)
								strcpy(modelhLRT,"K81uf+I");
							else
								strcpy(modelhLRT,"K81uf");			
						}
					}
				}
				else
				{			
					if (TestEqualSiteRates (HKY, HKYG) < alpha)
					{
						if (TestInvariableSites (HKYG, HKYIG) < alpha)
							strcpy(modelhLRT,"HKY+I+G");
						else
							strcpy(modelhLRT,"HKY+G");
					}
			
					else
					{
						if (TestInvariableSites (HKY, HKYI) < alpha)
							strcpy(modelhLRT,"HKY+I");
						else
							strcpy(modelhLRT,"HKY");
					}
				}
			}
		}
		else 
		{
			if (TestEqualSiteRates (F81, F81G) < alpha)
			{
				if (TestInvariableSites (F81G, F81IG) < alpha)
					strcpy(modelhLRT,"F81+I+G");					
				else
					strcpy(modelhLRT,"F81+G");	
			}
			else
			{
				if (TestInvariableSites (F81, F81I) < alpha)
					strcpy(modelhLRT,"F81+I");	
				else
					strcpy(modelhLRT,"F81");	
			}
		}
	}
	else 
	{
		if (TestTiequalsTv (JC, K80) < alpha) 
		{
			if (TestEqualTiRates (K80, TrNef) < alpha)	
			{
				if (TestEqualTvRates (TrNef, TIMef) < alpha)	
				{
					if (TestTwoTvRates (TIMef, SYM) < alpha)
					{
						if (TestEqualSiteRates (SYM, SYMG) < alpha)
						{
							if (TestInvariableSites (SYMG, SYMIG) < alpha)
								strcpy(modelhLRT,"SYM+I+G");	
							else
								strcpy(modelhLRT,"SYM+G");
						}
				
						else
						{
							if (TestInvariableSites (SYM, SYMI) < alpha)
								strcpy(modelhLRT,"SYM+I");
							else				
								strcpy(modelhLRT,"SYM");
						}
					}
					else
					{
						if (TestEqualSiteRates (TIMef, TIMefG) < alpha)
						{
							if (TestInvariableSites (TIMefG, TIMefIG) < alpha)
								strcpy(modelhLRT,"TIMef+I+G");	
							else
								strcpy(modelhLRT,"TIMef+G");
						}
				
						else
						{
							if (TestInvariableSites (TIMef, TIMefI) < alpha)
								strcpy(modelhLRT,"TIMef+I");
							else				
								strcpy(modelhLRT,"TIMef");
						}
					}
				}
				else
				{
					if (TestEqualSiteRates (TrNef, TrNefG) < alpha)
					{
						if (TestInvariableSites (TrNefG, TrNefIG) < alpha)
							strcpy(modelhLRT,"TrNef+I+G");	
						else
							strcpy(modelhLRT,"TrNef+G");
					}
				
					else
					{
						if (TestInvariableSites (TrNef, TrNefI) < alpha)
							strcpy(modelhLRT,"TrNef+I");
						else				
							strcpy(modelhLRT,"TrNef");
					}
				}
				
			}
			else		
			{	
				if (TestEqualTvRates (K80, K81) < alpha)
				{	 	
					if (TestTwoTvRates (K81, TVMef) < alpha)
					{
						if (TestEqualSiteRates (TVMef, TVMefG) < alpha)
						{
							if (TestInvariableSites (TVMefG, TVMefIG) < alpha)
								strcpy(modelhLRT,"TVMef+I+G");
							else
								strcpy(modelhLRT,"TVMef+G");
						}
					
						else
						{
							if (TestInvariableSites (TVMef, TVMefI) < alpha)
								strcpy(modelhLRT,"TVMef+I");
							else
								strcpy(modelhLRT,"TVMef");
						}
					}				
					else
					{
						if (TestEqualSiteRates (K81, K81G) < alpha)
						{
							if (TestInvariableSites (K81G, K81IG) < alpha)
								strcpy(modelhLRT,"K81+I+G");
							else
								strcpy(modelhLRT,"K81+G");
						}
					
						else
						{
							if (TestInvariableSites (K81, K81I) < alpha)
								strcpy(modelhLRT,"K81+I");
							else
								strcpy(modelhLRT,"K81");
						}
					}
				}
				else
				{
				
					if (TestEqualSiteRates (K80, K80G) < alpha)
					{
						if (TestInvariableSites (K80G, K80IG) < alpha)
							strcpy(modelhLRT,"K80+I+G");
						else
							strcpy(modelhLRT,"K80+G");
					}
				
					else
					{
						if (TestInvariableSites (K80, K80I) < alpha)
							strcpy(modelhLRT,"K80+I");
						else
							strcpy(modelhLRT,"K80");
					}
				}
			}
		}
		else 
		{
			if (TestEqualSiteRates (JC, JCG) < alpha)
			{
				if (TestInvariableSites (JCG, JCIG) < alpha)
					strcpy(modelhLRT,"JC+I+G");
				else
					strcpy(modelhLRT,"JC+G");
			}
			
			else
			{
				if (TestInvariableSites (JC, JCI) < alpha)
					strcpy(modelhLRT,"JC+I");
				else
					strcpy(modelhLRT,"JC");
			}
		}
	}
} /* end of method */



/********************* Sioux **********************/
#ifdef macintosh
	static void Sioux()
	{
		/* Don't exit the program after it runs or ask whether to save the
		window when the program exit */
		SIOUXSettings.autocloseonquit = 0;
		SIOUXSettings.asktosaveonclose = 1;

		/* Don't show (0) the status line */
		SIOUXSettings.showstatusline = 0;

		/* Make the window large enough to fit 1 line of text that
		contains 12 characters. */
		SIOUXSettings.columns = 85;
		SIOUXSettings.rows = 50;

		/* Place the window�s top left corner at (5,40). */
		SIOUXSettings.toppixel = 50;
		SIOUXSettings.leftpixel = 350;

		SIOUXSetTitle("\pModeltest v3.7 Console");

		/* Set the font to be 10-point, bold, italic Monaco. */
		SIOUXSettings.fontsize = 9;
		/*SIOUXSettings.fontface = 1 + 2;*/
		SIOUXSettings.fontid = 4;
	}
#endif

/********************* PrintTitle **********************/
static void PrintTitle (FILE *fp)
{
	#ifndef __MWERKS__
		fprintf(fp, "\n\n");
	#endif

	fprintf(fp, "Testing models of evolution - %s ", PROGRAM_NAME);
	fprintf(fp, "%s", VERSION_NUMBER);
	fprintf(fp, "\n(c) Copyright, 1998-2005 David Posada  (dposada@uvigo.es)");
	fprintf(fp, "\nFacultad de Biologia, Universidad de Vigo,\nCampus Universitario, 36310 Vigo, Spain");
	fprintf(fp, "\n_______________________________________________________________\n");
}


/********************* PrintDate ***********************/
static void PrintDate (FILE *fp)
{
	time_t now;
	char *date;
	  
	now=time(NULL);
	date= ctime(&now);
	fprintf(fp, "%s",date);
}


/********************* PrintOS ***********************/
static void PrintOS (FILE *fp)
{
	/* listing preprocessor macros for tried OS
	  	macintosh	:	macintosh GUI with Codewarrior
	  	__APPLE__	: 	macOS X Darwin
	  	__MACH__	:	macOS X Darwin
	  	_WIN32		:	winxp
		unix			:	debian linux
	*/

	#if macintosh
		fprintf(fp, "OS = Macintosh (Sioux console)\n");
	#elif __MACH__
			fprintf(fp, "OS = Macintosh (Darwin terminal)\n");
	#elif  _WIN32
			fprintf(fp, "OS = Windows\n");
	#elif  unix
			fprintf(fp, "OS = Unix-like\n");
	#else 
		fprintf(fp, "OS = unknown\n");
	#endif
}


/********************* CheckExpirationDate ***********************/
static void CheckExpirationDate (int month, int year, int different)
{
	time_t now;
	struct tm *ptime;
	
	 ptime = localtime(&now);
	/* fprintf (stderr, "\nmonth and time = %d %d",ptime->tm_mon,(ptime->tm_year+1900));*/
	/* note that months start at 0 */
	
	if (different == YES &&  (ptime->tm_mon != month-1  ||  ptime->tm_year+1900  != year))
		{}
	else if (different == NO &&  ptime->tm_year+1900 < year)
		{}
	else if (different == NO &&  (ptime->tm_mon < month-1  ||  ptime->tm_year+1900 == year))
		{
		fprintf(stderr, "\nProgram has expired .. now quitting.\n");
		exit(1);
		}
}


/************** CheckNA **********************/
/*	
	If value is NA prints "-"
*/

static char *CheckNA (double value)
{
	char *string;
	
	string = (char*) calloc (10, sizeof (char));       

	if (value == NA)
		return "  -  ";
	else
		{
		sprintf (string, "%8.4f", value);
		return string;
		}	
}


/************** PrintUsage **********************/
static void PrintUsage()
{
	    fprintf(stderr,"\n HELP \n");
	    fprintf(stderr,"\nModeltest is a program for comparing models of evolution using likelihood ");
	    fprintf(stderr,"ratio tests, AIC and BIC. The input are log likelihood scores. ");
	    fprintf(stderr,"You can input raw scores or a Paup scores file resulting from the execution ");
	    fprintf(stderr,"of the provided block of Paup commands (modelblock).\n"); 
	    fprintf(stderr,"\nThe program can also enter in a calculator mode for obtaining ");
	    fprintf(stderr,"the P-value associated with the log likelihood ratio statistic for two given ");
	    fprintf(stderr,"scores\n");
	    fprintf(stderr,"\nJC:    Jukes and Cantor 1969");
   	    fprintf(stderr,"\nK80:   Kimura 2-parameters, Kimura 1980 (also known as K2P)");
	    fprintf(stderr,"\nTrNef: Tamura-Nei 1993 with equal base frequencies");
	    fprintf(stderr,"\nK81:   Kimura 3-parameters, Kimura 1981 (also known as K3ST");
	    fprintf(stderr,"\nTIM:   Transitional model with equal base frequencies");  
	    fprintf(stderr,"\nTVM:   Transversional model with equal base frequencies");  
	    fprintf(stderr,"\nSYM:   Symmetrical model, Zharkikh 1994");  
	    fprintf(stderr,"\nF81:   Felsenstein 1981");
	    fprintf(stderr,"\nHKY:   Hasegawa-Kishino-Yano 1985");
		fprintf(stderr,"\nTrN:   Tamura-Nei 1993");  
	    fprintf(stderr,"\nK81uf: Kimura 3-parameters with unequal base frequencies");    
	    fprintf(stderr,"\nTIM:   Transitional model");  
	    fprintf(stderr,"\nTVM:   Transversional model");  
		fprintf(stderr,"\nGTR:   General time reversible, Rodriguez et al 1990 (also known as REV)");
	    fprintf(stderr,"\n\nI:   invariable sites     G: gamma distribution");
	    fprintf(stderr,"\n\nUsage:   -d : debug level (e.g. -d2)");
	    fprintf(stderr,"\n         -a : alpha level (e.g., -a0.01) (default is a=0.01)");
	    fprintf(stderr,"\n         -n : sample size (e.g., number of characters). Forces the use of AICc (e.g., -c345) (default is to use AIC)");
	    fprintf(stderr,"\n         -t : number of taxa. Forces to include branch lengths as parameters (e.g., -t28) (default is not to count them)");
	    fprintf(stderr,"\n         -w : confidence interval for averaging (e.g., -w0.95) (default is w=1.0)");
	    fprintf(stderr,"\n         -l : LRT calculator mode (e.g., -l)");
	    fprintf(stderr,"\n         -b : Use BIC instead of AIC for all calculations (default is to use AIC)");
	    fprintf(stderr,"\n         -? : help\n");       
 
		fprintf(stderr,"\nUNIX/WIN usage: modeltest3.7 [-d -a -n -t -w -l -i -f -?] < infile > outfile\n\n");
}


/************** Allocate **********************/
static void Allocate()
{ 
	modelhLRT=   (char*) calloc (10, sizeof (char));       
	modelIC=    (char*) calloc (10, sizeof (char)); 
	model=    (ModelSt*) calloc (NUM_MODELS, sizeof (ModelSt)); 
	order=    (ModelSt*) calloc (NUM_MODELS, sizeof (ModelSt)); 
}

/************** Free **********************/
void Free()
{ 
	free (modelhLRT);       
	free (modelIC); 
	free (model); 
	free (order); 
}



