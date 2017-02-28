/*==============================================================================
arg.c : parse command line arguments
(C) 2008-2010 A Annibale, L Fernandes, ACC Coolen, J Kleinjung, F Fraternali
Read the COPYING file for license information.
==============================================================================*/

#include "config.h"
#include "arg.h"

/*____________________________________________________________________________*/
/** print version */
static void print_version()
{
    fprintf(stdout, "\n%s %s\n", PACKAGE, VERSION);
}

/*____________________________________________________________________________*/
/** print version */
static void print_header()
{
    fprintf(stdout, "\nNETFEAT : Quantifying Complexity of Cellular Networks\n");
}

/*____________________________________________________________________________*/
/** print license */
static void print_license()
{
    fprintf(stdout, "(C) 2008-2011  Luis Fernandes, Alessia Annibale, Jens Kleinjung,\n"
			"\tAntony C C Coolen, Franca Fraternali\n"
			"netFeat is free software and comes with ABSOLUTELY NO WARRANTY.\n"
			"You are welcome to redistribute it under certain conditions.\n"
			"Read the LICENSE file for distribution details.\n\n");
}

/*____________________________________________________________________________*/
/** print citation */
static void print_citation()
{
    fprintf(stdout, "Mathematical methods:\n"
            "Annibale, A., Coolen, A.C.C., Fernandes, L.P., Fraternali, F., Kleinjung, J.\n"
            "Tailored graph ensembles as proxies or null models for real networks I:\n"
            "tools for quantifying structure.\n"
            "Journal of Physics A: Math. Theor. (2009) 42:485001 (25 pp)\n\n");
    fprintf(stdout, "Application to protein-protein interaction networks:\n"
            "Fernandes, L P, Annibale, A, Kleinjung, J, Coolen, A C C, Fraternali, F\n"
            "Protein networks reveal detection bias and species consistency\n"
			"when analysed by information-theoretic methods.\n"
            "PLoS ONE 5(8):e12083, 2010\n\n");
}

/*____________________________________________________________________________*/
/** set defaults */
static void set_defaults(Arg *arg)
{
	unsigned int n;

	for (n = 0; n < arg->nNetMax; ++ n) {
		arg->net[n].intsList = "";
		arg->net[n].protList = "";
		arg->net[n].matIn = "";
		arg->net[n].matOut = (char *)safe_malloc(128 * sizeof(char));
		sprintf(arg->net[n].matOut, "ints.%c.mat", (char)(n + 48));
		arg->net[n].listOut = (char *)safe_malloc(128 * sizeof(char));
		sprintf(arg->net[n].listOut, "ints.%c.dat", (char)(n + 48));
	}

    arg->matInFlag = 0;
    arg->loopDepth = 3;
    arg->selfInteraction = 0;
    arg->noNetProp = 0;
	arg->compare = 0;
	arg->consistency = 0;
	arg->smooth = 0;
	arg->fillval = 0.;
	arg->mtol = 1;
	arg->ltom = 1;
	arg->silent = 0;
	arg->nNet = 1;
}

/*____________________________________________________________________________*/
/** check input */
static void check_input(Arg *arg)
{
	unsigned int n;

	for (n = 0; n < arg->nNet; ++ n) {
		if (! arg->matInFlag) {
			assert(strcmp(arg->net[n].intsList, "") != 0);
			assert(strcmp(arg->net[n].protList, "") != 0);
		} else {
			assert(strcmp(arg->net[n].matIn, "") != 0);
		}
	}

	assert(arg->matInFlag == 0 || arg->matInFlag == 1);
	assert((arg->loopDepth > 2) && (arg->loopDepth < 10));
	assert(arg->selfInteraction == 0 || arg->selfInteraction == 1);
	assert(arg->noNetProp == 0 || arg->noNetProp == 1);
	assert(arg->compare == 0 || arg->compare == 1);
	assert(arg->consistency == 0 || arg->consistency == 1);
	assert(arg->smooth >= 0);
	assert(arg->fillval >= 0. && arg->fillval < 1.);
	assert(arg->silent == 0 || arg->silent == 1);
	assert(arg->nNet > 0);

	if (arg->compare && (arg->nNet < 2))
		Error("Network comparison requires input of 2 networks");
}

/*____________________________________________________________________________*/
/** parse command line long_options */
int parse_args(int argc, char **argv, Arg *arg)
{
	int c;
	const char usage[] = "Usage: \n\
     * input single interaction list and protein list:\n\
       netfeat [--intsList ...] [--protList ...] [OPTIONS ...]\n\
     * input single interaction matrix:\n\
       netfeat [--matIn ...] [OPTIONS ...]\n\
     * input two interaction lists and protein lists and compare networks:\n\
       netfeat --compare [--intsList ...] [--protList ...] [--intsList1 ...] [--protList1 ...] [OPTIONS ...]\n\
     * input two interaction matrices and compare networks:\n\
       netfeat --compare [--matIn ...] [--matIn1 ...] [OPTIONS ...]\n\
     INPUT\n\
       NETWORK 0\n\
       NOTE: choose either '--intsList and --protList' or '--matInt'\n\
       --intsList\t\t(mode: mandatory, type: char  , default: void)\n\
       --protList\t\t(mode: mandatory, type: char  , default: void)\n\
       --matIn\t\t\t(mode: instead  , type: char  , default: void)\n\
       NETWORK 1 (optional)\n\
       --intsList1\t\t(mode: optional , type: char  , default: void)\n\
       --protList1\t\t(mode: optional , type: char  , default: void)\n\
       --matIn1\t\t\t(mode: instead  , type: char  , default: void)\n\
     MODE\n\
       --loopDepth\t\t(mode: optional , type: int, default: 3)\n\
       --selfInteraction\t\t(mode: optional , type: int, default: 0)\n\
       --noNetProp\t\t(mode: optional , type: no_arg, default: on)\n\
       --compare\t\t(mode: optional , type: no_arg, default: off)\n\
       --smooth\t\t\t(mode: optional , type: int   , default: 0)\n\
       --fillval\t\t(mode: optional , type: double, default: 0. (if value 0, program uses internal value))\n\
     OUTPUT\n\
       --matOut\t\t\t(mode: mandatory, type: char  , default: 'ints.mat')\n\
       --listOut\t\t(mode: mandatory, type: char  , default: 'ints.dat')\n\
       --silent\n\
     INFO\n\
       --version\n\
       --cite\n\
       --help\n";

	if (argc < 2) {
		print_license();
		fprintf(stderr, "%s\n", usage);
		exit(1);
	}

	arg->net = safe_malloc(arg->nNetMax * sizeof(Net));
    set_defaults(arg);

	/*____________________________________________________________________________*/
    /** long option definition */
    static struct option long_options[] =
    {
		{"intsList", required_argument, 0, 1},
		{"protList", required_argument, 0, 2},
		{"matIn", required_argument, 0, 3},
		{"intsList1", required_argument, 0, 4},
		{"protList1", required_argument, 0, 5},
		{"matIn1", required_argument, 0, 6},
		{"matOut", required_argument, 0, 7},
		{"listOut", required_argument, 0, 8},
		{"loopDepth", required_argument, 0, 12},
		{"selfInteraction", no_argument, 0, 13},
		{"noNetProp", no_argument, 0, 14},
		{"compare", no_argument, 0, 15},
		{"consistency", no_argument, 0, 16},
		{"smooth", required_argument, 0, 17},
		{"fillval", required_argument, 0, 18},
		{"silent", no_argument, 0, 19},
		{"cite", no_argument, 0, 20},
		{"version", no_argument, 0, 21},
		{"help", no_argument, 0, 22},
		{0, 0, 0, 0}
    };

	/*____________________________________________________________________________*/
    /** assign parameters to long options */
    while ((c = getopt_long(argc, argv, "1:2:3:4:5:6:7:8:12: 13 14 15 16 17:18: 19 20 21 22 ", long_options, NULL)) != -1) {
        switch(c) {
            case 1:
                arg->net[0].intsList = optarg;
                break;
            case 2:
                arg->net[0].protList = optarg;
                break;
            case 3:
                arg->net[0].matIn = optarg;
				arg->matInFlag = 1;
                break;
            case 4:
                arg->net[1].intsList = optarg;
				arg->nNet = 2;
                break;
            case 5:
                arg->net[1].protList = optarg;
				arg->nNet = 2;
                break;
            case 6:
                arg->net[1].matIn = optarg;
				arg->nNet = 2;
                break;
            case 7:
                strcpy(arg->net[0].matOut, optarg);
                arg->mtol = 1;
                break;
            case 8:
                strcpy(arg->net[0].listOut, optarg); 
                arg->ltom = 1;
                break;
			case 12:
                arg->loopDepth = atoi(optarg);
                break;
            case 13:
                arg->selfInteraction = 1;
                break;
            case 14:
                arg->noNetProp = 1;
                break;
            case 15:
                arg->compare = 1;
				arg->nNet = 2;
                break;
            case 16:
                arg->consistency = 1;
                break;
            case 17:
                arg->smooth = atoi(optarg);
                break;
            case 18:
                arg->fillval = (double)atof(optarg);
                break;
            case 19:
				arg->silent = 1;
				break;
            case 20:
				print_citation();
                exit(1);
            case 21:
				print_version();
                exit(1);
            case 22:
				print_header();
				print_license();
                fprintf(stderr, "%s\n", usage);
                exit(1);
            default:
				print_header();
				print_license();
                fprintf(stderr, "%s\n", usage);
                exit(1);
        }
    }

	check_input(arg);

	/*____________________________________________________________________________*/
	/** print run modes */
	if (! arg->silent) {	
		fprintf(stdout, "Run mode\n");

		fprintf(stdout, "\tUsing loop depth of %d\n", arg->loopDepth);

		if (! arg->selfInteraction)
			fprintf(stdout, "\tEXcluding self-interactions\n");
		else
			fprintf(stdout, "\tINcluding self-interactions\n");

		if (arg->noNetProp)
			fprintf(stdout, "\tNO network property computation\n");
		else
			fprintf(stdout, "\tNetwork property computation\n");

		if (arg->compare)
			fprintf(stdout, "\tProcessing %d networks\n", arg->nNet);

		if (arg->consistency)
			fprintf(stdout, "\tPermforming consistency checks\n");

		if (arg->smooth)
			fprintf(stdout, "\tNumber of smoothening levels is %d\n",
				arg->smooth);

		if (arg->fillval)
			fprintf(stdout, "\tReplace '0' values with '%g'\n", arg->fillval);
	}

    return 0;
}

