#include "parameters.h"

Protocol
get_protocol(char *filename)
{
    FILE *fp;
    Protocol p;
    char line[200], key[200], val[200];
    p.T = 0.0; p.ns = 0;

    fp = fopen(filename, "r");
    if (fp == NULL) {
	fprintf(stdout, "Invalid parameter filename: %s."
		" Try again.\n", filename);
	exit(EXIT_FAILURE);
    }
    while (fgets(line, 199, fp) != NULL) {
    	sscanf(line, " %s = %s", key, val);
	if (strcmp(key, "#") == 0) {
	    continue;
	} else if (strcmp(key, "T") == 0) {
	    p.T = atof(val);
	} else if (strcmp(key, "ns") == 0) {
	    p.ns = atoi(val);
	}
    }

    int cl = fclose(fp);
    if (cl != 0) {
    	fprintf(stdout, "Failed to close protocol file %d.\n", cl);
    	exit(EXIT_FAILURE);
    }

    return p;
}

Parameters
get_parameters(char *filename)
{
    FILE *fp;
    Parameters p;
    char line[200], key[200], val[200];
    p.s0 = 0.0; p.s1 = 0.0; p.s2 = 0.0;
    p.g0 = 0.0; p.g1 = 0.0; p.g2 = 0.0;
    p.l0 = 0.0; p.l1 = 0.0; p.l2 = 0.0;
    p.l0 = 0.0; p.w1 = 0.0; p.w2 = 0.0;
    p.ti = 0.0; p.T = 0.0;
    p.cw = NULL;
    fprintf(stdout, "File name read in was %s.\n", filename);

    /* check filename to get ligand and spectral density ansatz.
     * assign to p.ligand here because functions.h isn't included*/
    if (strstr(filename, "CLA") != NULL
     || strstr(filename, "CHL") != NULL) {
    	fprintf(stdout, "Ligand name read as %s; using chlorophyll "
    		"spectral density\n", filename);
    	p.ligand = 1;
    } else if (strstr(filename, "CLC") != NULL) {
    	fprintf(stdout, "Ligand name read as %s; using ODO "
    		"spectral density\n", filename);
    	p.ligand = 2;
    } else if (strstr(filename, "A86") != NULL
    	    || strstr(filename, "LUT") != NULL
    	    || strstr(filename, "DD6") != NULL) {
    	fprintf(stdout, "Ligand name read as %s; using carotenoid "
    		"spectral density\n", filename);
	p.ligand = 0;
    } else {
    	fprintf(stdout, "Not sure what kind of spectral density"
    		" function to use. File name read in was %s. Try again.\n",
    		filename);
    	exit(EXIT_FAILURE);
    }

    fp = fopen(filename, "r");
    if (fp == NULL) {
	fprintf(stdout, "Invalid parameter filename: %s."
		" Try again.\n", filename);
	exit(EXIT_FAILURE);
    }
    while (fgets(line, 199, fp) != NULL) {
    	sscanf(line, " %s = %s", key, val);
	if (strcmp(key, "#") == 0) {
	    continue;
	} else if (strcmp(key, "s0") == 0) {
	    p.s0 = atof(val);
	} else if (strcmp(key, "s1") == 0) {
	    p.s1 = atof(val);
	} else if (strcmp(key, "s2") == 0) {
	    p.s2 = atof(val);
	} else if (strcmp(key, "g0") == 0) {
	    p.g0 = atof(val);
	} else if (strcmp(key, "g1") == 0) {
	    p.g1 = atof(val);
	} else if (strcmp(key, "g2") == 0) {
	    p.g2 = atof(val);
	} else if (strcmp(key, "l0") == 0) {
	    p.l0 = atof(val);
	} else if (strcmp(key, "l1") == 0) {
	    p.l1 = atof(val);
	} else if (strcmp(key, "l2") == 0) {
	    p.l2 = atof(val);
	} else if (strcmp(key, "w1") == 0) {
	    p.w1 = atof(val);
	} else if (strcmp(key, "w2") == 0) {
	    p.w2 = atof(val);
	} else if (strcmp(key, "T") == 0) {
	    p.T = atof(val);
	} else if (strcmp(key, "Aw_file") == 0) {
	    strcpy(p.aw_file, val);
	} else if (strcmp(key, "Fw_file") == 0) {
	    strcpy(p.fw_file, val);
	} else if (strcmp(key, "gt_file") == 0) {
	    strcpy(p.gt_file, val);
	} else if (strcmp(key, "lambda_file") == 0) {
	    strcpy(p.lambda_file, val);
	} 
    }

    int cl = fclose(fp);
    if (cl != 0) {
    	fprintf(stdout, "Failed to close parameter file %d.\n", cl);
    	exit(EXIT_FAILURE);
    }

    return p;
}
