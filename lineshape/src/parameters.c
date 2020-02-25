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

Parameters*
fortran_wrapper(int ligand)
{
  /* this is ugly af but i cannot for the life of me get
   * fortran to pass a string to get_parameters correctly.
   * hence, pick a ligand code in the fortran and then switch */
  Parameters *p;

  /* 
   * lookup table:
   * CLA (Chlorophyll a)  = 0
   * CHL (Chlorophyll b)  = 1
   * KC1 (Chlorophyll c1) = 2
   * KC2 (Chlorophyll c2) = 3
   * A86 (fucoxanthin)    = 4
   * DD6 (diadinoxanthin) = 5
   * LUT (lutein)         = 6
   */

  switch(ligand) {
  case 0:
    p = get_parameters("lineshape/in/CLA.def");
    break;
  case 1:
    p = get_parameters("lineshape/in/CHL.def");
    break;
  case 2:
    p = get_parameters("lineshape/in/KC1.def");
    break;
  case 3:
    p = get_parameters("lineshape/in/KC2.def");
    break;
  case 4:
    p = get_parameters("lineshape/in/A86.def");
    break;
  case 5:
    p = get_parameters("lineshape/in/DD6.def");
    break;
  case 6:
    p = get_parameters("lineshape/in/LUT.def");
    break;
  default:
    fprintf(stdout, "Unknown ligand code %d received from "
        "fortran. Try again.\n", ligand);
    exit(EXIT_FAILURE);
  }
  return p;
}

Parameters*
get_parameters(char *filename)
{
    FILE *fp;
    Parameters *p;
    char line[200], key[200], val[200];
    p->s0 = 0.0; p->s1 = 0.0; p->s2 = 0.0;
    p->g0 = 0.0; p->g1 = 0.0; p->g2 = 0.0;
    p->l0 = 0.0; p->l1 = 0.0; p->l2 = 0.0;
    p->l0 = 0.0; p->w1 = 0.0; p->w2 = 0.0;
    p->ti = 0.0; p->T = 0.0;
    p->cw = NULL;
    fprintf(stdout, "File name read in was %s.\n", filename);

    /* try and do something clever: check via the filename
     * whether we're simulating a chlorophyll or a carotenoid */
    if (strstr(filename, "CLA") != NULL
     || strstr(filename, "CHL") != NULL) {
    	fprintf(stdout, "Ligand name read as %s; using chlorophyll "
    		"spectral density\n", filename);
    	p->ligand = 1;
    } else if (strstr(filename, "CLC") != NULL) {
    	fprintf(stdout, "Ligand name read as %s; using ODO "
    		"spectral density\n", filename);
    	p->ligand = 2;
    } else if (strstr(filename, "A86") != NULL
    	    || strstr(filename, "LUT") != NULL
    	    || strstr(filename, "DD6") != NULL) {
    	fprintf(stdout, "Ligand name read as %s; using carotenoid "
    		"spectral density\n", filename);
	p->ligand = 0;
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
	    p->s0 = atof(val);
	} else if (strcmp(key, "s1") == 0) {
	    p->s1 = atof(val);
	} else if (strcmp(key, "s2") == 0) {
	    p->s2 = atof(val);
	} else if (strcmp(key, "g0") == 0) {
	    p->g0 = atof(val);
	} else if (strcmp(key, "g1") == 0) {
	    p->g1 = atof(val);
	} else if (strcmp(key, "g2") == 0) {
	    p->g2 = atof(val);
	} else if (strcmp(key, "l0") == 0) {
	    p->l0 = atof(val);
	} else if (strcmp(key, "l1") == 0) {
	    p->l1 = atof(val);
	} else if (strcmp(key, "l2") == 0) {
	    p->l2 = atof(val);
	} else if (strcmp(key, "w1") == 0) {
	    p->w1 = atof(val);
	} else if (strcmp(key, "w2") == 0) {
	    p->w2 = atof(val);
	} else if (strcmp(key, "T") == 0) {
	    p->T = atof(val);
	} else if (strcmp(key, "Aw_file") == 0) {
	    strcpy(p->aw_file, val);
	} else if (strcmp(key, "Fw_file") == 0) {
	    strcpy(p->fw_file, val);
	} else if (strcmp(key, "gt_file") == 0) {
	    strcpy(p->gt_file, val);
	} else if (strcmp(key, "lambda_file") == 0) {
	    strcpy(p->lambda_file, val);
	} 

	/* future file names */
	/* else if (strcmp(key, "mag") == 0) { */
	/*     strcpy(config.magFile, val); */
	/* } else if (strcmp(key, "stats") == 0) { */
	/*     strcpy(config.statsFile, val); */
	/* } */ 

	/* future pigment switch */
	/* else if (strcmp(key, "strain") == 0) { */
	/*     /1* don't think you can assign string to enum? *1/ */
	/*     if (strcmp(val, "NONE") == 0) { */
	/*     	config.strain = NONE; */
	/*     } else if (strcmp(val, "S100") == 0) { */
	/*     	config.strain = S100; */
	/*     } else if (strcmp(val, "S110") == 0) { */
	/*     	config.strain = S110; */
	/*     } else if (strcmp(val, "S111") == 0) { */
	/*     	config.strain = S111; */
	/*     } */
	/* } */

    }

    int cl = fclose(fp);
    if (cl != 0) {
    	fprintf(stdout, "Failed to close parameter file %d.\n", cl);
    	exit(EXIT_FAILURE);
    }

    return p;
}
