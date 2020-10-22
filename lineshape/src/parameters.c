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
	    fprintf(stdout, "tau = %lu\n", p.ns);
	} else if (strcmp(key, "chl_ansatz") == 0) {
	    p.chl_ansatz = atoi(val);
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
get_parameters(char *filename, ansatz chl_ansatz)
{
    FILE *fp;
    Parameters p;
    char s[4]; /* null char */
    char pref[3] = "gsw";
    char line[200], key[200], val[200];
    p.s0 = 0.0; p.s1 = 0.0; p.s2 = 0.0;
    p.g0 = 0.0; p.g1 = 0.0; p.g2 = 0.0;
    p.l0 = 0.0; p.l1 = 0.0; p.l2 = 0.0;
    p.l0 = 0.0; p.w1 = 0.0; p.w2 = 0.0;
    p.ti = 0.0; p.T = 0.0;  p.nu = 0.0;
    p.offset = 0.0; p.cw = NULL; p.cn = NULL;
    memset(p.gsw, 0., 3*48*sizeof(double));

    /* check filename to get ligand and spectral density ansatz.
     * assign to p.ligand here because functions.h isn't included*/
    if (strstr(filename, "CLA") != NULL
     || strstr(filename, "CHL") != NULL) {
    	fprintf(stdout, "Ligand name read as %s; assuming chlorophyll\n",
    	        filename);
    	p.ans = chl_ansatz;
    } else if (strstr(filename, "KC1") != NULL
    	    || strstr(filename, "KC2") != NULL) {
      /* NB: this is leftover from old FCP bits where I was using the
       * Renger form for chlorophyll a and b, but OBOs for chlorophyll
       * c1 and c2; right now we don't have to worry about it, but it
       * might be relevant later */
    	fprintf(stdout, "Ligand name read as %s; using OBO "
    		"spectral density\n", filename);
    	p.ans = OBO;
    } else if (strstr(filename, "A86") != NULL
    	    || strstr(filename, "LUT") != NULL
    	    || strstr(filename, "NEX") != NULL
    	    || strstr(filename, "XAT") != NULL
    	    || strstr(filename, "DD6") != NULL) {
    	fprintf(stdout, "Ligand name read as %s; using carotenoid "
    		"spectral density\n", filename);
	p.ans = CAR;
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
	} else if (strcmp(key, "nu") == 0) {
	    p.nu = atof(val);
	} else if (strcmp(key, "Aw_file") == 0) {
	    strcpy(p.aw_file, val);
	} else if (strcmp(key, "Fw_file") == 0) {
	    strcpy(p.fw_file, val);
	} else if (strcmp(key, "gt_file") == 0) {
	    strcpy(p.gt_file, val);
	} else if (strcmp(key, "offset_file") == 0) {
	    strcpy(p.offset_file, val);
	} else if (strcmp(key, "lambda_file") == 0) {
	    strcpy(p.lambda_file, val);
	} 
	/* if we're using the big Novoderezhkin/Mancal ansatz
	 * for the spectral density we have 3 * 48 = 144 extra
	 * parameters; generate the strings for these on the fly */
        for (unsigned int j = 0; j < 3; j++) {
          snprintf(s, 3, "%1s", pref + j);
          for (unsigned int i = 0; i < 48; i++) {
            snprintf(s + 1, 3, "%02d", i + 1);
            if (strcmp(key, s) == 0) {
              p.gsw[j][i] = atof(val);
            }
          }
        }
    }

    int cl = fclose(fp);
    if (cl != 0) {
    	fprintf(stdout, "Failed to close parameter file %d.\n", cl);
    	exit(EXIT_FAILURE);
    }

    return p;
}
