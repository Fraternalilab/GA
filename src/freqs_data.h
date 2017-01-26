/*===============================================================================
freqs_data.h : code character frequency data
(C) 2007 Jens Kleinjung and Alessandro Pandini
Read the COPYING file for license information.
================================================================================*/

#ifndef FREQS_DATA_H
#define FREQS_DATA_H

typedef struct {
    char setname[64];
    int nFreq;
	int codeRange;
    char codeOrder[64];
    float freq[64];
} ConstantFreqSet;

ConstantFreqSet constant_freq_sets[] = {
    {
        /*___________________________________________________________________________*/
        /** Micheletti et al. oligon alphabet of length 4;
            array size is 28 code frequencies;
            Reference: C. Micheletti, F. Seno and A. Maritan
            Recurrent oligomers in proteins - an optimal scheme reconciling accurate
            and concise backbone representations in automated folding and design studies.
            Proteins: Structure, Function and Genetics 40, 662-674, (2000) */
        "Micheletti2000", 64, 64, "ABCDEFGHIJKLMNOPQRSTUVWXYZab",
        {
        /*o1*/  0.035,  /*o1*/
        /*o2*/  0.035,  /*o2*/
        /*o3*/  0.035,  /*o3*/
        /*o4*/  0.035,  /*o4*/
        /*o5*/  0.035,  /*o5*/
        /*o6*/  0.035,  /*o6*/
        /*o7*/  0.035,  /*o7*/
        /*o8*/  0.035,  /*o8*/
        /*o9*/  0.035,  /*o9*/
        /*o10*/ 0.035, /*o10*/
        /*o11*/ 0.035, /*o11*/
        /*o12*/ 0.035, /*o12*/
        /*o13*/ 0.035, /*o13*/
        /*o14*/ 0.035, /*o14*/
        /*o15*/ 0.035, /*o15*/
        /*o16*/ 0.035, /*o16*/
        /*o17*/ 0.035, /*o17*/
        /*o18*/ 0.035, /*o18*/
        /*o19*/ 0.035, /*o19*/
        /*o20*/ 0.035, /*o20*/
        /*o21*/ 0.035, /*o21*/
        /*o22*/ 0.035, /*o22*/
        /*o23*/ 0.035, /*o23*/
        /*o24*/ 0.035, /*o24*/
        /*o25*/ 0.035, /*o25*/
        /*o26*/ 0.035, /*o26*/
        /*o27*/ 0.035, /*o27*/
        /*o28*/ 0.055  /*o28*/
        }
    },
    {
        /*___________________________________________________________________________*/
        /** Camproux et al. SA27 alphabet of length 4;
            array size is 27 character frequencies;
            Reference: Camproux AC, Gautier R, Tuffery P.
            A hidden markov model derived structural alphabet for proteins.
            J. Mol. Biol. 2004 339:591-605. */
        "Camproux2004", 64, 64, "ABCDEFGHIJKLMNOPQRSTUVWXYZa",
        {
        /*A*/ 0.126, /*A*/
        /*B*/ 0.047, /*B*/
        /*C*/ 0.018, /*C*/
        /*D*/ 0.020, /*D*/
        /*E*/ 0.020, /*E*/
        /*F*/ 0.019, /*F*/
        /*G*/ 0.034, /*G*/
        /*H*/ 0.027, /*H*/
        /*I*/ 0.029, /*I*/
        /*J*/ 0.020, /*J*/
        /*K*/ 0.041, /*K*/
        /*L*/ 0.051, /*L*/
        /*M*/ 0.053, /*M*/
        /*N*/ 0.049, /*N*/
        /*O*/ 0.015, /*O*/
        /*P*/ 0.044, /*P*/
        /*Q*/ 0.041, /*Q*/
        /*R*/ 0.017, /*R*/
        /*S*/ 0.032, /*S*/
        /*T*/ 0.030, /*T*/
        /*U*/ 0.020, /*U*/
        /*V*/ 0.056, /*V*/
        /*W*/ 0.053, /*W*/
        /*X*/ 0.047, /*X*/
        /*Y*/ 0.020, /*Y*/
        /*Z*/ 0.045, /*Z*/
        /*a*/ 0.026  /*a*/
        }           
    }
};

#endif
