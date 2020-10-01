/* Characteristics for the 62-step kinetic model for the Fuel n-Heptane */

#define spmax	3
#define ssmax	3
#define remax	1
#define Mmmax	1

#define nFUEL 	0
#define nN2	2
#define nO2	1

#define mM1	0


typedef enum SpeciesLabel {

        /* Computed species s.. */
        /* Steady-state species ss.. */
        sNXC7H16 = 0,
        sO2 = 1,
        sN2 = 2,
        sEnd
} SpeciesLabel;