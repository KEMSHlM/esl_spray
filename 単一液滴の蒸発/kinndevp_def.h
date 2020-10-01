/* Characteristics for the 62-step kinetic model for the Fuel n-Heptane */

#define spmax	3

#define nFUEL 	0
#define nO2 1
#define nN2	2


typedef enum SpeciesLabel {

        /* Computed species s.. */
        /* Steady-state species ss.. */
 
		sFUEL	=	0,
		sO2		=	1,
		sN2		=	2,
        sEnd
} SpeciesLabel;
