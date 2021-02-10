# Benceno Production Reactor
 Diseño de un reactor PFR para la obtención de benceno a partir de tolueno
 
 El proceso seleccionado para desarrollar el proyecto es la hidrodesalquilación térmica de tolueno (HDA).
Las reacciones implicadas, en fase gaseosa, son:

Reacción principal:                      C7H8 + H2 → C6H6 + CH4

Reacción secundaria:                      2 C6H6 ⇌ C12H10 + H2

En la reacción principal el tolueno e hidrógeno reaccionan, produciéndose la eliminación del grupo metilo y generando una molécula de benceno y otra de metano como productos de reacción. La reacción es estable en un rango de temperaturas de 870 – 1030 K y a una presión de 500 psi (34 –35 atm).

 
Para facilitar la comprensión del programa diseñado, se especifican las variables que en él intervienen:

	Variable independiente: L, longitud del reactor.
	Variables dependientes: y, un total de 7 que engloban: 
 	  n_j  (kmol/s): Un total de 5 elementos, cada uno referido al flujo molar de cada elemento que interviene en el proceso.
 	  T (K): Temperatura.
 	  P (Pa): Presión.
  
Se debe definir una corriente de alimento al reactor como inicio de la simulación; dicha corriente tiene las siguientes características:

 	Temperatura de entrada: 922 K.
 	Presión de entrada: 3,447·10^6 Pa.
	Se alimenta al reactor un total de 356.287 tm/año.

