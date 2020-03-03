## Funding

1. DFG Antrag: MG für variable Koeffizienten und bewegte Gitter

   * Vorarbeiten
       - lineare p-MG für kondensierte CG-SEM
       - lineare p-MG für hybride DG-SEM
       - p-MG und Strömungslöser für variable Koeffizienten
       - effiziente Operatoren, hybride Parallelisierung
       - HiSPEET, künftig public domain

   * Ziele
       - lineare Löser für variable Koeffizienten?
       - Anwendung auf Phasenfeldmethoden
       - Ursachen für Probleme bei kurzwellige Koeffizienten finden
       - Lösung durch alternative Techniken: FAS, pCG/Schwarz-Glätter, ... ?
       - bessere Skalierbarkeit durch Kombination p/h-MG

## Theory

1. SDC for advection-diffusion equation
   *

## Code

1. Linear IPH solvers for elliptic equations with constant coefficients
   * residual evaluation
   * Schwarz method on vertex-centered subdomains
   * Schwarz/MG-CG
   * integration into flow solver

1. More convenient and safer TPOs
   * use assumed shape instead of explicit one
   * argument checking switched on in DEBUG mode
