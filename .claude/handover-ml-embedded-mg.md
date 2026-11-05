# Übergabe: eingebetteter Mehrgitter-Drucklöser im Multilevel-BDF2

Branch: `ml-ins-coupled-embedded-mg`. Alle Änderungen liegen uncommittet im
Arbeitsverzeichnis. **Nichts davon wurde bisher übersetzt oder ausgeführt.**

## Arbeitsweise

- Jörg geht offene Punkte Schritt für Schritt durch: Vor größeren Änderungen
  fragen, nicht selbstständig committen.
- Sprache mit Jörg: Deutsch. Code, Kommentare, Commit-Messages: Englisch.
- Geänderte Fortran-Dateien folgen dem Skill `hispeet-fortran-style`.
- Umgebung: `module load hispeet/gnu` (MPI aus `/opt/homebrew`, `ParMETIS_ROOT`).
- Build-Ablauf noch nicht verifiziert: `build/` ist vorhanden und mit
  `/opt/homebrew/bin/mpifort` konfiguriert, Ziel `ml_ins_solver_3d`.
  Vor dem ersten Lauf den tatsächlichen Ablauf prüfen.

## Ziel

Fall `coupled = T`: `StokesCycle` ruft auf jedem Level `ins_op(l) % StokesSolver`.
Der Druck soll dort mit `'MG'` oder `'MGCG'` gelöst werden können. Dafür erzeugt
und besitzt `ML_INS_Operator_3D` einen gemeinsamen Löser `ml_solver_p`
(Option 2 der Diskussion).

## Entscheidungen

- Konstruktor `ml_ins = ML_INS_Operator_3D(...)` bleibt (funktioniert laut Jörg
  auch mit ifx und nvfortran). Kein Umbau auf `Init`, kein `target` an `this`.
- `ml_solver_p` und `ml_op_p` bleiben `allocatable`. Eine normale Komponente läge
  im Funktionsergebnis, die Zeiger in `ins_op(l)` wären danach ungültig.
- Lokal verfeinerte Gitter sind Zielanwendung: keine Rückkehr bei `part < 0`.
  Leere Ränge auf `l_top_` können auf gröberen Levels aktiv sein.
- OpenMP-Races (fehlende Barrieren, `bv_l` als `save`-Zeiger) bleiben hier
  unbearbeitet, laut Jörg wohl im OpenMP-Zweig gelöst.
- `n_parts == 0` auf `l_top_` ist ausgeschlossen.

## Geänderte Dateien

- `src/ins_flow/multilevel/ml__ins__operator__3d.f`: Option
  `ml_ins_opt%ml_solver_p` mit Bcast, Komponente `ml_solver_p`, nur bei
  `'MG'`/`'MGCG'` alloziert und an alle `ins_op(l)` übergeben.
- `src/ins_flow/operator/ins__operator__3d__mp_pressuresolver.f`: fehlendes
  `i_max` im Mixed-Order-Aufruf von `ML_PressureSolver` ergänzt.
- `src/ins_flow/operator/ins__operator__3d.f`: Fehlerabbruch bei
  `'MG'`/`'MGCG'` ohne `ml_solver_p`.
- `src/base/multilevel_3d/elliptic/ml__dg__elliptic_solver__3d__mp_cs_mg_solver.f`
  und `..._mp_cs_mgcg_solver.f`:
  - `active = part >= 0` auf `l_top_`; Skalarprodukte und `CalibrateArray` nur
    aktiv über `comm_top` (`comm_parts`).
  - `converged` und `r_2` per Bcast über `comm_all` (`comm_world`) mit
    `root = proc_part(0)`, nicht 0: Weltrang 0 kann auf `l_top_` leer sein.
  - MG: `r % Init` nur bis `l_top_`.
  - MGCG: `singular = F` auf leeren Rängen, `ni = min(i, i_max_)`,
    `r_new = r_old` nach dem Startresiduum.

## Erwartete Verhaltensänderung

Mixed-Order-Läufe mit MG nutzen jetzt `i_max_p` statt `ml_solver_p%i_max`, auch
Single-Level. Ohne `r_red`/`r_max` im Löser gibt es keine Konvergenzprüfung,
es laufen also bis zu `i_max_p` Zyklen (Default 1000).

## Nächste Schritte

1. Übersetzen. Möglicher Stolperstein: nicht alloziertes `this % ml_solver_p`
   als fehlendes optionales Argument (F2008), bei nvfortran ungeprüft.
2. Referenzlauf `develop/ins_flow/multilevel/variable_viscosity/test-coupled-4g.prm`
   mit dem bisherigen Drucklöser.
3. Kopie der Parameterdatei mit:

   ```fortran
   ml_ins_opt%ins_op%pressure_solver = 'MGCG'
   ml_ins_opt%ins_op%i_max_p         = 5
   ml_ins_opt%ml_solver_p%i_crs      = 1000
   ml_ins_opt%ml_solver_p%r_red      = 1e-08
   ml_ins_opt%ml_solver_p%r_max      = 1e-12
   ```

   Der Fall ist vollperiodisch, daher `MGCG` (`MG` kalibriert nicht).
   `i_crs` ist wichtig: Default 1 ist auf Level 1 zu wenig.
4. Mit mehr Rängen als `n_parts_root = 8` laufen lassen, damit es inaktive
   Ränge gibt. Deadlocks prüfen, Ergebnis mit dem Referenzlauf vergleichen.
5. Fall mit lokaler Verfeinerung testen.

## Offene Fragen

- Eingefrorene Elemente: `CS_MG_Solver_X` addiert die Korrektur auf allen
  Elementen (`MergeArrays`), `FAS_MG_Cycle` behandelt sie gesondert. Nicht geprüft.
- `penalty`, `schwarz`, `interior_bc` existieren in `ins_op` und `ml_solver_p`
  getrennt. Ob sie übereinstimmen müssen, ist nicht geprüft.
- Kommentare beider CS-Löser nennen noch „global refinement“.
- Jörg hatte noch weitere Punkte angekündigt, die er durchgehen möchte.
