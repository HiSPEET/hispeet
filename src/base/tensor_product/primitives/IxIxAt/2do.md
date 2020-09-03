## Anpassung der subop_1 Routinen

### Ausgangsdateien

      ../../tpo_aaa/subop_1__llll.f
      ../../tpo_aaa/subop_1__llu2l_r.f
      ../../tpo_aaa/subop_1__llb4u4_r.f

### Aktionen

* Header anpassen

* Intent von `v` zu `inout` ändern

* Dimension von `u` und `v` setzen

    - zweite auf `nb`
    - dritte auf  `nc`
    
* Deklaration von Input-Argumenten einfügen

      integer,   intent(in)    :: nb                   !< 2nd dimension of u,v
      integer,   intent(in)    :: nc                   !< 3rd dimension of u,v
      real(RNP), intent(in)    :: alpha                !< factor α
      real(RNP), intent(in)    :: beta                 !< factor β

* Ersetzungen

      __NA1__  →  _NA2_
      __NA2__  →  _NA1_
      At       →  A
      z1       →  v
      subroutine SubOp_1(At, u, z1)
        → subroutine PROC(IxIxAt__,_NA1_,_NA2_)(A, nb, nc, alpha, beta, u, v)
      end subroutine SubOp_1
        →  end subroutine PROC(IxIxAt__,_NA1_,_NA2_)

* Zuweisung(en) anpassen
  
      z1(i,j,k) = tmp  →  v(i,j,k) = alpha * tmp + beta * v(i,j,k)
  
  o.ä.

### Testprogram

0. Zu testende Routine durch Anpassen der `#include` Direktive in `vec_test.F` einfügen

1. Preprozessor-Ausgabe prüfen

        ifort -free -E -P vec_test.fpp > vec_test.f

2. Korrektheit prüfen 

        ifort -free -O2 -check all -o vec_test vec_test.fpp
        ./vec_test

3. Vektorisierung anschauen

        ifort -free -O2 -o vec_test -qopt-report=2 -qopt-report-phase=vec vec_test.fpp
   * Optimierungsreport `vec_test.optrpt` im Editor öffnen
   * Mit `vec_test.fpp` sowie `vec_test.f` vergleichen
   * Hinweisen nachgehen
   
4. Leistungstest

        ifort -free -O2 -o vec_test vec_test.fpp
        ./vec_test

5. Änderungen vornehmen und mit 1.  fortsetzen

